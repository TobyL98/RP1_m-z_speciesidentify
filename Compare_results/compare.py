#####################
## Compare.py
####################

# function gets the theoretical peptides generated from the sequences
# and filtered by the LCMSMS data
# and compares how many match an actual PMF within a certain tolerance
# tolerance is usually 0.2
# returns the resuls with the number of matches for each species
# ordered by highest number of matches
# highest number of matches should be the species most closely related
# to species in the database

import pandas as pd
import glob
import argparse

################
# FUNCTIONS
################

# function does the comparison between one set of theoretical peptides
# and the PMF within a certain allowance
# theor_peaks are the theoretical peaks
# act_peaks are the actual peaks from PMF

def compare(theor_peaks, act_peaks):
    
    #creates columns that are the same in each
    #can then be merged
    act_peaks["join"] = 1
    theor_peaks["join"] = 1

    # generating allowance for comparison
    # with theoretical peaks
    act_peaks["MZ_plus"] = act_peaks["MZ"] + 0.2
    act_peaks["MZ_minus"] = act_peaks["MZ"] - 0.2
    
    # merges the datasets so we have all the data
    # then filters where there is a MATCH within the tolerance
    # number of rows left will be number of matches
    merged_df = act_peaks.merge(theor_peaks, how = "outer", on = ["join"])
    merged_df = merged_df.loc[(merged_df["mass1"] >= merged_df["MZ_minus"]) &
                              (merged_df["mass1"] <= merged_df["MZ_plus"])]
    
    # remove duplicates if MZ value has matched more than once
    merged_df = merged_df.drop_duplicates(subset = ["MZ"])

    # assigns number of matches to count
    match_count = merged_df.shape[0]
                
    # turns match results to a df
    result_df = pd.DataFrame([match_count], columns = ["Match"])

    # combines with the taxon information to identify which species it is
    taxon_df = theor_peaks.loc[[0] ,["CLASS", "SUBCLASS", "INFRACLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"]]
    final_df = pd.concat([taxon_df , result_df], axis = 1)
    return(final_df)
     

################
# Main code
################
# set up argparse
parser = argparse.ArgumentParser()
# add arguments
# adds the folder where the input theoretical peak spectrums are
parser.add_argument("-it", "--inputTheor",
                    help = "the folder that contains the theoretical peptides csv files to compare against PMF",
                    default= "C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/integrate_results")
# adds where the output file should be saved
parser.add_argument("-o", "--output",
                    help= "The file name for the output results file of number of matches",
                    default = "Results/matches.csv")
# adds where the input peptide mass fingerprint
parser.add_argument("-ip", "--inputPMF",
                    help= "The input Peptide mass fingerprint (PMF) from an unknown organism.")
args = parser.parse_args()

# read in txt file of PMF values from data
input_PMF = args.inputPMF
dtype= {"MZ": 'float32', "intensity": 'float32'}
act_peaks_df = pd.read_table(input_PMF, sep = "\t", header = None, names = ["MZ", "intensity"], dtype= dtype)

# get all csv files
file_path = "{0}/*.csv".format(args.inputTheor)
csv_files = glob.glob(file_path)

# read in all csvs of theoretical peptide peaks
dtype = {"mass1": 'float32', "CLASS": 'category', "SUBCLASS": 'category', "ORDER": 'category',
          "FAMILY": 'category', "GENUS": 'category', "SPECIES": 'category'}
usecols = ["mass1", "CLASS", "SUBCLASS", "INFRACLASS", "ORDER", "FAMILY", "GENUS", "SPECIES"]
results_list = []
for csv in csv_files:
    theor_peaks_df = pd.read_csv(csv, sep = ",", dtype = dtype, 
                                 usecols = usecols)

    # run the compare function
    result_df = compare(theor_peaks_df, act_peaks_df)
    # add all results to a list
    results_list.append(result_df)

# put all results in one dataframe
final_results_df = pd.concat(results_list)
final_results_df = final_results_df.sort_values(by =["Match"], ascending = False)
print(final_results_df.head())
output_path = args.output
final_results_df.to_csv(output_path)


