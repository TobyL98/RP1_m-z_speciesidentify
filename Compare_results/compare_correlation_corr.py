#####################
## Compare_correlation.py
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
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

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


# function does the actual cross correlation
# using scipy.signal.correlate
# then takes out the correlation at 0 lag
# and the mean of the correlation at lags
# and returns the differnece between the two
def cross_corr(act_peaks, theor_peaks):
    corrs = signal.correlate(theor_peaks, act_peaks)
    lags = signal.correlation_lags(len(theor_peaks), len(act_peaks))
    corr0 = corrs[np.where(lags == 0)]
    corr_mean = np.mean(corrs)
    corr_diff = corr0 - corr_mean
    return(corr_diff)


# function does the comparison between one set of theoretical peptides
# and the PMF using cross correlation
# similar to the orginal SEQUEST paper
# splits the peaks list up into 9 smaller sets
# returns the means of the corr_dif from those 9 sets
def compare_corr(theor_peaks, act_peaks):
    
    MZ_theor_peaks = theor_peaks["mass1"]
    MZ_act_peaks = act_peaks["MZ"]

    # peak list of how we split the data for comparison
    peak_list = [(800, 1100), (1100, 1400), (1400, 1700), (1700, 2000),
                 (2000, 2300), (2300, 2600), (2600, 2900), (2900, 3200),
                 (3200, 3500)]

    corr_diff_list = []
    # goes through the peak list filters and filter both datasets
    # then performs the cross correlation function
    for x in peak_list:
        filter_act_peaks = MZ_act_peaks[MZ_act_peaks.between(x[0], x[1], "left")]

        filter_theor_peaks = MZ_theor_peaks[MZ_theor_peaks.between(x[0], x[1], "left")]
        filter_theor_peaks = filter_theor_peaks.sort_values()

        corr_diff = cross_corr(filter_theor_peaks, filter_act_peaks)
        corr_diff_list.append(corr_diff)

    # finds the mean for corr_diff
    corr_diff_array = np.concatenate(corr_diff_list)
    corr_diff_mean = np.mean(corr_diff_array)
    return(corr_diff_mean)
     

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
match_results_df = pd.concat(results_list)
match_results_df = match_results_df.sort_values(by =["Match"], ascending = False)
#print(final_results_df.head())


###### STARTING CROSS CORRELATION

# picking out top 20 species for cross correlation
# as cross correlation computationall expensive
species_file_list = []
for species in match_results_df["SPECIES"].head(10):
    species = species.replace(" ", "_")
    species_file = "{0}_filtered.csv".format(species)
    species_file_list.append(species_file)

# filter the CSV files by the top 20 species
# calculates the cross correlation score
# score is value in dictionary
# species name is the key
cross_corr_dict = {}
for csv in csv_files:
    csv_name = csv.split("\\")[-1]
    if csv_name in species_file_list:
        theor_peaks_df = pd.read_csv(csv, sep = ",", dtype = dtype, 
                                 usecols = usecols)
        
        # run the compare correlation function
        # returns the corr_mean value
        # adds a value to the dictionary
        corr_mean = compare_corr(theor_peaks_df, act_peaks_df)
        csv_names = csv_name.split("_")
        species_name = csv_names[0] + " " + csv_names[1]
        cross_corr_dict[species_name] = corr_mean

# converts that dictionary into a dataframe
cross_corr_df = pd.DataFrame.from_dict(cross_corr_dict, orient = "index", columns = ["Corr_Score"])

# merges the initial match dataframe with the cross correlation dataframe
final_df = match_results_df.merge(cross_corr_df, how="inner", left_on= "SPECIES", right_index= True)

# normalises both the match_score and Corr_Score to 1
# by dividing by the maximum value
final_df["match_score"] = final_df["Match"] / final_df["Match"].max()
final_df["Corr_Score"] = final_df["Corr_Score"] / final_df["Corr_Score"].max()

# calculates a final ID_score which is the final ranking
# includes 0.3 weight for the Corr_Score and 0.7 for the match_score
#final_df["ID_score"] = (final_df["Corr_Score"] * 0.3) + (final_df["match_score"] * 0.7)
#final_df = final_df.drop(["Corr_Score", "match_score"], axis = 1)
#final_df = final_df.sort_values(by= ["ID_score"], ascending= False)

# outputs the final_df to the specified file path by user
output_path = args.output
final_df.to_csv(output_path)
print(final_df.head(10))



