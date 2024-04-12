#####################
## scoring_system_NCBI.py
####################

# does the same but with NCBI data downloaded by user
# function gets the theoretical peptides generated from the sequences
# and filtered by the LCMSMS data
# and compares how many match an actual PMF within a certain tolerance
# tolerance is usually 0.2
# takes the top 20 from that match
# and calculates with scoring system based on MOWSE
# final results are filtered based on MOWSE score

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
    taxon_df = theor_peaks.loc[[0] ,["GENUS", "SPECIES"]]
    final_df = pd.concat([taxon_df , result_df], axis = 1)
    return(final_df)


# does a score similar to the MOWSE score in mascot
# although adjusted slightly as collagen sequences
# tend to be very similar
# creates a weighting system that puts greater weight 
# at matches that have a higher m/z value
# or where you are less likely to get a false positive
def mowse(counts, act_peaks, theor_peaks, range_list):
    #creates columns that are the same in each
    #can then be merged
    act_peaks["join"] = 1
    theor_peaks["join"] = 1

    # generating allowance for comparison
    # with theoretical peaks
    act_peaks["MZ_plus"] = act_peaks["MZ"] + 0.2
    act_peaks["MZ_minus"] = act_peaks["MZ"] - 0.2

    score_list = []
    index = 0
    # loops through rows in count
    for index, row in counts.iterrows():
        # gets the start and end range for m/z value matches
        start_range = row["start"]
        end_range = row["end"]

        # filters both sets of peak values by this range
        filter_act_peaks = act_peaks[act_peaks["MZ"].between
                                                    (start_range, end_range, "left")]
        filter_theor_peaks = theor_peaks[theor_peaks["mass1"].between
                                                    (start_range, end_range, "left")]
        
        # merges the datasets so we have all the data
        # then filters where there is a MATCH within the tolerance
        # number of rows left will be number of matches
        merged_df = filter_act_peaks.merge(filter_theor_peaks, how = "outer", on = ["join"])
        merged_df = merged_df.loc[(merged_df["mass1"] >= merged_df["MZ_minus"]) &
                                  (merged_df["mass1"] <= merged_df["MZ_plus"])]
        
        # assigns number of matches to count
        match_count = merged_df.shape[0]

        # creates the weight by inverting the normalised counts
        # means more weight will be given to ranges with less counts
        if row["Normalised_Count"] != 0:
            weight = 1 / row["Normalised_Count"]
        else:
            weight = 0
        # multiplies the count by the weight
        # weight is the specific weight for that range
        # to get a final result
        match_result = weight * match_count

        score_list.append(match_result)
        index += 1
    # sums the results for each range to get final score
    final_score = sum(score_list)
    final_score = round(final_score, 2)
    return(final_score)

################
# Main code
################
# set up argparse
parser = argparse.ArgumentParser()
# add arguments
# adds the folder where the input theoretical peak spectrums are
parser.add_argument("-it", "--inputTheor",
                    help = "the folder that contains the theoretical peptides csv files to compare against PMF",
                    default= "C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/results_NCBI")
# adds where the output file should be saved
parser.add_argument("-o", "--output",
                    help= "The file name for the output results file of number of matches",
                    default = "Results/NCBI/matches.csv")
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
dtype = {"mass1": 'float32', "GENUS": 'category', "SPECIES": 'category'}
usecols = ["mass1", "GENUS", "SPECIES"]
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

##############################
## MOWSE like scoring system
##############################


# counts the number of masses at a particular range in actual (experimental) peaks
# for species in the top 20
# and finds the average for each range
# allows matches to be weighted by likelihood
# of incorrect matches by 
# counts are then used in mowse function
range_list = [800, 1100, 1400, 1700, 2000,  2300,
              2600,  2900, 3200,  3500]

act_peaks = act_peaks_df["MZ"]

count_list = []
index = 0
for num in range_list:
    start_range = range_list[index]
    if start_range != 3500:
        end_range = range_list[index + 1]

        #filters the experimental peaks by the range
        filter_act_peaks = act_peaks[act_peaks.between
                                                    (start_range, end_range, "left")]
        count = len(filter_act_peaks)
        count_list.append((start_range, end_range, count))
    index += 1
count_df = pd.DataFrame(count_list, columns= ["start", "end", "count"])
#normalises the count by dividing by the sum of all counts
count_df["Normalised_Count"] = count_df["count"] / count_df["count"].sum()

# picking out top 20 species for mowse like scoring system
# as scoring system computationally expensive
species_file_list = []
for species in match_results_df["SPECIES"].head(20):
    species = species.replace(" ", "_")
    species_file = "{0}_filtered.csv".format(species)
    species_file_list.append(species_file)

# filter the CSV files by the top 20 species
# calculates the mowse like score
# score is value in dictionary
# species name is the key
score_dict = {}
for csv in csv_files:
    csv_name = csv.split("\\")[-1]
    if csv_name in species_file_list:
        theor_peaks_df = pd.read_csv(csv, sep = ",", dtype = dtype, 
                                 usecols = usecols)

        mowse_score = mowse(count_df, act_peaks_df, theor_peaks_df, range_list)

        csv_names = csv_name.split("_")
        species_name = csv_names[0] + " " + csv_names[1]
        # if it has three words in name
        # and third word doesn't contain .csv
        if csv_names[2].endswith(".csv") == False:
            species_name += " " + csv_names[2]
        score_dict[species_name] = mowse_score

# converts that dictionary into a dataframe
score_df = pd.DataFrame.from_dict(score_dict, orient = "index", columns = ["Score"])

# merges the initial match dataframe with the mowse score dataframe
final_df = match_results_df.merge(score_df, how="inner", left_on= "SPECIES", right_index= True)
final_df["Score"] = final_df["Score"] / final_df["Score"].max()

# Sorts by the mowse like acor first and the match second
final_df = final_df.sort_values(by= ["Score", "Match"], ascending= False)

# outputs the top ten results
# and outputs the results to the filepath chosen
print(final_df.head(10))
output_path = args.output
match_results_df.to_csv(output_path)









