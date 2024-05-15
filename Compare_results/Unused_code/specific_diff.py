########################
## Specific_diff.py
########################

import pandas as pd
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
    return(merged_df)

# read in txt file of PMF values from data
input_PMF = "PMF_samples/Carnivora/Mirounga_leonina_sample.txt"
dtype= {"MZ": 'float32', "intensity": 'float32'}
act_peaks_df = pd.read_table(input_PMF, sep = "\t", header = None, names = ["MZ", "intensity"], dtype= dtype)

# read in all csvs of theoretical peptide peaks
dtype = {"mass1": 'float32', "GENUS": 'category', "SPECIES": 'category'}
usecols = ["mass1", "GENUS", "SPECIES"]

csv_list = ["C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/results_NCBI/Mirounga_leonina_filtered.csv",
            "C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/results_NCBI/Mirounga_angustirostris_filtered.csv"]
results_list = []

for csv in csv_list:
    theor_peaks_df = pd.read_csv(csv, sep = ",", dtype = dtype, 
                                 usecols = usecols)

    result_df = compare(theor_peaks_df, act_peaks_df)
    # add all results to a list
    results_list.append(result_df)

result_df1 = results_list[0]
result_df2 = results_list[1]

compare_df = result_df1.merge(result_df2, on = 'mass1', indicator = True, how = "outer")
compare_df1 = compare_df[compare_df["_merge"] == "left_only"]
compare_df2 = compare_df[compare_df["_merge"] == "right_only"]

print(compare_df1)
print(compare_df2)