#####################
## Compare_one.py
####################

# SAME as compare.py but only for one theoretical peptide
# function gets the theoretical peptides generated from the sequences
# and filtered by the LCMSMS data
# and compares how many match an actual PMF within a certain tolerance
# tolerance is is usually 0.2
# returns the resuls with the number of matches for each species
# ordered by highest number of matches
# highest number of matches should be the species most closely related
# to species in the database

import pandas as pd
import glob
from scipy import signal
import matplotlib.pyplot as plt
import numpy as np

################
# FUNCTIONS
################

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
def compare(theor_peaks, act_peaks):
    
    MZ_theor_peaks = theor_peaks["mass1"]
    MZ_act_peaks = act_peaks["MZ"]

    peak_list = [(800, 1100), (1100, 1400), (1400, 1700), (1700, 2000),
                 (2000, 2300), (2300, 2600), (2600, 2900), (2900, 3200),
                 (3200, 3500)]

    corr_diff_list = []
    for x in peak_list:
        filter_act_peaks = MZ_act_peaks[MZ_act_peaks.between(x[0], x[1], "left")]

        filter_theor_peaks = MZ_theor_peaks[MZ_theor_peaks.between(x[0], x[1], "left")]
        filter_theor_peaks = filter_theor_peaks.sort_values()

        corr_diff = cross_corr(filter_theor_peaks, filter_act_peaks)
        corr_diff_list.append(corr_diff)
    corr_diff_array = np.concatenate(corr_diff_list)
    corr_diff_mean = np.mean(corr_diff_array)
    print(corr_diff_mean)
    
    

################
# Main code
################

# read in csv of theoretical peptide peaks
dtype = {"mass1": float}
theor_peaks_df = pd.read_csv("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/integrate_results/Felis_catus_filtered.csv", 
                          sep = ",", dtype= dtype)


# read in txt file of PMF values from data
dtype= {"MZ": float}
act_peaks_df = pd.read_table("PMF_samples/Felis_catus_sample.txt", sep = "\t", header = None, names = ["MZ", "intensity"], dtype= dtype)

# run the function
compare(theor_peaks_df, act_peaks_df)

