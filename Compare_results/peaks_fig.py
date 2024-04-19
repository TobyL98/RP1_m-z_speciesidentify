import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path


# reads in the experimental PMF text file with m/z values
def read_exp_PMF(input_PMF_path):

    # read in txt file of PMF values from data
    dtype= {"MZ": 'float32', "intensity": 'float32'}
    act_peaks_df = pd.read_table(input_PMF_path, sep = "\t", header = None, names = ["MZ", "intensity"], dtype= dtype)
    act_peaks_df["intensity"] = 1
    act_peaks_df = act_peaks_df.loc[(act_peaks_df["MZ"] > 3000) & (act_peaks_df["MZ"] < 3150)]
    return(act_peaks_df)

def read_theor_PMF(theor_PMF_path):

    dtype = {"mass1": 'float32'}
    usecols = ["mass1"]
    theor_peaks_df = pd.read_csv(theor_PMF_path, sep = ",", dtype = dtype, 
                                     usecols = usecols)
    theor_peaks_df["intensity"] = 1
    theor_peaks_df = theor_peaks_df.loc[(theor_peaks_df["mass1"] > 3000) & (theor_peaks_df["mass1"] < 3150)]
    print(theor_peaks_df.head())
    return(theor_peaks_df)


def plot_peaks(exp_PMF, theor_PMF):
    plt.style.use('_mpl-gallery')

    fig, (ax1, ax2) = plt.subplots(1, 2, sharex= True)
    #fig.tight_layout()
    ax1.stem(exp_PMF["MZ"], exp_PMF["intensity"], markerfmt= " ", basefmt = "black")
    ax1.set_xlabel("m/z", fontsize = 25)
    ax1.set_title("Experimental ZooMS Spectrum", fontsize = 28)
    ax1.set_yticks([])
    ax1.tick_params(axis = 'x', labelsize=20)
    ax1.grid(False)

    
    ax2.stem(theor_PMF["mass1"], theor_PMF["intensity"], markerfmt= " ", basefmt = "black")
    ax2.set_xlabel("m/z", fontsize = 25)
    ax2.set_title("Theoretical Spectrum", fontsize = 28)
    ax2.set_yticks([])
    ax2.tick_params(axis = 'x', labelsize=20)
    ax2.grid(False)
    
    

    plt.subplots_adjust(bottom= 0.1, right = 0.8, left = 0.1, top = 0.9)
    plt.show()


# Read in sample PMF (goat)
goat_exp_PMF = Path(r"PMF_samples/Artiodactyla/Capra_hircus_sample.txt")
exp_PMF = read_exp_PMF(goat_exp_PMF)

p = Path().absolute()
goat_theor_path = p.parent / r"theoretical_peptides_pipeline/PTM_rules/Integration_code/results_NCBI/Capra_hircus_filt.csv"
theor_PMF = read_theor_PMF(goat_theor_path)

plot_peaks(exp_PMF, theor_PMF)




