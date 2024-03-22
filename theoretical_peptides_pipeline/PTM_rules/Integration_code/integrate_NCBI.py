#########################
# Code - integrate_NCBI.py
#########################

# same code but for NCBI data
# so no taxonomic information
# combines the data from LCMSMS peptide_rules.py
# and theoretical peptides from sequence from run_parseseq_NCBI.R
# to create theoretical peptides that can be used
# to compare against PMFs to identify species
# theoretical peptides from run_parseseq.R need to be in
# folder in_silico_res_NCBI as csv files


## Loading modules
import pandas as pd
import matplotlib.pyplot as plt
import glob

##########################
## FUNCTIONS
#########################

def integrate(output_path):
    print("Filter theoretical peptides by LCMSMS data.")
    print("Final Output files are in integrate_results folder. Names are\n")

    # loading in the theoretical peptides from sequence
    # load in data
    # load in all csv files in in_silico_res folder
    csv_files = glob.glob('C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/in_silico_res_NCBI/*.csv')

    for csv in csv_files:
        predict_df = pd.read_csv(csv, sep = ',')
    
    
        #changing column name
        new_names = {"seq": "pep_seq", "nhyd": "hyd_count", "nglut": "deam_count", "seq_start": "pep_start", "seq_end": "pep_end"}
        predict_df = predict_df.rename(columns = new_names)

        # creating a pep_diff so merge excludes correct rows
        predict_df["pep_diff"] = predict_df["pep_end"] - predict_df["pep_start"]

        # loading the LC-MS/MS data from peptide_rules.py
        LCMSMS_df = pd.read_csv("C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/peptide_rules/sequence_masses.csv", sep = ',')

        # creating columns for a range around pep_start and pep_end
        # can then filter between these ranges
        LCMSMS_df["pep_start_max"] = LCMSMS_df["pep_start"] + 4
        LCMSMS_df["pep_start_min"] = LCMSMS_df["pep_start"] - 4
        LCMSMS_df["pep_end_max"] = LCMSMS_df["pep_end"] + 4
        LCMSMS_df["pep_end_min"] = LCMSMS_df["pep_end"] - 4

        # creating pep_diff to includ in merge
        LCMSMS_df["pep_diff"] = LCMSMS_df["pep_end"] - LCMSMS_df["pep_start"]

        # drop columns we don't need for merge
        LCMSMS_df = LCMSMS_df.drop(columns = ["pep_start", "pep_end", "Unnamed: 0", "index", "pep_seq",
                                              "prot_acc", "PMF_predict", "pep_exp_mr", "pep_miss"])

        # merging two datasets
        predict_LC_df = pd.merge(predict_df, LCMSMS_df, 
                                 on = ["hyd_count", "deam_count", "pep_diff"], how = 'inner')

        # filtering the values so that we only have rows where the pep_starts and pep_ends are within +- 4
        predict_LC_df = predict_LC_df[(predict_LC_df["pep_start"] >= predict_LC_df["pep_start_min"]) &
                                      (predict_LC_df["pep_start"] <= predict_LC_df["pep_start_max"]) &
                                      (predict_LC_df["pep_end"] >= predict_LC_df["pep_end_min"]) &
                                      (predict_LC_df["pep_end"] <= predict_LC_df["pep_end_max"])]

        predict_LC_df = predict_LC_df.sort_values(by = ["pep_id", "pep_score"], ascending = False)

        # drop duplicates with same sequence and hyd and deam count
        subset_list = ["pep_seq", "hyd_count", "deam_count", "mass1"]
        predict_LC_df = predict_LC_df.drop_duplicates(subset= subset_list)

        # drop columns not needed after filtering
        predict_LC_df = predict_LC_df.drop(columns = ["pep_start_max", "pep_start_min",
                                                      "pep_end_min", "pep_end_max", "Unnamed: 0"])

        predict_LC_df.reset_index(inplace= True, drop = True)

        # naming csv
        species_name = predict_LC_df.loc[0, "SPECIES"]
        species_name = species_name.replace(" ", "_")
        print("{0}_filtered.csv".format(species_name))
        # output_path is input into function
        csv_name =  "{0}/{1}_filtered.csv".format(output_path, species_name)
        predict_LC_df.to_csv(csv_name)
    print("######################################")

