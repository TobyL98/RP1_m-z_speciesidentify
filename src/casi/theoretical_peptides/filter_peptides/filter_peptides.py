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
import glob
import sys
from pathlib import Path

import pandas as pd

##########################
## FUNCTIONS
#########################


def integrate(output_path):
    print("Filter theoretical peptides by LCMSMS data.")
    print("Final Output csv files are in filtered_peptides folder.\n")

    # loading in the theoretical peptides from sequence
    # load in data
    # load in all csv files in in_silico_res folder
    folder_path = output_path / "unfiltered_peptides"
    csv_files = folder_path.glob("*.csv")

    for csv in csv_files:
        predict_df = pd.read_csv(csv, sep=",")

        # changing column name
        new_names = {
            "seq": "pep_seq",
            "nhyd": "hyd_count",
            "ndeam": "deam_count",
            "seq_start": "pep_start",
            "seq_end": "pep_end",
        }
        predict_df = predict_df.rename(columns=new_names)

        # creating a pep_diff so merge excludes correct rows
        predict_df["pep_diff"] = predict_df["pep_end"] - predict_df["pep_start"]

        # loading the LC-MS/MS data from peptide_rules.py
        lcmsms_input_fp = output_path / "lcmsms_masses.csv"
        lcmsms_df = pd.read_csv(lcmsms_input_fp, sep=",")

        # creating columns for a range around pep_start and pep_end
        # can then filter between these ranges
        lcmsms_df["pep_start_max"] = lcmsms_df["pep_start"] + 4
        lcmsms_df["pep_start_min"] = lcmsms_df["pep_start"] - 4
        lcmsms_df["pep_end_max"] = lcmsms_df["pep_end"] + 4
        lcmsms_df["pep_end_min"] = lcmsms_df["pep_end"] - 4

        # creating pep_diff to include in merge
        lcmsms_df["pep_diff"] = lcmsms_df["pep_end"] - lcmsms_df["pep_start"]

        # drop columns we don't need for merge
        lcmsms_df = lcmsms_df.drop(
            columns=[
                "pep_start",
                "pep_end",
                "Unnamed: 0",
                "index",
                "pep_seq",
                "prot_acc",
                "PMF_predict",
                "pep_exp_mr",
                "pep_miss",
            ]
        )
        # merging two datasets
        predict_lc_df = pd.merge(
            predict_df,
            lcmsms_df,
            on=["hyd_count", "deam_count", "pep_diff"],
            how="inner",
        )

        # filtering the values so that we only have rows where the pep_starts
        # and pep_ends are within +- 4
        predict_lc_df = predict_lc_df[
            (predict_lc_df["pep_start"] >= predict_lc_df["pep_start_min"])
            & (predict_lc_df["pep_start"] <= predict_lc_df["pep_start_max"])
            & (predict_lc_df["pep_end"] >= predict_lc_df["pep_end_min"])
            & (predict_lc_df["pep_end"] <= predict_lc_df["pep_end_max"])
        ]

        predict_lc_df = predict_lc_df.sort_values(
            by=["pep_id", "pep_score"], ascending=False
        )

        # drop duplicates with same sequence and hyd and deam count
        subset_list = ["pep_seq", "hyd_count", "deam_count", "mass1"]
        predict_lc_df = predict_lc_df.drop_duplicates(subset=subset_list)

        # drop columns not needed after filtering
        predict_lc_df = predict_lc_df.drop(
            columns=[
                "pep_start_max",
                "pep_start_min",
                "pep_end_min",
                "pep_end_max",
                "Unnamed: 0",
            ]
        )

        predict_lc_df.reset_index(inplace=True, drop=True)

        # naming csv
        species_name = predict_lc_df.loc[0, "SPECIES"]
        species_name = species_name.replace(" ", "_")
        # output_path is input into function
        output_folder = output_path / "filtered_peptides"
        output_folder.mkdir(exist_ok=True)
        csv_name = f"{species_name}_col1peptides_filt.csv"
        csv_filepath = output_folder / csv_name
        predict_lc_df.to_csv(csv_filepath)
    print("######################################")


if __name__ == "__main__":
    sys.exit()
