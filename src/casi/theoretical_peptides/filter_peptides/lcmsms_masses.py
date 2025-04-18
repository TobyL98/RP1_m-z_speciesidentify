#########################
# Code - peptide_rules.py
#########################

# Code uses the LC-MSMS data from multiple species
# collagen A1 and A2
# groups the data by peptides that have the same start and end positions
# presumably the same or very similar pepetides due to lack variation
# in collagen
# uses these groups to work out common oxidation and deamidations in the peptide
# these can then be applied to theoretically (in silico) generated peptides
# to predict what likely PTMs are at that peptide
# OUTPUTS: sequence_masses.csv that is used in integrate.py

## Loading modules
import glob
import re
from pathlib import Path
import sys

import pandas as pd

###################
###  Functions  ###
###################


# This function uses the Pep-var_mod column
# to work out the number of PTMS in targeted residue (P, K, N/Q)
# that there is in each peptide fragment sequenced by LC-MS/MS
def mod_count(mods, res):
    mods = str(mods)

    mod_count = 0
    mods_list = []
    if re.search(r";", mods):
        # splits if it has ;
        # then seperates target residue modifications from other mods
        # example input: Oxidation (K); 2 Oxidation (P)
        all = mods.split(";")
        for x in all:
            if re.search(res, x):
                mod = x.strip(" ")
                mods_list.append(mod)
    # if only one mod still appends to a list
    else:
        mod = mods.strip(" ")
        mods_list.append(mod)

    for mod in mods_list:
        # if it contains target residue
        if re.search(res, mod):
            mod_count += 1  # assumes it will be 1
            mod_num = mod.split(" ")[0]  # splits to get number
            # if the number matches a number between 2 and 10
            # assigns that number to count, overwrites 1
            for num in range(2, 11):
                if mod_num == str(num):
                    mod_count += num
                    mod_count -= 1  # removes 1 so initial 1 not counted

    return mod_count


# Loads in all the LCMSMS mascot csv files from the LCMSMS folder
# and concatenates into one single dataframe
def data_load(file_path):
    csv_files = file_path.glob("*.csv")

    # loops through all CSV files
    # creates dataframe with columns required
    # adds df to list
    df_list = []
    dtypes = {"pep_score": "float32", "pep_exp_mr": "float32"}
    use_cols = [
        "pep_seq",
        "pep_score",
        "pep_start",
        "pep_end",
        "pep_exp_mr",
        "pep_miss",
        "pep_var_mod",
        "prot_acc",
    ]
    for csv in csv_files:
        new_df = pd.read_csv(csv, sep=",", dtype=dtypes, usecols=use_cols)

        # drop all empty rows
        new_df.dropna(how="all", inplace=True)
        # change start and end columns to int
        new_df = new_df.astype({"pep_start": int, "pep_end": int})

        df_list.append(new_df)

    # combines all dataframes in the list into single df
    df = pd.concat(df_list)

    pos_sorted_df = df.sort_values(by=["pep_start", "pep_end"])
    return pos_sorted_df


# reads in the dataframe with all possible PTMs
# groups the data by peptides that have the same start and end positions
# presumably the same or very similar pepetides due to lack variation
# in collagen
# uses these groups to work out common oxidation and deamidations in the peptide
def df_filter(df):
    # create a startend list
    # means don't repeat sequences in dataframe
    pep_startend_list = []

    pep_seq_df_list = []
    seq_count = 0

    # This code loops through rows
    #  finds specific start and end and sequence length
    #  filters the dataframe by the specific start and end
    for index, row in df.iterrows():
        pep_start = row["pep_start"]
        # creates a start list based on peptide start
        # allows two behind and two after to account for frameshifts
        start_list = list(range(pep_start - 4, pep_start + 4, 1))
        pep_end = row["pep_end"]
        # creates an end list based on peptide start
        end_list = list(range(pep_end - 4, pep_end + 4, 1))
        # gets the sequences length
        seq_length = pep_end - pep_start
        # puts into a tuple which accounts for all three factors
        startend = (pep_start, pep_end, seq_length)

        # if tuple value is not already in list
        if startend not in pep_startend_list:
            seq_count += 1

            # filters df by start_list, end list and seq_length
            pep_seq_df = df.loc[
                df["pep_start"].isin(start_list)
                & df["pep_end"].isin(end_list)
                & (df["pep_end"] - df["pep_start"] == seq_length)
            ]

            # dropping duplicates so only keep the top pep_score
            pep_seq_df = pep_seq_df.sort_values(by=["pep_score"], ascending=False)
            pep_seq_df = pep_seq_df.drop_duplicates(
                subset=[
                    "pep_seq",
                    "pep_start",
                    "pep_end",
                    "pep_miss",
                    "pep_var_mod",
                    "prot_acc",
                ]
            )
            pep_seq_df = pep_seq_df.sort_values(by=["prot_acc"])
            pep_seq_df = pep_seq_df.reset_index(drop=True)
            pep_seq_df["pep_id"] = seq_count  # add a unique identifier

            # use apply on function to count number of hydroxylations (M, P and K)
            # and demidations (N and Q)
            # from LC-MS/MS data
            pep_seq_df["hyd_count"] = pep_seq_df["pep_var_mod"].apply(
                mod_count, res=r"[MPK]"
            )
            pep_seq_df["deam_count"] = pep_seq_df["pep_var_mod"].apply(
                mod_count, res=r"NQ"
            )

            # group by used to ensure there are no duplicate PTMs are included
            # hyd_count and deam_count are important for this
            # summarises prot_acc to include all the species that have this PTM
            # summarises other numbers
            pep_seq_df = pep_seq_df.groupby(
                ["pep_seq", "pep_miss", "hyd_count", "deam_count", "pep_id"],
                as_index=False,
                dropna=False,
            ).agg(
                {
                    "prot_acc": ", ".join,
                    "pep_score": "max",
                    "pep_exp_mr": "mean",
                    "pep_start": "min",
                    "pep_end": "min",
                }
            )

            # add 1 to calculate PMF value
            pep_seq_df["PMF_predict"] = pep_seq_df["pep_exp_mr"] + 1
            # adds the df to a dictionary
            pep_seq_df_list.append(pep_seq_df)

            # generates the pep_startend_list
            # adds all peptide start, end and sequence length combinations that have been done
            # means they are not done more than once
            num_count = 0
            for num in start_list:
                startend_poss = (start_list[num_count], end_list[num_count], seq_length)
                pep_startend_list.append(startend_poss)
                num_count += 1

    # add all the different peptide dfs into one df
    all_peps_df = pd.concat(pep_seq_df_list)
    return all_peps_df


def mass_lcsmsms(lcmsms_dir, output_folder):
    print("Generating LCMSMS Filter Rules")
    if not lcmsms_dir.is_dir():
        raise FileNotFoundError(f"""The directory containing the LCMSMS data does not exist {lcmsms_dir}.
Ensure the following directory is created and put the LCMSMS data in it""")
    # read LCMSMS CSV files and merge to one dataframe
    all_df = data_load(lcmsms_dir)

    # filter to obtain desired output of likely PTMs
    final_peps_df = df_filter(all_df)

    # reorder for presentation
    correct_order = [
        "pep_id",
        "pep_seq",
        "pep_start",
        "pep_end",
        "pep_exp_mr",
        "hyd_count",
        "deam_count",
        "pep_miss",
        "prot_acc",
        "pep_score",
        "PMF_predict",
    ]
    final_peps_df = final_peps_df.reindex(columns=correct_order)
    final_peps_df = final_peps_df.reset_index()

    # save as a csv file
    output_file = output_folder / "lcmsms_masses.csv"
    final_peps_df.to_csv(output_file, sep=",")
    print(f"Output: {output_file}")
    print("######################################")


if __name__ == "__main__":
    sys.exit()
