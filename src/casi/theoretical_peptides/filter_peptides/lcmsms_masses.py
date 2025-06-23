"""
peptide_rules.py
Code uses the LC-MSMS data from multiple species collagen A1 and A2 
It groups the data by peptides that have the same start and end positions.
It uses these groups to work out common oxidation and deamidations in the peptide.
These can then be applied to theoretically (in silico) generated peptides
to predict what likely PTMs are in that peptide
OUTPUTS: sequence_masses.csv that is used in integrate.py
"""

import re
from pathlib import Path
import sys
from collections import namedtuple

import pandas as pd

Positions = namedtuple(
    "Positions",
    [
        "start",
        "end",
        "length",
        "start_list",
        "end_list",
    ]
)
 
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


def data_load(dirpath: Path) -> pd.DataFrame:
    """
    Loads in all the LCMSMS mascot csv files from the LCMSMS folder and 
    concatenates them into one single dataframe

    args
        filepath (Path): the path to the directory that contains the LCMSMS csv files

    returns
        lcmsms_df (pd.Dataframe): contains all the LCMSMS data in one dataframe

    """
    csv_files = dirpath.glob("*.csv")

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

    lcmsms_df = df.sort_values(by=["pep_start", "pep_end"])
    return lcmsms_df

def find_positions(row: pd.DataFrame) -> namedtuple:
    """
    Finds the start, end position and length for each peptide
    fragment generated in the LC-MS/MS experiment. The start and end position
    are the start and end positions in the full protein sequence.
    Each row is a different peptide fragment. Creates start and end lists
    that allow for frameshifts in different peptides.

    args
        lcmsms_df (pd.DataFrame) -> dataframe of LC-MS/MS peptide fragments

    returns
        positions (namedtuple) -> Contains the start position, end positions,
        fragment sequence length, start list and end list.
    """

    
    pep_start = row["pep_start"]
    # allows four behind and four after to account for frameshifts
    start_list = list(range(pep_start - 4, pep_start + 4, 1))
    # creates an end list
    pep_end = row["pep_end"]
    end_list = list(range(pep_end - 4, pep_end + 4, 1))
    seq_length = pep_end - pep_start
    
    positions = Positions(
        pep_start,
        pep_end,
        seq_length,
        start_list,
        end_list
    )

    return positions

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

        # NEW FUNCTION

        # if tuple value is not already in list
        if startend not in pep_startend_list:
            seq_count += 1

            # filters df by start_list, end list and seq_length
            pep_seq_df = df.loc[
                df["pep_start"].isin(start_list)
                & df["pep_end"].isin(end_list)
                & (df["pep_end"] - df["pep_start"] == seq_length)
            ]

            # dropping duplicates so only keep the top pep_score for each modification
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

        # NEW FUNCTION

            # use apply on function to count number of hydroxylations (M, P and K)
            # and demidations (N and Q)
            # from LC-MS/MS data
            pep_seq_df["hyd_count"] = pep_seq_df["pep_var_mod"].apply(
                mod_count, res=r"[MPK]"
            )
            pep_seq_df["deam_count"] = pep_seq_df["pep_var_mod"].apply(
                mod_count, res=r"NQ"
            )

        # NEW FUNCTION

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

        # NEW FUNCTION

            # add 1 to calculate PMF value
            pep_seq_df["PMF_predict"] = pep_seq_df["pep_exp_mr"] + 1
            # adds the df to a dictionary
            pep_seq_df_list.append(pep_seq_df)

            # generates the pep_startend_list
            # adds all peptide start, end and sequence length combinations that have been done
            # means they are not done more than once
            num_count = 0 #TODO code needs to be moved if we are refactoring
            for num in start_list:
                startend_poss = (start_list[num_count], end_list[num_count], seq_length)
                pep_startend_list.append(startend_poss)
                num_count += 1         

    # add all the different peptide dfs into one df
    all_peps_df = pd.concat(pep_seq_df_list)
    return all_peps_df

def remove_duplicates(lcmsms_df, positions):
    """
    """

    pep_fragment = lcmsms_df.loc[
        lcmsms_df["pep_start"].isin(positions.start_list)
        & lcmsms_df["pep_end"].isin(positions.end_list)
        & (lcmsms_df["pep_end"] - lcmsms_df["pep_start"] == positions.length)
    ]
    # dropping duplicates so only keep the top pep_score
    pep_fragment = pep_fragment.sort_values(by=["pep_score"], ascending=False)
    pep_fragment = pep_fragment.drop_duplicates(
        subset=[
            "pep_seq",
            "pep_start",
            "pep_end",
            "pep_miss",
            "pep_var_mod",
            "prot_acc",
        ]
    )
    pep_fragment = pep_fragment.sort_values(by=["prot_acc"])
    pep_fragment = pep_fragment.reset_index(drop=True)

    return pep_fragment

def filter_lcmsms(lcmsms_df: pd.DataFrame):
    """
    """

    for index, row in lcmsms_df.iterrows():

        positions = find_positions(row)

        remove_duplicates = filter_lcmsms(lcmsms_df, positions)



    return None


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
