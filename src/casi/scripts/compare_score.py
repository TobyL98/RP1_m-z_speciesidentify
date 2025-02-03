"""
Compare_correlation_NCBI.py

 Gets the theoretical peptides generated from the sequences
 and filtered by the LCMSMS data
 and compares how many match an actual PMF within a certain tolerance
 tolerance is usually 0.2
 returns the results with the number of matches for each species
 ordered by highest number of matches
 highest number of matches should be the species most closely related
 to species in the database
"""

import sys
from pathlib import Path
import argparse

import pandas as pd

################
# FUNCTIONS
################

def file_test(arg):
    """Tests if the input file exists"""
    p = Path(arg)
    if p.is_file():
        return p
    else:
        raise FileNotFoundError(arg)

def directory_test(arg):
    """Test if the input directory exists"""
    p = Path(arg)
    if p.is_dir():
        return p
    else:
        raise Exception("The input directory does not exist {0}".format(p))

def output_test(arg):
    """Test if directory of new output file exists"""
    p = Path(arg)
    par = p.parent
    if par.is_dir():
        return p
    else:
        raise Exception(
            "The directory of the new output file does not exist {0}".format(p)
        )


# test if --top5 input is 0 or 1 only
def test_01(arg):
    arg = int(arg)
    if arg == 0 or arg == 1:
        return arg
    else:
        raise Exception("The input -m5 --top5 should be the integer 0 or 1 only")


def parse_args(argv):
    description = """This code compares the experimental PMF peaks from a sample
    and the theoretical peaks generated from the all the species COL1 theoretical 
    peptides generated in the theoretical peptides pipeline. It will score a match 
    if the theoretical peptide m/z valus is within a specified threshold (+- 0.2 is default).
    The output to the command line will be the top 10 match scores (number of matches) species.
    All species match scores wil be outputted to a CSV file. """

    # set up argparse
    parser = argparse.ArgumentParser(description=description)
    # add arguments
    # adds the folder where the input theoretical peak spectrums are
    parser.add_argument(
        "-it",
        "--inputTheor",
        help="the folder that contains the theoretical peptides csv files to compare against PMF",
        type=directory_test,
    )
    # adds where the output file should be saved
    parser.add_argument(
        "-o",
        "--output",
        help="The file name for the output results file of number of matches",
        type=output_test,
    )
    # adds where the input peptide mass fingerprint
    parser.add_argument(
        "-ip",
        "--inputPMF",
        help="The input Peptide mass fingerprint (PMF) from an unknown organism.",
        type=file_test,
    )
    parser.add_argument(
        "-t",
        "--threshold",
        help="The threshold for matches between the experimental and theoretical spectrum. Default is 0.2 Da",
        default=0.2,
        type=float,
    )
    parser.add_argument(
        "-m5",
        "--top5",
        help="""If 1 is inputted will provide an excel file of m/z peak matches for the top 5 match counts as an xlsx (excel) file.
                        Default is 0""",
        default=0,
        type=test_01,
    )
    args = parser.parse_args()
    return args

def read_exp_PMF(input_PMF):
    """reads in the experimental PMF text file with m/z values"""
    # read in txt file of PMF values from data
    dtype = {"MZ": "float32", "intensity": "float32"}
    act_peaks_df = pd.read_table(
        input_PMF, sep="\t", header=None, names=["MZ", "intensity"], dtype=dtype
    )
    return act_peaks_df

def read_theor_csv(input_theor_path):
    """reads in all the csv files theoretical peptide m/z values
    and assigns all as dataframes to a list"""
    # get all csv files
    csv_files = input_theor_path.glob("*.csv")

    # read in all csvs of theoretical peptide peaks
    dtype = {
        "mass1": "float32",
        "GENUS": "category",
        "SPECIES": "category",
        "pep_seq": "category",
    }
    usecols = [
        "mass1",
        "GENUS",
        "SPECIES",
        "pep_seq",
        "pep_start",
        "pep_end",
        "hyd_count",
        "deam_count",
        "missed_cleaves",
    ]
    theor_peaks_df_list = []
    for csv in csv_files:
        theor_peaks_df = pd.read_csv(csv, sep=",", dtype=dtype, usecols=usecols)
        theor_peaks_df_list.append(theor_peaks_df)
    return theor_peaks_df_list


# function does the comparison between one set of theoretical peptides
# and the PMF within a certain allowance
# theor_peaks are the theoretical peaks
# act_peaks are the actual peaks from PMF
def compare(theor_peaks, act_peaks, threshold):
    """Function does the comparison between one set of theoretical peptides
    and the PMF within a certain allowance theor_peaks are the 
    theoretical peaks act_peaks are the actual peaks from PMF"""
    # creates columns that are the same in each
    # can then be merged
    act_peaks["join"] = 1
    theor_peaks["join"] = 1

    # generating allowance for comparison
    # with theoretical peaks
    act_peaks["MZ_plus"] = act_peaks["MZ"] + threshold
    act_peaks["MZ_minus"] = act_peaks["MZ"] - threshold

    # merges the datasets so we have all the data
    # then filters where there is a MATCH within the tolerance
    # number of rows left will be number of matches
    merged_df = act_peaks.merge(theor_peaks, how="outer", on=["join"])
    merged_df = merged_df.loc[
        (merged_df["mass1"] >= merged_df["MZ_minus"])
        & (merged_df["mass1"] <= merged_df["MZ_plus"])
    ]

    # remove duplicates if MZ value has matched more than once
    matches_df = merged_df.drop_duplicates(subset=["MZ"])

    # assigns number of matches to count
    match_count = matches_df.shape[0]

    # turns match results to a df
    result_df = pd.DataFrame([match_count], columns=["Match"])

    # combines with the taxon information to identify which species it is
    taxon_df = theor_peaks.loc[[0], ["GENUS", "SPECIES"]]
    final_df = pd.concat([taxon_df, result_df], axis=1)
    return (final_df, matches_df, match_count)

def peaks_comparison(theor_peaks_list, act_peaks_df, thresh, output):
    """Reads in the theoretical peaks and actual peaks
    and runs compare function organises correct output"""
    results_list = []
    matches_dict = {}
    for theor_peaks in theor_peaks_list:
        # run the compare function
        result_df, matches_df, match_count = compare(theor_peaks, act_peaks_df, thresh)
        # add all results to a list
        results_list.append(result_df)
        # adds all the matches df to a dictionary
        # match count is key
        # value is the df of all the matches
        matches_dict[match_count] = matches_df

    # put all results in one dataframe
    match_results_df = pd.concat(results_list)
    match_results_df = match_results_df.sort_values(by=["Match"], ascending=False)
    match_results_df = match_results_df.reset_index(drop=True)

    # outputs the top matches results
    # and saves to csv
    final_output_df = match_results_df.head(10)
    print("RESULTS:")
    print(final_output_df.to_markdown())
    match_results_df.to_csv(output)
    return matches_dict

def top_5(match_dict, option_match, output):
    """Saves top 5 matches_df as xlxs to same output 
    location as matches if option inputted in command"""
    if option_match == 1:
        top_5_match = sorted(match_dict.keys(), reverse=True)[:5]
        output_parent = output.parent
        top5_path = output_parent / r"top5_matches.xlsx"
        writer = pd.ExcelWriter(top5_path, engine="xlsxwriter")

        count = 0
        for match in top_5_match:
            count += 1
            df = match_dict[match]
            df = df.drop(columns=["join", "MZ_plus", "MZ_minus"])
            df = df.rename(
                columns={
                    "MZ": "Exp MZ",
                    "intensity": "Exp intensity",
                    "mass1": "Theor MZ",
                }
            )
            df = df.reset_index(drop=True)
            df.to_excel(writer, sheet_name="match{0}".format(count))

        writer.close()
        print("\n")
        print(
            "5 species with highest matches m\\z values for all matches have been outputted. File called:"
        )
        print("'top5_matches.xlsx'")

def main(argv=sys.argv[1:]):
    """Main method and logic"""
    print("""####################
Program: Compare_NCBI.py
#####################""")

    args = parse_args(argv)

    input_PMF = Path(args.inputPMF)
    # reads in experimental PMF csv
    actual_peaks_df = read_exp_PMF(input_PMF)

    input_theor_folder = args.inputTheor
    # reads all csvs for species theoretical PMFs
    theoretical_peaks_df_list = read_theor_csv(input_theor_folder)

    output_path = Path(args.output)
    threshold = args.threshold
    print("\nThreshold for match is +- {0}".format(threshold))
    # compares experimental and theoretical PMFs withins a threshold
    matches_dictionary = peaks_comparison(
        theoretical_peaks_df_list, actual_peaks_df, threshold, output_path
    )

    match_opt = args.top5
    # if required outputs the experimental and theoretical peaks that matches
    # for the top 5 matches
    top_5(matches_dictionary, match_opt, output_path)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
