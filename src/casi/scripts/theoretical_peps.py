##########################
# theoretical_peps_NCBI.py
##########################

# Same as original process but for files COL1A1 and COl1A2
# files directly from NCBI
# This runs the whole process of generating
# theoretical peptides that can then
# be used to compare against PMF
# starting files need to be:
# fasta file of COL1A1 with taxonomic information
# fasta file of COL1A2 with taxonomic information

import argparse
import sys
from pathlib import Path
from importlib.resources import files

import pandas as pd

from casi.theoretical_peptides.sort_sequences.fasta_col_clean import run_clean_col
from casi.theoretical_peptides.sort_sequences.merge_cola1a2 import col1a1a2_combine
from casi.theoretical_peptides.sort_sequences.fasta_to_csv import FastaToCSV
from casi.theoretical_peptides.generate_peptides.cleave_all_sequences import collagen_peptide_mass
from casi.theoretical_peptides.filter_peptides.lcmsms_masses import mass_lcsmsms
from casi.theoretical_peptides.filter_peptides.filter_peptides import integrate


def file_test(arg):
    """tests if the input file exists"""
    p = Path(arg)
    if p.is_file():
        return p
    else:
        raise FileNotFoundError(arg)


def directory_test(arg):
    """test if the input directory exists"""
    p = Path(arg)
    if p.is_dir():
        return p
    else:
        raise Exception("The directory does not exist: {0}".format(p))
    
def import_lcsmsms(lcmsms_arg):
    """Imports the default Mascot LC-MS/MS output csv files from the package.
    if the user has defined a directory path, this will be used instead

    Args:
        lcsmsms_input_dir: Absolute path to directory or None
    Returns:
        lcmsms_dir (directory path): Absolute path to the directory
        containing LC-MS/MS data
    """
    if lcmsms_arg is not None:
        lcmsms_dir = Path(lcmsms_arg)
    else: # get default from package
        lcmsms_dir = Path(files("casi.theoretical_peptides.input_data.lcmsms"))

    return lcmsms_dir



def parse_args(argv):
    """Command line arguments"""
    info = """Code generates theoretical COL1 peptides for species
    from there COL1A1 and COL1A2 sequences in a database such as NCBI.
    The inputs are COL1A1 and COL1A2 sequence fasta files downloaded from the database. 
    The output are the theoretical peptide csvs. The ouput csvs can the be used to 
    compare against ZooMS experimental PMFs in compare_NCBI.py"""
    # set up argparse
    parser = argparse.ArgumentParser(description=info)
    # arguments
    parser.add_argument(
        "-ia1",
        "--inputa1",
        help="""The input COL1A1 peptide sequences fasta file to be used.
        Needs to contain taxonomic informations.""",
        type=file_test,
        required=True
    )  # adds the input COL1A1 file
    parser.add_argument(
        "-ia2",
        "--inputa2",
        help="""The input COL1A2 peptide sequences fasta file to be used.
        Needs to contain taxonomic informations.""",
        type=file_test,
        required=True
    )  # adds the input COL1A2 file
    # adds folder where theoretical peptide results should be outputted
    parser.add_argument(
        "-o",
        "--output",
        help="The output folder for where all the outputs will be saved.",
        type=directory_test,
        required=True
    )
    parser.add_argument(
        "-lc",
        "--lcmsms",
        help="""Folder containing the mascot LC-MS/MS csv output files.
        Default LC-MS/MS csv files from the Buckley lab will be used 
        to identify mammals (if no input given).""",
        type=directory_test
    )
    parser.add_argument(
        "-sc",
        "--species_class",
        help="""The class that your species is in. There are currently only two answers
        to this: 'birds' or 'mammals' as the pipeline currently works for birds and mammals.""",
        default="mammals",
        choices=["birds", "mammals"]
    )
    args = parser.parse_args()  # parse arguments
    return args


def main(argv=sys.argv[1:]):
    args = parse_args(argv)

    # cleans the COL1A1 sequences provided
    print("STEP 1:")
    a1_file = Path(args.inputa1)
    output_folder = Path(args.output)
    run_clean_col(a1_file, output_folder, "COL1A1")

    # cleans the COL1A2 sequences provided
    print("STEP 2:")
    class_input = str(args.species_class)
    a2_file = Path(args.inputa2)
    run_clean_col(a2_file, output_folder, "COL1A2", class_input)

    # Combines COLA1 and COL1A2 and adds taxonomic information
    # Outputs as Sequences/COL1A1A2_combined_seqs.fasta
    print("STEP 3:")
    col1a1a2_combined = col1a1a2_combine(output_folder)

    # Generates all possible theoretical peptides and their masses
    print("STEP 4:")
    collagen_peptide_mass(col1a1a2_combined, output_folder)

    # formatting possible LCMSMS masses into one document
    # used to then filter theoretical peptides
    print("Step 5:")
    lcmsms_dir = import_lcsmsms(args.lcmsms)
    mass_lcsmsms(lcmsms_dir, output_folder)

    # integrates the theoretical peptides generated
    # with the LCMSMS data
    # to generate final theoretical peptides
    print("STEP 6:")
    output_folder = args.output
    integrate(output_folder)


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
