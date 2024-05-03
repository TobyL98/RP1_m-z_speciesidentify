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


import pandas as pd
import glob
import subprocess
import argparse
import time
from pathlib import Path
import sys

## importing scripts used in pipeline
from Sequences.fasta_A1_clean_NCBI import cleanA1
from Sequences.fasta_A2_clean_NCBI import cleanA2
from Sequences.fasta_seq_amendA1A2_NCBI import COLA1A2combine
from Sequences.fasta_to_csv_NCBI import FastaToCSV
from PTM_rules.Integration_code.integrate_NCBI import integrate

# tests if the input file exists
def file_test(arg):
    p = Path(arg)
    if p.is_file():
        return p
    else:
        raise FileNotFoundError(arg)

# test if the input directory exists
def directory_test(arg):
    p = Path(arg)
    if p.is_dir():
        return p
    else:
        raise Exception("The input directory does not exist {0}". format(p))


def parse_args():
    description = """"Code generates theoretical COL1 peptides for species from there COL1A1 and COL1A2 sequences in a database such as NCBI.
    The inputs are COL1A1 and COL1A2 sequence fasta files downloaded from the database. The output are the theoretical peptide 
    csvs. The ouput csvs can the be used to compare against ZooMS experimental PMFs in compare_NCBI.py"""
    # set up argparse
    parser = argparse.ArgumentParser()
    # arguments
    parser.add_argument("-iA1", "--inputA1",
                        help = """The input COL1A1 peptide sequences fasta file to be used. Needs to contain taxonomic informations.
                        Default is 'Sequences/COL1A1_seqs_NCBI.fasta'""",
                        default= 'Sequences/COL1A1_seqs_NCBI.fasta',
                        type = file_test) # adds the input COL1A1 file
    parser.add_argument("-iA2", "--inputA2",
                        help = """The input COL1A2 peptide sequences fasta file to be used. Needs to contain taxonomic informations.
                        Default is 'Sequences/COL1A2_seqs_NCBI.fasta'"""
                        , default= 'Sequences/COL1A2_seqs_NCBI.fasta'
                        , type = file_test) # adds the input COL1A2 file
    # adds folder where theoretical peptide results should be outputted
    parser.add_argument("-o", "--output",
                        help = """The output folder for where the theoretical peptide csvs will be saved.
                        Default is 'PTM_rules/Integration_code/results_NCBI'"""
                        , default= "PTM_rules/Integration_code/results_NCBI"
                        , type = directory_test)
    parser.add_argument("-r", "--Rscript",
                        help = "The file path of Rscript.exe. Default is 'C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe'",
                        default= "C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe",
                        type = file_test)
    args = parser.parse_args()  # parse arguments
    return(args)

def main(argv):
    args = parse_args()

    # cleans the COL1A1 sequences provided
    print("STEP 1:")
    A1_file = Path(args.inputA1)
    cleanA1(A1_file)

    # cleans the COL1A2 sequences provided
    print("STEP 2:")
    A2_file = Path(args.inputA2)
    cleanA2(A2_file)

    # Combines COLA1 and COL1A2 and adds taxonomic information
    # Outputs as Sequences/COL1A1A2_combined_seqs.fasta
    print("STEP 3:")
    COLA1A2combine()

    # Converts the fasta file with COL1A1 COL1A2 combined into csv
    print("STEP 4:")
    FastaToCSV()

    # run the R script run_parseseq.R
    # this generates all possible theoretical peptides
    # slow part of process!
    print("STEP 5:")
    command = Path(args.Rscript)
    path2Script = Path(r"insilico_digest/run_parseseq_NCBI.R")
    cmd = [command, path2Script]

    # tries to run the R subprocess but will return the error for what is going wrong
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        print(e)

    # integrates the theoretical peptides generated from run_parseseq.R
    # with the LCMSMS data
    # to generate final theoretical peptides
    print("STEP 6:")
    output_folder = args.output
    integrate(output_folder)

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))


