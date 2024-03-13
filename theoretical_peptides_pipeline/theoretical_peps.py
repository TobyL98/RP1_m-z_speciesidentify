##########################
# theoretical_peps.py
##########################

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

## importing scripts used in pipeline
from Sequences.fasta_A1_clean import cleanA1
from Sequences.fasta_A2_clean import cleanA2
from Sequences.fasta_seq_amendA1A2 import COLA1A2combine
from Sequences.fasta_to_csv import FastaToCSV
from PTM_rules.Integration_code.integrate import integrate


#####################
## MAIN CODE
#####################

# set up argparse
parser = argparse.ArgumentParser()
# arguments
parser.add_argument("-iA1", "--inputA1",
                    help = "The input COL1A1 peptide sequences fasta file to be used. Needs to contain taxonomic informations",
                    default= 'Sequences/COL1A1_seqs.fasta') # adds the input COL1A1 file
parser.add_argument("-iA2", "--inputA2",
                    help = "The input COL1A2 peptide sequences fasta file to be used. Needs to contain taxonomic informations"
                    , default= 'Sequences/COL1A2_seqs.fasta') # adds the input COL1A2 file
# adds folder where theoretical peptide results should be outputted
parser.add_argument("-o", "--output",
                    help = "The output folder for where the theoretical peptide csvs will be saved."
                    , default= "C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/PTM_rules/Integration_code/integrate_results")
args = parser.parse_args()  # parse arguments

start_time = time.time()
# cleans the COL1A1 sequences provided
print("STEP 1:")
A1_file = args.inputA1
cleanA1(A1_file)

# cleans the COL1A2 sequences provided
print("STEP 2:")
A2_file = args.inputA2
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
command = 'C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe'
path2Script = 'C:/Users/tobyl/OneDrive - The University of Manchester/Bioinformatics Masters/Research project 1/Git_repositories/RP1_m-z_speciesidentify/theoretical_peptides_pipeline/insilico_digest/run_parseseq.R'
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

end_time = time.time()
overall_time = end_time - start_time
time_mins = overall_time / 60

print(overall_time)
print(time_mins)




