#!/bin/bash

input_folder=data/inputs/fasta_input
input_col1a1=$input_folder/COL1A1_seqs_NCBI.fasta
input_col1a2=$input_folder/COL1A2_seqs_NCBI.fasta

output_folder=data/outputs/theor_pep_outputs
python_file=source_code/species_inference/theoretical_peps.py

python $python_file -iA1 $input_col1a1 -iA2 $input_col1a2 -o $output_folder

