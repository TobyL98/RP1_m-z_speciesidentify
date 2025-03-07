#!/bin/bash

input_folder=data/inputs
input_col1a1=$input_folder/AvesCOL1A1_240225.txt
input_col1a2=$input_folder/AvesCOL1A2.txt

output_folder=data/bird_outputs

theoretical_peps -ia1 $input_col1a1 -ia2 $input_col1a2 -o $output_folder -sc birds

