#!/bin/bash


for file in inputs/PinHoleDataset_SIMPER/Bovidae/peak_list/*
do
    echo "Running code with input file: $file"
    output_name=$(basename "$file" peaklist.txt)result.csv
    python compare_NCBI.py -ip $file -o "outputs/PinHoleDataset_SIMPER/Bovidae/$output_name"
done