#!/bin/bash


for file in inputs/PinHoleDataset_SIMPER/Cervidae/peak_list/*
do
    echo "Running code with input file: $file"
    output_name=$(basename "$file" peaklist.txt)result.csv
    python compare_NCBI.py -ip $file -o "outputs/PinHoleDataset_SIMPER_output_run2/Cervidae/$output_name"
done