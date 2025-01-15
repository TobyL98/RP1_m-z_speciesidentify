#!/bin/bash

for folder in inputs/PinHoleDataset_SIMPER/*
do
    output_folder=outputs/PinHoleDataset_SIMPER_output/$(basename "$folder")
    mkdir -p "$output_folder"

    for file in "$folder"/peak_list/*
    do
        echo "Running code with input file: $file"
        output_name=$(basename "$file" peaklist.txt)result.csv
        python compare_NCBI.py -ip "$file" -o "$output_folder"/"$output_name"
    done
done