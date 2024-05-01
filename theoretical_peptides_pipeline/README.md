# Generating Theoretical Peptide Spectra from Sequences
## Overview

This pipeline generates theoretical collagen 1 (COL1) PMF spectra from COL1 sequences in a database. The theoretical COL1 PMF spectra can then be compared to experimentally generated spectra

## Usage
```
usage: theoretical_peps_NCBI.py [-h] [-iA1 INPUTA1] [-iA2 INPUTA2] [-o OUTPUT] [-r RSCRIPT]

options:
  -h, --help            show this help message and exit
  -iA1 INPUTA1, --inputA1 INPUTA1
                        The input COL1A1 peptide sequences fasta file to be used. Needs to contain taxonomic informations
  -iA2 INPUTA2, --inputA2 INPUTA2
                        The input COL1A2 peptide sequences fasta file to be used. Needs to contain taxonomic informations
  -o OUTPUT, --output OUTPUT
                        The output folder for where the theoretical peptide csvs will be saved.
  -r RSCRIPT, --Rscript RSCRIPT
                        The file path of Rscript.exe
```
