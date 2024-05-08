# Generating the LCMSMS Rules from the LCMSMS Mascot Data files
## Overview
This code uses the Liquid Chromatography Mass Spectrometry Mass Spectrometry (LCMSMS) data from multiple species for COL1 and groups the data by peptides that are in the same position in the sequence in different species.
The code uses these grouops to work out the likely number of oxidations and demidations for each group of peptides. The output can subsequently be used to filter the possible PTMs
in the theoretical peptides pipeline

## Usage
```
python peptide_rules.py
```
### Inputs
The inputs are a set of mascot CSV output files from LCMSMS data for species similar to the species trying to be identified (if searching for mammals there are already 7 mammalian csv files).
The columns required are: pep_seq, pep_score, pep_start, pep_end, pep_exp_mr, pep_miss, pep_var_mod and prot_acc. Column headers are needed. For an example:

