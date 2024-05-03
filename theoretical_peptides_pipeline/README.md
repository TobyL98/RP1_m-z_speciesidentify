# Generating Theoretical Peptide Spectra from Sequences
## Overview

This pipeline generates theoretical collagen 1 (COL1) PMF spectra from COL1A1 and COL1A2 sequences in a database. The theoretical COL1 PMF spectra can then be compared to experimentally generated spectra

## Usage
```
usage: theoretical_peps_NCBI.py [-h] [-iA1 INPUTA1] [-iA2 INPUTA2] [-o OUTPUT] [-r RSCRIPT]

options:
  -h, --help            show this help message and exit
  -iA1 INPUTA1, --inputA1 INPUTA1
                        The input COL1A1 peptide sequences fasta file to be used. Needs to contain taxonomic informations.
                        Default is 'Sequences/COL1A1_seqs_NCBI.fasta'
  -iA2 INPUTA2, --inputA2 INPUTA2
                        The input COL1A2 peptide sequences fasta file to be used. Needs to contain taxonomic informations
                        Default is 'Sequences/COL1A2_seqs_NCBI.fasta'
  -o OUTPUT, --output OUTPUT
                        The output folder for where the theoretical peptide csvs will be saved.
                        Default is 'PTM_rules/Integration_code/results_NCBI'
  -r RSCRIPT, --Rscript RSCRIPT
                        The file path of Rscript.exe. Default is 'C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe'
```

### Example usage
*Running with default inputs and output*. 

This will generate theoretical peptide spectra for 202 different mammalian species
```
python theoretical_peps_NCBI.py 
```

*Running with different inputs and outputs*
```
python theoretical_peps_NCBI.py -iA1 Sequences/COL1A1_seqs_NCBI.fasta -iA2 Sequences/COL1A2_seqs_NCBI.fasta -o PTM_rules/Integration_code/results_NCBI -r "C:/Program Files/R/R-4.3.2/bin/x64/Rscript.exe" 
```
- The fasta input files iA1 and iA2 must have the correct format. See input section for example
- To run the R script the correct file location for Rscript.exe is required. If you do not know this see Known issues' section

### Input Format
The input FASTA files should have the following format
```
>CAH6775934.1 Col1a1 [Phodopus roborovskii]
MFSFVDLRLLLLLGATALLTHGQEDIPEITCIHDGVKIPNGETWKADTCLICICHNGTAVCDAVVCQEKL
DCVNPQTREGECCPFCPEEFVSPDQELIGVEGPKGDRGPQGPRGPAGPPGKDGIPGQPGLPGPPGPPGPP
```
The species name needs to be at the end of the header in square brackets. This is the default output from NCBI protein

### Output
The output is a csv file for each species with the generated peptides and their respective mass values. 
An example of the start in excel is below. Only the columns mass1, GENUS and SPECIES are required for comparing to experimental spectra

<img width="625" alt="theor_out_example" src="https://github.com/TobyL98/RP1_m-z_speciesidentify/assets/158182593/a758da8c-da45-464d-b657-df312d919189">

### Changing LCMSMS Data
A key part of this code is to filter all the possible Post-translational modifications (PTMs) with results from LCMSMS analysis to only include the possible PTMs in the LCMSMS data.
The LCMSMS data is the results from a Mascot query of the LCMSMS data. The current LCMSMS results used are in PTM_rules/LCMSMS. They include the mammalian species: a bat species from genus Myotis, Rattus norvegicus, Felis catus, a species from Bos genus, Canis lupus familiaris, a species from Elephantidae family and Mus musculus.

If you're doing analysis for a different class (e.g., birds) or very divergent/rare mammalian species then you can change the LCMSMS data by removing the mascot result csv files and replacing them with new csv files. You will potentially need several different species to obtain all the possible PTMS. Once, the LCMSMS files have been replaces run the code *peptide_rules.py* to integrate them into the pipeline

## Common errors
**1. Rscript.exe file cannot be found**

One section of the code using an R script modified from the package 'Baccolite' (Hickinbotham et al., 2020) to generate the theoretical peptides and their corresponding masses.
To run this in python using subprocess requires the Rscript.exe which runs R. This is entered in the comman line option -r. 
If you do not know where this file is you can find it using R studio:
- start



