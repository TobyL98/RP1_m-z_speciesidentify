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
- To run the R script the correct file location for Rscript.exe is required. If you do not know where this is see the 'Known issues' section

### Input Format
The input FASTA files should have the following format
```
>CAH6775934.1 Col1a1 [Phodopus roborovskii]
MFSFVDLRLLLLLGATALLTHGQEDIPEITCIHDGVKIPNGETWKADTCLICICHNGTAVCDAVVCQEKL
DCVNPQTREGECCPFCPEEFVSPDQELIGVEGPKGDRGPQGPRGPAGPPGKDGIPGQPGLPGPPGPPGPP
```
The species name needs to be at the end of the header in square brackets. This is the default output from NCBI protein

### Output
The output is a csv file for each species with the generated peptides and their respective mass values. By default the output csv files are saved to PTM_rules/Integration_code/results_NCBI.
An example of the start in excel is below.

<img width="625" alt="theor_out_example" src="https://github.com/TobyL98/RP1_m-z_speciesidentify/assets/158182593/a758da8c-da45-464d-b657-df312d919189">

### Changing LCMSMS Data
A key part of this code is to filter all the possible Post-translational modifications (PTMs) with results from LCMSMS analysis to only include the possible PTMs from the LCMSMS data in the final results.
The LCMSMS data is the results from a Mascot query of the LCMSMS experiment. The current LCMSMS results used are in PTM_rules/LCMSMS. They include the mammalian species: a bat species from genus Myotis, Rattus norvegicus, Felis catus, a species from Bos genus, Canis lupus familiaris, a species from Elephantidae family and Mus musculus.

If you're doing analysis for a different class (e.g., birds) or very divergent/rare mammalian species then you can change the LCMSMS data by removing the mascot result csv files and replacing them with new csv files. You will potentially need several different species to obtain all the possible PTMs. Once the LCMSMS files have been replaced, run the code *peptide_rules.py* to integrate them into the pipeline

## Common errors
**1. Rscript.exe file cannot be found**

One section of the code uses an R script modified from the package 'Baccolite' (Hickinbotham et al., 2020) to generate the theoretical peptides and their corresponding masses.
To run this in python using subprocess requires the Rscript.exe which runs R. This is entered in the comman line option -r. 
If you do not know where this file is you can find it:
- On the command line enter the command "where /r "c:\Program Files" Rscript.exe". Then pick the correct Rscript.exe location (if you did not install R into Program files then search in the folder you installed it into)
- In R studio enter the command: "file.path(R.home(), "bin", "Rscript.exe")"

## Code that makes up the pipeline
There are several functions that make up this pipeline which are briefly explained. For more information look at the individual code:

**1. fasta_A1_clean_NCBI.py** - cleans the fasta format COL1A1 sequences downloaded from a database. Output is a fasta file called *COL1A1_seqs_clean_NCBI.fasta* in the 'Sequences' folder

**2. fasta_A2_clean_NCBI.py** - cleans the fasta format COL1A2 sequences downloaded from a database. Output is a fasta file called *COL1A2_seqs_clean_NCBI.fasta* in the 'Sequences' folder

**3. fasta_seq_amendA1A2_NCBI.py** - merges the COL1A1 and COL1A2 sequences into one COL1 sequence and reformats the fasta header. Output is a fasta file called *COL1A1A2_combined_seqs_NCBI.fasta* in the 'Sequences' folder

**4. fasta_to_csv_NCBI.py** - reads in the previous COL1 fasta file and outputs as csv file with a column for genus name, species name and the COL1 sequence for each species. Output is called *sequences_taxon_NCBI.csv* in the 'Sequences' folder

**5. run_parse_seq_NCBI.R** - The fifth code function is in R and is based off parse_seq.R from the package "Baccolite" (Hickinbotham et al., 2020). For each species COL1 sequence, this code simulates trypsin digestion to produce the peptides and then calculates the masses of the peptides. Within the masses it has to account for the common COL1 PTMs. These common PTMs are oxidation (+16) and deamidation (+1). It has to account for all the combinations of PTMs. For example, if there were three possible oxidations it would account for 4 masses (0 oxidations, 1 oxidation (+16), 2 oxidations (+32) and 3 oxidations (+48). The outputs are csv files for each species of the peptide sequence, PTMs and masses. The outputs are in the folder "insilico_digest/in_silico_res_NCBI"

**6. integrate_NCBI.py** - This code reads in the previous csv files and the LCMSMS Mascot result files and uses the known PTMs (oxidations and deamidations) from the LCMSMS data to filter the possible PTMs to reduce the number of masses in the output. The output is a filtered csv file for each species with the peptide sequences, PTMs and masses. This final output is outputted to the folder specified in -o in the command line. 



## References
Simon Hickinbotham, Sarah Fiddyment, Timothy L Stinson, Matthew J Collins (2020) How to Get Your Goat: Automated Identification of Species from MALDI-ToF Spectra Bioinformatics, March 2020


