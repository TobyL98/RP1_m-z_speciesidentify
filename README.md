# Sequence Based Species Identification of ZooMS Peptide Mass Fingerprints 

## Overview
This project uses Collagen 1 sequence data downloaded from databases (e.g., NCBI protein) to semi-automate the identification process of Zoo Archaelology by Mass Spectrometry (ZooMS) Peptide Mass Fingerprint (PMF) data. In order to do this there are two pipeline for this process:
1. Theoretical peptide m/z values for species are generated from COL1α1 and COL1α2 sequences downloaded as FASTA files from a database. This pipeline is run from *theoretical_peps_NCBI.py* within theoretical_peptides_folder. A README describes the process within the folder
2. The experimental m/z value file (created in software such as mMass) is compared to the theoretical peptide m/z values for matches. This pipeline is run from *compare_NCBI.py* within Compare_results folder. A README describes the process within the folder

To briefly explain the pipelines, the COL1α1 and COL1α2 sequences are downloaded as a fasta file from a database for a group of species (e.g., mammals). The COL1α1 and COL1α2 sequences are cleaned and merged to form COL1. The COL1 sequences undergo *in silico* trypsin digest to form peptides allowing for up to one missed cleavage. The m/z values for the peptides are calculated based of their contistuent amino acids masses. These m/z values need to account for all combinations of post-translational modifications (PTMs). These PTMs are filtered based on the likely PTMs from previous liquid chromatography mass spectrometry mass spectrometry ZooMS Mascot results to reduce the number of possible m/z values. The output for each species is a csv file containing the species name and the m/z values for all the peptides. Peak picking from the ZooMS experimental PMF leads to a text file with experimental m/z values (this is done in external software such as mMass). The experimental m/z values and the species theoretical m/z values are compared for matches within a certain threshold (e.g., +- 0.2). The species with the most matches is the most closely related species to the sample



![flow_diagram_3](https://github.com/TobyL98/RP1_m-z_speciesidentify/assets/158182593/fe4de66b-4cf2-497b-b4f6-657ec5526320)

## Getting Started - Users
1. Clone this repository
2. Download all the required dependencies
3. To generate theoretical peptide spectra from species COL1 sequences look at README in theoretical_peptides_pipeline folder
4. To compare the theoretical and experimental peptide spectra follow the README in Compare_results

## Getting Started - Developers
Please clone this repository and use it as you would like. The information and code for the theoretical peptide spectra generation are in the theoretical_peptides_pipeline_folder. 
The information and code for the theoretical and experimental peptide spectra comparison is in the Compare_results folder

## Dependencies
**Python**
- using version 3.12.1
- Pandas v 2.1.4
- Numpy v 1.26.3

**R**
- Using version 4.3.2
- dplyr v 1.1.4
- stringr v 1.5.1
- Bacollite v 1.0.1 (https://github.com/bioarch-sjh/bacollite/blob/master/INSTALL.md) (Hickinbotham et al., 2021)

# References

Simon Hickinbotham, Sarah Fiddyment, Timothy L Stinson, Matthew J Collins (2020) How to Get Your Goat: Automated Identification of Species from MALDI-ToF Spectra Bioinformatics, March 2020


