# Sequence Based Species Identification of ZooMS Peptide Mass Fingerprints 

## Overview
This project uses Collagen 1 sequence data downloaded from databases (e.g., NCBI protein) to semi-automate the identification process of Zoo Archaelology by Mass Spectrometry (ZooMS) Peptide Mass Fingerprint (PMF) data. In order to do this there are two pipeline for this process:
1. Theoretical peptide m/z values for species are generated from COL1α1 and COL1α2 sequences downloaded as FASTA files from a database. This pipeline is run from *theoretical_peps_NCBI.py* within theoretical_peptides_folder. A README describes the process within the folder
2. The experimental m/z value file (created in software such as mmass) is compared to the theoretical peptide m/z values for matches. This pipeline is run from *compare_NCBI.py* within Compare_results folder. A README describes the process within the folder

Briefly, the COL1α1 and COL1α2 sequences are downloaded as a fasta file from a database for a group of species (e.g., mammals). The COL1α1 and COL1α2 sequences are cleaned and merged to form COL1. The COL1 sequences undergo *in silico* trypsin digest to form peptides allowing for up to one missed cleavage. The m/z values for the peptides are calculated based of their contistuent amino acids masses. These masses need to account for all combinations of Post-translationa modifications (PTMs).


