# Sequence Based Species Identification of ZooMS Peptide Mass Fingerprints 

## Overview
This project uses Collagen 1 sequence data downloaded from databases (e.g., NCBI protein) to semi-automate the identification process of Zoo Archaelology by Mass Spectrometry (ZooMS) Peptide Mass Fingerprint (PMF) data. In order to do this there are two pipelines for this process:
1. Theoretical peptide m/z values for species are generated from COL1α1 and COL1α2 sequences downloaded as FASTA files from a database. A README describes the process within the folder
2. The experimental m/z value file (created in software such as mMass) is compared to the theoretical peptide m/z values for matches.  A README describes the process within the folder

To briefly explain the pipelines, the COL1α1 and COL1α2 sequences are downloaded as a fasta file from a database for a group of species (e.g., mammals). The COL1α1 and COL1α2 sequences are cleaned and merged to form COL1. The COL1 sequences undergo *in silico* trypsin digest to form peptides allowing for up to one missed cleavage. The m/z values for the peptides are calculated based off their contistuent amino acids masses. These m/z values need to account for all combinations of post-translational modifications (PTMs). These PTMs are filtered based on the likely PTMs from previous liquid chromatography mass spectrometry mass spectrometry ZooMS Mascot results to reduce the number of possible m/z values. The output for each species is a csv file containing the species name and the m/z values for all the peptides. Peak picking from the ZooMS experimental PMF leads to a text file with experimental m/z values (this needs to be done with external software). The experimental m/z values and the species theoretical m/z values are compared for matches within a certain threshold (e.g., +- 0.2). The species with the most matches is the most closely related species to the sample.



![flow_diagram_3](https://github.com/TobyL98/RP1_m-z_speciesidentify/assets/158182593/fe4de66b-4cf2-497b-b4f6-657ec5526320)

## Installation
Currently still in development and not published as a package on PyPi. Therefore, currently the repository has to be cloned before installing the package.

**Clone this repository**
```
git clone https://github.com/TobyL98/RP1_m-z_speciesidentify.git
```
If you do not have git installed you can instead in the github repository select 'Code' and 'Download Zip'. Then save and unzip the foldr in the location you choose.

**Option 1 - Pip**

As the package is not currently published you must use pip within the directory you now have called 'RP1_m-z_speciesidentify'.
```
pip install casi
```
**Option 2 - Conda**

Unfortunately the package is not available in conda. However a conda environment with the correct dependencies is installed using the yaml file. You will still need to use pip as recommended above to install the package.
```
conda env create -f  species_identify.yml
```

## Further Instructions

For further instructions on how to generate the theoretical peptide m/z values from COL1 peptide sequences and to use these to predict the identity of species from its ZooMS fingerprint look at:
**docs/generate_peptides.md**  
and  
**docs/predict_species.md**

## Contributions
Contributions are welcomed for this repository!

Please raise any issues that you encounter



