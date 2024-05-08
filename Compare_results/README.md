# Comparing the experimental ZooMS m/z values to the species theoretical m/z values
## Overview
This code compares the m/z values from an experimental ZooMS spectra to the theoretical m/z values that have been generated for a variety of species in the theoretical_peptides_pipeline.
There is a match if the m/z value from the experimental spectra is within a certain threshold (e.g., +- 0.2) of the theoretical m/z value. The species with the most matches is most closely
related to the experimental spectra.

## Usage
```
usage: compare_NCBI.py [-h] [-it INPUTTHEOR] [-o OUTPUT] [-ip INPUTPMF] [-t THRESHOLD] [-m5 TOP5]

This code compares the experimental PMF peaks from a sample and the theoretical peaks generated from the all the species COL1 theoretical peptides generated in   
the theoretical peptides pipeline. It will score a match if the theoretical peptide m/z valus is within a specified threshold (+- 0.2 is default). The output to  
the command line will be the top 10 match scores (number of matches) species. All species match scores wil be outputted to a CSV file.

options:
  -h, --help            show this help message and exit
  -it INPUTTHEOR, --inputTheor INPUTTHEOR
                        the folder that contains the theoretical peptides csv files to compare against PMF
  -o OUTPUT, --output OUTPUT
                        The file name for the output results file of number of matches
  -ip INPUTPMF, --inputPMF INPUTPMF
                        The input Peptide mass fingerprint (PMF) from an unknown organism.
  -t THRESHOLD, --threshold THRESHOLD
                        The threshold for matches between the experimental and theoretical spectrum. Default is 0.2 Da
  -m5 TOP5, --top5 TOP5
                        If 1 is inputted will provide an excel file of m/z peak matches for the top 5 match counts as an xlsx (excel) file. Default is 0
```

### Example usage
*Default. Still need to include the input theoretical PMF*
```
python compare_NCBI.py -ip PMF_samples/Rodentia/Rattus_rattus_sample.txt
```

*All input and output options included*
```
python compare_NCBI.py -ip PMF_samples/Rodentia/Rattus_rattus_sample.txt -it "C:\Users\tobyl\OneDrive - The University of Manchester\Bioinformatics Masters\Research project 1\Git_repositories\RP1_m-z_speciesidentify\theoretical_peptides_pipeline\PTM_rules\Integration_code\results_NCBI" -o Results/rattus_rattus_matches.csv -t 0.2 -m5 1  
```

### Inputs
