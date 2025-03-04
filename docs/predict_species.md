# Identifies a Species by Comparing the ZooMS Peak List with The Theoretical Peptides

This guide gives instructions on how to predict the identity of the species from the ZooMS Peptide Mass Fingerprint peak list. This is done by comparing matches between the m/z values generated from the sequences for each species and the PMF peak list m/z values. The species with the most matches will be identified as the top species.

## Step 1 - Checks Before Starting

Before running the script, please ensure the package 'casi' is correctly installed in your environment. 
This can be checked in the command line with:
```
pip show casi
```
Alternatively use 'pip list' and see if 'casi' is one of the packages listed.

To check with conda:
```
conda list casi
```
Alternatively use 'conda list' and see if 'casi' is one of the packages listed.

## Input - Input Theoretical Peptides CSV Files (-it)

This is the theoretical peptides with m/z values that were generated using the 'theoretical_peps' script. Please see the 'generate_peptides.md' README for instructions on how to generate this folder. The folder required will be called 'filtered_peptides' within the output folder you specified in the 'theoretical_peps' script.

## Input - ZooMS Peptide Mass Fingerprint Peak List (-ip)

This a peak list .txt file that will have needed to be generated from the mzxml file. The peak list will have two columns with m/z values in the first column and intensity in the second column. There is an example located here: data/inputs/example_peaklist.txt. The start of the example is below:
```
785.368371779	22.9861030893
809.402735853	22.2748415047
810.375320918	18.7618098992
833.03300322	72.3902620525
836.406724125	405.027417201
837.407612611	289.907062763
```

## Optional Input - Threshold for a Match (-t)

The threshold for a match is the tolerance that counts as a match between the Peptide Mass Fingerprint (PMF) m/z value and the theoretical peptide m/z value. The default is +-0.2. This means that if the m/z value in the PMF is 1500.0 then the matching theoretical m/z value could be between 1499.8 and 1500.2.

## Optional Input - Top 5 Matches (-m5)

This has two options (1 or 0). If 1 is inputted, it will provide an xlsx (excel) file of m/z peak matches for the top 5 species in the output folder. Default is 0. This is useful if you want to interrogate the matches manually in more detail.

## Step 2 - Running the Script

The script 'compare_score.py' is used to compare matches between the m/z values generated from the sequences for each species and the PMF peak list to identify the species. There are three required inputs: the input species theoretical peptides csv files folder (-it), the ZooMS PMF peak list (-ip) and the output folder (-o). There are two optional inputs: threshold for a match (-t) and the option to output all of the matching m/z values for the top 5 matches (-m5).
The first example is without the optional input:
```
compare_score -ip example_peak_list.mzXML -it theoretical_results/filtered_peptides -o results_folder
```
The second example is with the optional inputs
```
compare_score -ip example_peak_list.mzXML -it theoretical_results/filtered_peptides -o results_folder -t 0.2 -m5 1
```
Any further instructions or detail required please look at:
```
compare_score --help
```

## Step 3 - Interpreting the Outputs

The species with the top 10 highest match scores are outputted to the command line. For more detail, the matches for every species are outputted to the results csv file.
For classification, if the top species has at least one greater match than the second top it is classified as that species. This is shown in the first example where *Bubalus bubalis* is correctly identified:
```
RESULTS:
|    | GENUS   | SPECIES                  |   Match |
|---:|:--------|:-------------------------|--------:|
|  0 | Bubalus | Bubalus bubalis          |     131 |
|  1 | Bubalus | Bubalus kerabau          |     130 |
|  2 | Bos     | Bos javanicus            |     123 |
|  3 | Bos     | Bos indicus x Bos taurus |     123 |
|  4 | Bos     | Bos taurus               |     123 |
|  5 | Bison   | Bison bison bison        |     123 |
|  6 | Moschus | Moschus berezovskii      |     122 |
|  7 | Bos     | Bos indicus              |     122 |
|  8 | Bos     | Bos mutus                |     122 |
|  9 | Dama    | Dama dama                |     118 |
```
If there are multiple species that are the top match and have the same number of matches, then it is matched to the taxonomic level that all those species would be included in. For example, these are all matched to genus *Panthera*
```
RESULTS:
|    | GENUS        | SPECIES                  |   Match |
|---:|:-------------|:-------------------------|--------:|
|  0 | Panthera     | Panthera leo             |      76 |
|  1 | Panthera     | Panthera onca            |      76 |
|  2 | Panthera     | Panthera uncia           |      76 |
|  3 | Felis        | Felis catus              |      74 |
|  4 | Leopardus    | Leopardus geoffroyi      |      73 |
|  5 | Prionailurus | Prionailurus viverrinus  |      73 |
|  6 | Lynx         | Lynx rufus               |      73 |
|  7 | Lynx         | Lynx canadensis          |      73 |
|  8 | Neofelis     | Neofelis nebulosa        |      73 |
|  9 | Prionailurus | Prionailurus bengalensis |      73 |
```
Finally, these can all be matched to the family *Mustelidae*
```
RESULTS:
|    | GENUS       | SPECIES               |   Match |
|---:|:------------|:----------------------|--------:|
|  0 | Mustela     | Mustela lutreola      |      74 |
|  1 | Mustela     | Mustela nigripes      |      74 |
|  2 | Neogale     | Neogale vison         |      74 |
|  3 | Lutra       | Lutra lutra           |      74 |
|  4 | Mustela     | Mustela erminea       |      74 |
|  5 | Mustela     | Mustela putorius furo |      74 |
|  6 | Meles       | Meles meles           |      72 |
|  7 | Lontra      | Lontra canadensis     |      71 |
|  8 | Eumetopias  | Eumetopias jubatus    |      70 |
|  9 | Callorhinus | Callorhinus ursinus   |      70 |
```
