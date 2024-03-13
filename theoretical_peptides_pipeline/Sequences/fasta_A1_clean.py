#####################
## fasta_A1_clean.py
#####################

# reads in the fasta file taken from e.g, ensembl or NCBI entrez
# for COL1A1
# and outputs a clean fasta file
# this can then be used in fasta_seq_amendA1A2.py
# based on sequences Sakim took from ensembl
# NOT A GENERAL CLEANING TOOL straight from ensembl

import pandas as pd

def cleanA1(fileName):
    fileObj = open(fileName, 'r')
    sequences = {}   #  a dict, to contain all our sequences ...
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            name = line[1:].rstrip('\n')  # remove newline, ignore '>'
            
            sequences[name] = ''
        else:
            # all to uppercase
            line = line.upper()
            # replace '-' with X when residue unknown
            # otherwise will break downstream code
            line = line.replace("-", "X")
            sequences[name] += line.rstrip('\n')  # concat to growing dict value
    #print(len(sequences))

    # filters the A1 sequences by sequence length
# cleans the dataset
    clean_sequences = {}
    for key, value in sequences.items():
        if len(value) >= 1055 and len(value) <= 1058:
            clean_sequences[key] = value
    clean_num = len(clean_sequences)

    # Converting back to fasta format
    fh = open("Sequences/COL1A1_seqs_clean.fasta", "w")
    for key in clean_sequences:
        # first line is key
        print(">{0}".format(key), file = fh)
        # add lines for sequence
        print(clean_sequences[key], file = fh)

        fileObj.close()
    print("COL1A1 sequences have been cleaned")
    print("Number of clean COL1A1 sequences = {0}".format(clean_num))
    print("Output is COL1A1_seqs_clean.fasta")
    print("######################################")
    

#file = 'COL1A1_seqs.fasta'
#cleanA1(file)

