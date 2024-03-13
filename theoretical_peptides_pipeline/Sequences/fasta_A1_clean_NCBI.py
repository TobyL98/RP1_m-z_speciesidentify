#####################
## fasta_A1_clean_NCBI.py
#####################

# reads in the fasta file taken from NCBI 
# for COL1A1
# and outputs a clean fasta file
# this can then be used in fasta_seq_amendA1A2.py
# based on sequences downloaded from NCBI
# can be used generally for Mammalia class
# but may miss some sequences based on exclusion rules
# PLEASE CHECK RULES IF SPECIES MISSING

import pandas as pd
import re

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

    # filters the A1 sequences
    # cleans the dataset
    clean_sequences = {}
    for key, value in sequences.items():
        # finds the correct start position
        match1 = re.search(r"Q([A-Z]{2})YG", value)
        if match1:
            start = match1.start()

        # finds the correct end position
        match2 = re.search(r"YRA", value)
        if match2:
            end = match2.end()
        
        # if sequences has both start and end position
        # will slice correct section of sequence
        if match1 and match2:
            value = value[start: end]
        
        # check sequence is correct length
        # if it is will add to dictionary
        if len(value) >= 1055 and len(value) <= 1058:
            clean_sequences[key] = value
    clean_num = len(clean_sequences)

    ## Converting back to fasta format
    fh = open("Sequences/COL1A1_seqs_clean_NCBI.fasta", "w")
    for key in clean_sequences:
        # first line is key
        print(">{0}".format(key), file = fh)
        # add lines for sequence
        print(clean_sequences[key], file = fh)

        fileObj.close()
    print("COL1A1 sequences have been cleaned")
    print("Number of clean COL1A1 sequences = {0}".format(clean_num))
    print("Output is COL1A1_seqs_clean_NCBI.fasta")
    print("######################################")
    


