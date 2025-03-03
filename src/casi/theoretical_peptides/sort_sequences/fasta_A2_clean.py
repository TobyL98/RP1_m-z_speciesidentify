#####################
## fasta_A2_clean.py
#####################

# reads in the fasta file taken from NCBI 
# for COL1A1
# and outputs a clean fasta file
# this can then be used in fasta_seq_amendA1A2.py
# based on sequences downloaded from NCBI
# can be used generally for Mammalia class
# but may miss some sequences based on exclusion rules
# PLEASE CHECK RULES IF SPECIES MISSING

import re
import sys
from pathlib import Path

import pandas as pd

def cleanA2(file_name, output_dir, species_class):
    file_obj = open(file_name, 'r', encoding="utf-8")
    sequences = {}   #  a dict, to contain all our sequences ...
    for line in file_obj:
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

    if species_class == "mammals": # for mammals
        # filters the A2 sequences
        # cleans the dataset
        clean_sequences = {}
        for key, value in sequences.items():
            # finds the correct start position
            # was getting some false positive matches with match2 in some cases (mouse as example)
            # so do match1 first to avoid this
            match1 = re.search(r"Q[YF][DSLA][A-Z][KG]", value)
            match2 = re.search(r"Q[YF][DSLA][A-Z][KGS]", value)
            if match1:
                start = match1.start()
            elif match2:
                start = match2.start()


            # finds the correct end position
            match3 = re.search(r"YRA", value)
            if match3:
                end = match3.end()


            # if sequences has both start and end position
            # will slice correct section of sequence
            if match1 and match3:
                value = value[start: end]

            elif match2 and match3:
                value = value[start: end]

            #check sequence is correct length
            # if it is will add to dictionary
            if len(value) >= 1038 and len(value) <= 1042:
               clean_sequences[key] = value
        clean_num = len(clean_sequences)

    elif species_class == "birds": # for birds
        # filters the A2 sequences
        # cleans the dataset
        clean_sequences = {}
        for key, value in sequences.items():
            # finds the correct start position
            # was getting some false positive matches with match2 in some cases (mouse as example)
            # so do match1 first to avoid this
            match1 = re.search(r"QYD[PG]S[K]", value)
            if match1:
                start = match1.start()

            # finds the correct end position
            match2 = re.search(r"[FY]RA", value)
            if match2:
                end = match2.end()


            # if sequences has both start and end position
            # will slice correct section of sequence
            if match1 and match2:
                value = value[start: end]

            #check sequence is correct length
            # if it is will add to dictionary
            if len(value) >= 1038 and len(value) <= 1042:
               clean_sequences[key] = value
        clean_num = len(clean_sequences)   

    ## Converting back to fasta format
    output_filepath = Path(output_dir / "COL1A2_seqs_clean_NCBI.fasta")
    fh = open(output_filepath, "w", encoding="utf-8")
    for key in clean_sequences:
        # first line is key
        print(">{0}".format(key), file = fh)
        # add lines for sequence
        print(clean_sequences[key], file = fh)

    file_obj.close()
    print("COL1A2 sequences have been cleaned")
    print("Number of clean COL1A2 sequences = {0}".format(clean_num))
    print("Output is COL1A2_seqs_clean_NCBI.fasta")
    print("######################################")

if __name__ == "__main__":
    sys.exit(cleanA2(sys.argv[1:]))
    


