#####################
## fasta_seq_amendA1A2_NCBI.py
#####################

# reads in the cleaned COL1A1 and COL1A2 fasta file
# does this for both COL1A1 and COL1A2
# creates dictionaries where key is species
# and sequence is COL1A1 or COL1A2 sequence
# no taxonomic information as straight from NCBI
# combines the COL1A1 and COL1A2 dictionary values (sequences) by keys
# outputs as a fasta file

import sys
from pathlib import Path

import pandas as pd


# function reads the fasta file and outputs as dictionary
# also standardises species information
def readCOLFastaFile(file_name):
    file_obj = open(file_name, "r")
    sequences = {}  #  a dict, to contain all our sequences ...
    for line in file_obj:
        if line.startswith(">"):  # Ah ha ! a new sequence?
            header = line[1:].rstrip("\n")  # remove newline, ignore '>
            # gets out the sepcies name
            if header.endswith("]"):
                species_name = header.split("[")[1]
                species_name = species_name.strip("]")
                species_ID = True

                # assigns standard format for species name
                genus = species_name.split(" ")[0]
                name = "GENUS=" + genus
                name += "|SPECIES=" + species_name
            # if no species name enters false so not used further
            else:
                species_ID = False

            ## will only create new dictionary item if
            ## if there is a species name
            if species_ID == True:
                sequences[name] = ""
        # adds the sequences lines as the dictionary value
        elif species_ID == True:
            sequences[name] += line.rstrip("\n")  # concat to growing dict value
        else:
            continue

    file_obj.close()
    return sequences


# function combines the COL1A1 and COL1A2 sequences to one sequence
def COLA1A2combine(output_dir):
    # formats the A1 sequences
    file = output_dir / "COL1A1_seqs_clean_NCBI.fasta"
    output_COL1A1_dict = readCOLFastaFile(file)

    # formats COL1A2 sequence
    file2 = output_dir / "COL1A2_seqs_clean_NCBI.fasta"
    output_COL1A2_dict = readCOLFastaFile(file2)

    ### Combining COL1A1 and COL1A2 dictionaries
    # now one sequence with R (Arginine) in between
    # combine if they have same key (i.e., same species)
    A1A2_combined_dict = {}
    for key in output_COL1A1_dict:
        if key in output_COL1A2_dict:
            A1A2_combined_dict[key] = (
                output_COL1A1_dict[key] + "R" + output_COL1A2_dict[key]
            )

    ## Converting back to fasta format
    output_file = output_dir / "COL1A1A2_combined_seqs_NCBI.fasta"
    fh = open(output_file, "w")
    for key in A1A2_combined_dict:
        # first line is key
        print(">{0}".format(key), file=fh)
        # add lines for sequence
        print(A1A2_combined_dict[key], file=fh)

    fh.close()
    print(
        "COL1A1 and COL1A2 sequences have been combined and taxonomic information added."
    )
    print(
        "Number of complete COL1A1 and COL1A2 Sequences = {0}".format(
            len(A1A2_combined_dict)
        )
    )
    print("Output is COL1A1A2_combined_seqs_NCBI.fasta")
    print("######################################")


if __name__ == "__main__":
    sys.exit()
