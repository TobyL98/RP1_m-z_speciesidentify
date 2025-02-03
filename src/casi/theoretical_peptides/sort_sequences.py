"""sort_sequences.py
Contains all the code to have as input  COL1A1 and
COL1A2 files for different species.. The code cleans the sequences to ensure
the sequences are correct and then merges the sequences.
Output is a csv file with a row for each sequence"""

import pandas as pd
import re
import sys


def cleana1(fileName):
    """reads in the fasta file taken from NCBI for COL1A1
    and outputs a clean fasta file this can then be used in
    fasta_seq_amendA1A2 based on sequences downloaded from NCBI
    can be used generally for Mammalia class but may miss
    some sequences based on exclusion rules.
    PLEASE CHECK RULES IF SPECIES MISSING"""
    fileObj = open(fileName, "r")
    sequences = {}  #  a dict, to contain all our sequences ...
    for line in fileObj:
        if line.startswith(">"):  # Ah ha ! a new sequence?
            name = line[1:].rstrip("\n")  # remove newline, ignore '>'

            sequences[name] = ""
        else:
            # all to uppercase
            line = line.upper()
            # replace '-' with X when residue unknown
            # otherwise will break downstream code
            line = line.replace("-", "X")
            sequences[name] += line.rstrip("\n")  # concat to growing dict value
    # print(len(sequences))

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
            value = value[start:end]

        # check sequence is correct length
        # if it is will add to dictionary
        if len(value) >= 1055 and len(value) <= 1058:
            clean_sequences[key] = value
    clean_num = len(clean_sequences)

    ## Converting back to fasta format
    fh = open("Sequences/COL1A1_seqs_clean_NCBI.fasta", "w")
    for key in clean_sequences:
        # first line is key
        print(">{0}".format(key), file=fh)
        # add lines for sequence
        print(clean_sequences[key], file=fh)

        fileObj.close()
    print("COL1A1 sequences have been cleaned")
    print("Number of clean COL1A1 sequences = {0}".format(clean_num))
    print("Output is COL1A1_seqs_clean_NCBI.fasta")
    print("######################################")


def cleana2(filename):
    """Reads in the fasta file taken from NCBI for COL1A1 and
    outputs a clean fasta file this can then be used in
    fasta_seq_amendA1A2.py based on sequences downloaded from NCBI
    can be used generally for Mammalia class but may miss some
    sequences based on exclusion rules.
    PLEASE CHECK RULES IF SPECIES MISSING"""
    fileobj = open(filename, "r")
    sequences = {}  #  a dict, to contain all our sequences ...
    for line in fileobj:
        if line.startswith(">"):  # Ah ha ! a new sequence?
            name = line[1:].rstrip("\n")  # remove newline, ignore '>'

            sequences[name] = ""
        else:
            # all to uppercase
            line = line.upper()
            # replace '-' with X when residue unknown
            # otherwise will break downstream code
            line = line.replace("-", "X")
            sequences[name] += line.rstrip("\n")  # concat to growing dict value

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
        if match1 and match2:
            value = value[start:end]
        # print(key)
        # print(len(value))
        # print(value)

        # check sequence is correct length
        # if it is will add to dictionary
        if len(value) >= 1038 and len(value) <= 1042:
            clean_sequences[key] = value
    clean_num = len(clean_sequences)

    ## Converting back to fasta format
    fh = open("Sequences/COL1A2_seqs_clean_NCBI.fasta", "w")
    for key in clean_sequences:
        # first line is key
        print(">{0}".format(key), file=fh)
        # add lines for sequence
        print(clean_sequences[key], file=fh)

    fileObj.close()
    print("COL1A2 sequences have been cleaned")
    print("Number of clean COL1A2 sequences = {0}".format(clean_num))
    print("Output is COL1A2_seqs_clean_NCBI.fasta")
    print("######################################")


# function reads the fasta file and outputs as dictionary
# also standardises species information
def readCOLFastaFile(fileName):
    fileObj = open(fileName, "r")
    sequences = {}  #  a dict, to contain all our sequences ...
    names = []  #  ... and a list to contain their names
    for line in fileObj:
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

    fileObj.close()
    return sequences


# function combines the COL1A1 and COL1A2 sequences to one sequence
def COLA1A2combine():
    # formats the A1 sequences
    file = "Sequences/COL1A1_seqs_clean_NCBI.fasta"
    output_COL1A1_dict = readCOLFastaFile(file)

    # formats COL1A2 sequence
    file2 = "Sequences/COL1A2_seqs_clean_NCBI.fasta"
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
    fh = open("Sequences/COL1A1A2_combined_seqs_NCBI.fasta", "w")
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


def FastaToCSV():
    """Reads in the fasta file (in this case combined COL1A1
    and COL1A2 file for NCBI data)
    Outputs a csv with columns for the taxonomic information
    and the sequence."""
    fileName = "Sequences/COL1A1A2_combined_seqs_NCBI.fasta"

    fileObj = open(fileName, "r")
    known_taxons = []  #  ... and a list to contain their names
    sequences = {"sequence": []}  # creates sequence dictionary
    sequence_count = 0

    for line in fileObj:
        if line.startswith(">"):  # Ah ha ! a new sequence?
            if sequence_count > 0:  # ensures first sequence created before append
                sequences["sequence"].append(
                    sequence
                )  # appends sequences to dictionary
            sequence = ""  # resets the sequence when we reach new entry
            sequence_count = 1

            header = line[1:].rstrip("\n")  # remove newline, ignore '>'

            ## gets taxonomic information format:
            # CLASS=Mammalia|SUBCLASS=Theria|INFRACLASS=Prototheria|ORDER=Monotremata
            # assign each taxonomic level to list in dictionary
            taxons = header.split("|")
            for info in taxons:
                taxon, name = info.split("=")[0:2]

                if (
                    taxon in known_taxons
                ):  # if taxon already a key appends new name to list
                    sequences[taxon].append(name)
                else:  # otherwise adds the name to list of known taxon key
                    sequences[taxon] = [name]
                    known_taxons.append(taxon)

        else:
            # adds lines to create sequence
            sequence += line.rstrip("\n")

    # ensure last sequence is appended
    sequences["sequence"].append(sequence)

    # convert dictionary to pandas dataframe
    df = pd.DataFrame(data=sequences)

    fileObj.close()
    df.to_csv("Sequences/sequences_taxon_NCBI.csv")
    print("Fasta file converted to CSV")
    print("Output is sequences_taxon_NCBI.csv")
    print("######################################")


if __name__ == "__main__":
    sys.exit()
