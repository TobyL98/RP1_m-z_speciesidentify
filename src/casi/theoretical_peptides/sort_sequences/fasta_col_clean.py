"""fasta_a1_clean
Reads in the fasta file for COL1A1/COL1A2 and outputs a clean fasta file.

Usage: Used on fasta files of COL1A1/COl1A2 sequences in NCBI format.
Can be used generally for Mammalia/Aves class but may miss some
sequences based on exclusion rules.
PLEASE CHECK RULES IF SPECIES MISSING

The file can also be imported as a module and contains the following functions:

    * read_fasta - reads the fasta file and converts to dictionary.
    * clean_a1 - Cleans the COL1A1 sequences
    * clean_a2_birds - Cleans the COl1A2 sequences for birds
    * clean_a2_mammals - Cleans the COl1A2 sequences for mammals
    * convert_fasta - Converts to fasta format and saves the file
    * run_clean_col - Runs the full pipeline cleaning the cola1 sequences in the fasta file.
"""

import re
import sys
from pathlib import Path

def read_fasta(file_name):
    """reads the fasta file and converts to dictionary.

    Parameters
    ----------
    file_name : str
        The location of the fasta file

    Returns
    -------
    sequences : dict
        Dictionary with fasta definition as name and sequence
        as the value
    """
    file_obj = open(file_name, "r", encoding="utf-8")
    sequences = {}
    for line in file_obj:
        if line.startswith(">"):  # definition for new sequence
            name = line[1:].rstrip("\n")
            sequences[name] = ""
        else:
            line = line.upper()
            # replace '-' with X when residue unknown
            line = line.replace("-", "X")
            sequences[name] += line.rstrip("\n")
    return sequences


def clean_a1(sequence_dict):
    """Cleans the COL1A1 sequences.

    Removes any sequences which are too short, too long or inaccurate.
    And cleaves to get the final protein product.

    Parameters
    ----------
    sequence_dict : dict
        The dictionary of sequence names(key) and sequence(value)

    Returns
    -------
    clean_sequences : dict
        Dictionary of cleaned sequences with names(key) and sequences(value)
    """
    clean_sequences = {}
    for key, value in sequence_dict.items():
        # finds the correct start position
        match1 = re.search(r"Q([A-Z]{2})YG", value)
        if match1:
            start = match1.start()
        # finds the correct end position
        match2 = re.search(r"YRA", value)
        if match2:
            end = match2.end()
        # if sequences has both start and end position
        if match1 and match2:
            value = value[start:end]
        # check sequence is correct length
        if len(value) >= 1055 and len(value) <= 1058:
            clean_sequences[key] = value
    return clean_sequences


def clean_a2_birds(sequence_dict):
    """Cleans the COL1A2 sequences.

    Removes any sequences which are too short, too long or inaccurate.
    And cleaves to get the final protein product. Cleavage rules for COl1A2
    in mammals is different to birds. This is for birds.

    Parameters
    ----------
    sequence_dict : dict
        The dictionary of sequence names(key) and sequence(value)

    Returns
    -------
    clean_sequences : dict
        Dictionary of cleaned sequences with names(key) and sequences(value)
    """
    clean_sequences = {}
    for key, value in sequence_dict.items():
        # finds the correct start position
        match1 = re.search(r"QYD[PG]S[K]", value)
        if match1:
            start = match1.start()
        # finds the correct end position
        match2 = re.search(r"[FY]RA", value)
        if match2:
            end = match2.end()
        # if sequences has both start and end position
        if match1 and match2:
            value = value[start:end]
        # check sequence is correct length
        if len(value) >= 1038 and len(value) <= 1042:
            clean_sequences[key] = value
    return clean_sequences

def clean_a2_mammals(sequence_dict):
    """Cleans the COL1A2 sequences.

    Removes any sequences which are too short, too long or inaccurate.
    And cleaves to get the final protein product. Cleavage rules in COl1A2
    for mammals is different to birds. This is for mammals.

    Parameters
    ----------
    sequence_dict : dict
        The dictionary of sequence names(key) and sequence(value)

    Returns
    -------
    clean_sequences : dict
        Dictionary of cleaned sequences with names(key) and sequences(value)
    """
    clean_sequences = {}
    for key, value in sequence_dict.items():
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
        match3 = re.search(r"[FY]RA", value)
        if match3:
            end = match3.end()
        # if sequences has both start and end position
        if match1 and match3:
            value = value[start: end]
        elif match2 and match3:
            value = value[start: end]
        # check sequence is correct length
        if len(value) >= 1038 and len(value) <= 1042:
            clean_sequences[key] = value
    return clean_sequences


def convert_fasta(clean_sequences, output_dir, col_type):
    """Converts to fasta format and saves the file

    Parameters
    ----------
    clean_sequence_dict : dict
        Dictionary of cleaned sequences with names(key) and sequences(value)
    output_dir : str
        filepath for directory where the output file will be saved
    collagen_type : str
        whether it is 'COL1A1' or 'COL1A2' collagen

    Returns
    -------
    col1a1_seq_clean.fasta : file
        Fasta file of the cleaned sequences
    """
    output_filepath = Path(output_dir / f"{col_type}_seqs_clean_NCBI.fasta")
    fh = open(output_filepath, "w", encoding="utf-8")
    for key, value in clean_sequences.items():
        print(f">{key}", file=fh)
        print(value, file=fh)

    fh.close()

def print_outputs(clean_sequences, col_type):
    """Prints the outputs to the command line.
    
    Parameters
    ----------
    clean_sequences : dict
        Dictionary of cleaned seqeunces with names(key) and sequences(value)
    col_type : str
        Whether it is 'COL1A1' or 'COL1A2'
    """
    print(f"{col_type} sequences have been cleaned")
    clean_num = len(clean_sequences)
    print(f"Number of clean {col_type} sequences = {clean_num}")
    print(f"Output is {col_type}_seqs_clean_NCBI.fasta")
    print("######################################")

def run_clean_col(file_name, output_dir, collagen_type, class_input=None):
    """Runs the full pipeline cleaning the cola1/a2 sequences in the fasta file.

    Saves the fasta file in the output directory

    Parameters
    ----------
    file_name : str
        filepath for location of fasta file
    output_dir : str
        filepath for directory where the output file will be saved
    collagen_type : str
        whether it is 'COL1A1' or 'COL1A2' collagen
    class_input : str
        whether the COl1A2 sequences are from mammals or birds
    """
    sequence_dict = read_fasta(file_name)
    if collagen_type == "COL1A1":
        clean_sequence_dict = clean_a1(sequence_dict)
    if collagen_type == "COL1A2":
        if class_input == "birds":
            clean_sequence_dict = clean_a2_birds(sequence_dict)
        elif class_input == "mammals":
            clean_sequence_dict = clean_a2_mammals(sequence_dict)
    convert_fasta(clean_sequence_dict, output_dir, collagen_type)
    print_outputs(clean_sequence_dict, collagen_type)


if __name__ == "__main__":
    sys.exit(run_clean_col(sys.argv[1:]))
