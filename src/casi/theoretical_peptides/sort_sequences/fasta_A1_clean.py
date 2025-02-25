"""fasta_a1_clean
Reads in the fasta file for COL1A1 and outputs a clean fasta file.

Usage: Used on fasta files of col1a1 sequences downloaded from 
NCBI. Can be used generally for Mammalia class but may miss some 
sequences based on exclusion rules. 
PLEASE CHECK RULES IF SPECIES MISSING

The file can also be imported as a module and contains the following functions:

    * read_fasta - reads the fasta file and converts to dictionary.
    * clean_a1 - Cleans the COL1A1 sequences
    * convert_fasta - Converts to fasta format and saves the file
    * run_clean_a1 - Runs the full pipeline cleaning the cola1 sequences in the fasta file.
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
    file_obj = open(file_name, 'r', encoding="utf-8")
    sequences = {} 
    for line in file_obj:
        if line.startswith('>'): # definition for new sequence  
            name = line[1:].rstrip('\n')         
            sequences[name] = ''
        else:
            line = line.upper()
            # replace '-' with X when residue unknown
            line = line.replace("-", "X")
            sequences[name] += line.rstrip('\n')
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

def convert_fasta(clean_sequences, output_dir):
    """Converts to fasta format and saves the file
    
    Parameters
    ----------
    clean_sequence_dict : dict
        Dictionary of cleaned sequences with names(key) and sequences(value)
        
    Returns
    -------
    col1a1_seq_clean.fasta : file
        Fasta file of the cleaned sequences
    """
    output_filepath = Path(output_dir / "COL1A1_seqs_clean_NCBI.fasta")
    fh = open(output_filepath, "w", encoding="utf-8")
    for key, value in clean_sequences.items():
        print(f">{key}", file = fh)
        print(value, file = fh)

    fh.close()
    print("COL1A1 sequences have been cleaned")
    clean_num = len(clean_sequences)
    print(f"Number of clean COL1A1 sequences = {clean_num}")
    print("Output is COL1A1_seqs_clean_NCBI.fasta")
    print("######################################")

def run_clean_a1(file_name, output_dir):
    """Runs the full pipeline cleaning the cola1 sequences in the fasta file.

    Saves the fasta file in the output directory
    
    Parameters
    ----------
    file_name : str
        filepath for location of fasta file
    output_dir : str
        filepath for directory where the output file will be saved 
    """
    sequence_dict = read_fasta(file_name)
    clean_sequence_dict = clean_a1(sequence_dict)
    convert_fasta(clean_sequence_dict, output_dir)

if __name__ == "__main__":
    sys.exit(run_clean_a1(sys.argv[1:]))