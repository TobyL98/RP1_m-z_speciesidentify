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
import re
from pathlib import Path
from collections import namedtuple

import pandas as pd
import taxopy

RankLineage = namedtuple(
    "RankLineage",
    [
    "species",
    "genus",
    "subfamily",
    "family",
    "order",
    ]
)

def read_col_fasta(file_name):
    """Reads the fasta file"""
    file_obj = open(file_name, "r")
    sequences = {}
    for line in file_obj:
        if line.startswith(">"):
            name = line[1:].rstrip("\n")
            sequences[name] = ""
        # adds the sequences lines as the dictionary value
        else:
            sequences[name] += line.rstrip("\n")

    file_obj.close()
    return sequences

def get_taxa(sequences: dict) -> dict:
    """
    Gets the taxonomic information from the key in the sequences dictionary. The key need to
    have the species name in the format of 'GENUS=genus|SPECIES=species_name'.
    The taxonomic information is then taken from the NCBI taxon database and added as
    a RankLineage object which becomes the key in the new sequences dictionary

    Args:
        sequences (dict): Dictionary with species names as key and sequence as value

    Returns:
        new_sequences (dict): Dictionary with RankLineage object as key and sequence as value
    """
    
    # get taxopy database
    print("\nGetting taxon information. Takes a while to dowmload NCBI taxon database")
    new_sequences: dict = {}

    # downloads NCBI taxonomy database
    taxdb = taxopy.TaxDb()

    for key, value in sequences.items():

        # extracts species name and gets taxon information
        key_list = re.split("[=|]", key)
        species_name = key_list[1]
        genus_name = key_list[-1]

        try:
            taxon_id = taxopy.taxid_from_name(genus_name, taxdb)[0]
        except IndexError:
            print(f"Genus {genus_name} not found in NCBI taxon database")
            continue

        taxon_object = taxopy.Taxon(taxon_id, taxdb)
        lineage = taxon_object.ranked_name_lineage

        # create dictionary with taxonomic ranks information
        ranks: dict = {}
        for rank in lineage:
            ranks[rank[0]] = rank[1]

        # create RankLineage object
        lineage = RankLineage(
            species=species_name,
            genus=ranks.get("genus"),
            subfamily=ranks.get("subfamily"),
            family=ranks.get("family"),
            order=ranks.get("order"),
        )
        print(lineage)

        new_sequences[lineage] = value
    return new_sequences
    
def change_header(sequences):
    """Rearranges the fasta file header to obtain the species, genus and
    the family and outptuts into a readable header."""

    new_sequences = {}
    for key, value in sequences.items():
    # gets out the species name
        if key.endswith("]"):
            species_name = key.split("[")[1]
            species_name = species_name.strip("]")
            # assigns standard format for species name
            genus = species_name.split(" ")[0]
            name = "SPECIES=" + species_name
            name += "|GENUS=" + genus
            # create new dictionary with new header
            new_sequences[name] = value
    del sequences
    return new_sequences

# function combines the COL1A1 and COL1A2 sequences to one sequence
def COLA1A2combine(output_dir):
    # formats the A1 sequences
    file = output_dir / "COL1A1_seqs_clean_NCBI.fasta"
    col1a1_dict = read_col_fasta(file)

    # formats COL1A2 sequence
    file2 = output_dir / "COL1A2_seqs_clean_NCBI.fasta"
    col1a2_dict = read_col_fasta(file2)

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

def col1a2_combine():
    """Combines COl1a1 and COL1A2"""
    output_dir = Path("data/outputs")
    # formats the A1 sequences
    file = output_dir / "COL1A1_seqs_clean_NCBI.fasta"
    col1a1_dict = read_col_fasta(file)
    col1a1_dict = change_header(col1a1_dict)
    col1a1_dict = get_taxa(col1a1_dict)


if __name__ == "__main__":
    sys.exit(col1a2_combine())
