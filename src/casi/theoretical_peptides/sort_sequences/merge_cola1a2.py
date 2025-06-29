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
from collections import namedtuple

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

def read_col_fasta(file_name: Path) -> dict:
    """
    Reads the fasta file

    Args:
        file_name (str): The location of the fasta file

    Returns:
        sequences (dict): Dictionary with header as key and sequence as value
    """
    file_obj = open(file_name, "r", encoding="utf-8")
    sequences = {}
    for line in file_obj:
        if line.startswith(">"):
            # gets out the species name from the header
            name = line[1:].rstrip("\n")
            sequences[name] = ""
        # adds the sequences lines as the dictionary value
        else:
            sequences[name] += line.rstrip("\n")

    file_obj.close()
    return sequences

def get_species(sequences: dict) -> dict:
    """
    Reads in the sequences dictionary and creates a new sequences dictionary
    with the species name and genus in a list as the key. Ectracts the species name 
    from the current key.

    Args:
        sequences (dict): Dictionary with header as key and sequence as value

    Returns:
        new_sequences (dict): Dictionary with species name, genus name list 
        as a key and sequence as value
    """
    new_sequences = {}
    for key, value in sequences.items():
        # gets out the species name
        if key.endswith("]"):
            species_name = key.split("[")[1]
            species_name = species_name.strip("]")
            # assigns standard format for species name
            #name = "SPECIES=" + species_name
            #name += "|GENUS=" + genus
            # create new dictionary with new header
            new_sequences[species_name] = value
    del sequences
    return new_sequences

def merge_col(col1a1_sequences: dict, col1a2_sequences: dict) -> dict:
    """
    Merges the COl1A1 and COl1A2 sequences on the dictionary key (species name).

    Args:
        col1a1_sequences (dict): Dictionary with species names as key and COl1A1
                                 sequence as value
        col1a2_sequences (dict): Dictionary with species names as key and COl1A2
                                 sequence as value

    Returns:
        a1a2_combined (dict): Dictionary with species names as key and merged
                              COl1A1 and COl1A2 sequences as value
    """
    a1a2_combined = {}
    for key in col1a1_sequences:
        if key in col1a2_sequences:
            a1a2_combined[key] = (
                col1a1_sequences[key] + "R" + col1a2_sequences[key]
            )

    return a1a2_combined

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

        genus = key.split(" ")[0].strip()

        # gets taxonomy taxon id from species or genus name
        try:
            taxon_id = taxopy.taxid_from_name(key, taxdb)[0]
        except IndexError:
            try:
                taxon_id = taxopy.taxid_from_name(genus, taxdb)[0]
            except IndexError:
                print("No taxon ID found for", key)
                taxon_id = None

        taxon_object = taxopy.Taxon(taxon_id, taxdb)
        lineage = taxon_object.ranked_name_lineage

        # create dictionary with taxonomic ranks information
        ranks: dict = {}
        for rank in lineage:
            ranks[rank[0]] = rank[1]

        # create RankLineage object
        lineage = RankLineage(
            species=key,
            genus=ranks.get("genus"),
            subfamily=ranks.get("subfamily"),
            family=ranks.get("family"),
            order=ranks.get("order"),
        )

        new_sequences[lineage] = value
    return new_sequences

def create_fasta(sequences: dict, output_file: Path) -> None:
    """
    Creates a fasta file with the sequences as values and the taxonomic information
    as the header. The header is in the format of:
    >SPECIES=species_name|GENUS=genus_name|SUBFAMILY=subfamily_name|
    FAMILY=family_name|ORDER=order_name

    Args:
        sequences (dict): Dictionary containing sequences with taxonomic information
        output_file (Path): Path to the output fasta file
    """
    
    with open(output_file, "w", encoding="utf-8") as f:
        for key, value in sequences.items():
            # Write the header with taxonomic information
            print(f">SPECIES={key.species}|"
                  f"GENUS={key.genus}|"
                  f"SUBFAMILY={key.subfamily}|"
                  f"FAMILY={key.family}|"
                  f"ORDER={key.order}",
                  file=f)
            # Write the sequence value
            print(value, file=f)

    print(
        "COL1A1 and COL1A2 sequences have been combined and taxonomic information added."
    )
    print(f"Number of complete COL1A1 and COL1A2 Sequences = {len(sequences)}")
    print("Output is COL1A1A2_combined_seqs_NCBI.fasta")
    print("######################################")

def col1a1a2_combine(output_dir: Path) -> dict:
    """
    The function reads the cleaned COl1A1 and COL1A2 sequences from the
    output directory, formats them, and merges them by species. It then
    retrieves the taxonomic information for the combined sequences and
    creates a fasta file with the taxonomic headers
    
    Args:
        output_dir (path): the path to the output directory
    """
    # Read and format the COl1A1 sequences
    col1a1_file = output_dir / "COL1A1_seqs_clean_NCBI.fasta"
    col1a1_dict = read_col_fasta(col1a1_file)
    # Get the species information for the COl1A1 sequences
    col1a1_dict = get_species(col1a1_dict)

    # Read and format the COL1A2 sequences
    col1a2_file = output_dir / "COL1A2_seqs_clean_NCBI.fasta"
    col1a2_dict = read_col_fasta(col1a2_file)
    # Get the species information for the COL1A2 sequences
    col1a2_dict = get_species(col1a2_dict)

    # Merge COl1A1 and COL1A2 sequences by species
    col1a2_combined = merge_col(col1a1_dict, col1a2_dict)

    # Retrieve taxonomic information for the combined sequences
    col1a2_combined = get_taxa(col1a2_combined)
    
    # Create and save the combined fasta file with taxonomic headers
    output_file = output_dir / "COL1A1A2_combined_seqs_NCBI.fasta"
    create_fasta(col1a2_combined, output_file)

    return col1a2_combined

     


if __name__ == "__main__":
    sys.exit()
