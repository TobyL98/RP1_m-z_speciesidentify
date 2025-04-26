"""
cleave_all_sequences.py

This code will split a COL1 sequence input into
its respetive peptides allowing up to 1 cleavage.
Also will calculate the mass of the peptide and
the mass of all the combinations of hydroxylations
and deamidations
"""

from pathlib import Path
import sys
import time
from collections import namedtuple

import pandas as pd
from tqdm import tqdm

from casi.theoretical_peptides.generate_peptides import cleave_mass


def run_cleave_mass(collagen_seq: str,
                    species_info: namedtuple) -> pd.DataFrame:
    """
    In silico trypsin digest of the COL1 sequence
    and mass calculation. Includes common
    post-translational modifications.

    args
        collagen_seq (str): collagen sequence
        species_info (namedtuple): taxonomic inffromation of the species

    returns
        collagen_pep_df (pd.DataFrame): dataframe with peptide sequence, mass,
        num of hydroxylations, num of deamidations and taxonomic information.
    """

    # cleave and calculate mass
    collagen_pep_df = cleave_mass.cleave_and_mass(collagen_seq, "trypsin", 1)
    # add taxonomic information
    collagen_pep_df["species"] = species_info.species
    collagen_pep_df["genus"] = species_info.genus
    collagen_pep_df["subfamily"] = species_info.subfamily
    collagen_pep_df["family"] = species_info.family
    collagen_pep_df["order"] = species_info.order

    return collagen_pep_df


def collagen_peptide_mass(col_dict: dict, output_folder: Path):
    """
    Processes each collagen sequence in the provided dictionary,
    performing in silico digestion and mass calculation for each peptide.

    Args:
        col_dict (dict): Dictionary with species taxonomic information as keys
                         and collagen sequences as values.
        output_folder (Path): Path to the folder where the output files will be saved.

    Generates:
        CSV files for each species containing the theoretical peptides
        and their calculated masses.
    """

    print("Generating possible peptides from trypsin digest and their possible masses. This code will take longer to run.")

    output_folder = output_folder / "unfiltered_peptides"
    output_folder.mkdir(exist_ok=True)

    # Total number of iterations for progress tracking
    total_iterations = len(col_dict)

    # Iterate through each entry in the collagen dictionary
    for key, value in tqdm(col_dict.items(), total=total_iterations, desc="Generating Theoretical Peptides"):
        time.sleep(0.01)

        # Perform trypsin digest and calculate peptide masses
        col_pep_df = run_cleave_mass(value, key)

        # Construct the output filename using species name
        species_name = key.species.replace(" ", "_")
        output_name = f"{species_name}_col1_peptides.csv"
        output_filename = output_folder / output_name

        col_pep_df.to_csv(output_filename)

if __name__ == "__main__":
    sys.exit()
