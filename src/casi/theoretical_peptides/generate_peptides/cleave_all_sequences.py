"""cleave_all_sequences.py

This code will split a COL1 sequence input into
its respetive peptides allowing up to 1 cleavage.
Also will calculate the mass of the peptide and
the mass of all the combinations of hydroxylations
and deamidations

Input: csv with species name and COL sequence
Output: Dataframe with peptide sequence, mass,
num of hydroxylations, num of deamidations."""

from pathlib import Path
import sys
import time
from collections import namedtuple

import pandas as pd
from tqdm import tqdm

from casi.theoretical_peptides.generate_peptides import cleave_mass
from casi.theoretical_peptides.sort_sequences import merge_cola1a2


def run_cleave_mass(collagen_seq: str, species_info: namedtuple) -> pd.DataFrame:
    """In silico trypsin digest of the COL1 sequence
    and mass calculation. Includes common
    post-translational modifications."""

    # cleave and calculate mass
    collagen_pep_df = cleave_mass.cleave_and_mass(collagen_seq, "trypsin", 1)
    # add species and genus information
    collagen_pep_df["species"] = species_info.species
    collagen_pep_df["genus"] = species_info.genus
    collagen_pep_df["subfamily"] = species_info.subfamily
    collagen_pep_df["family"] = species_info.family
    collagen_pep_df["order"] = species_info.order

    return collagen_pep_df


def collagen_peptide_mass(col_dict, output_folder):
    """Runs throughs each input csv. Extracts the sequence.
    Sequence is then in silico cleaved and masses for
    each peptide is calculated

    Input: dataframe with collagen sequence and species information
    Output: csv file with generated collagen peptides for each peptide"""

    print("""Generating possible peptides from trypsin digest and
their possible masses. This code will take longer to run.""")

    output_folder = output_folder / "unfiltered_peptides"
    output_folder.mkdir(exist_ok=True)
    total_iterations = len(col_dict)
    for key, value in tqdm(
        col_dict.items(),
        total=total_iterations,
        desc="Generating Theoretical Peptides",
    ):
        time.sleep(0.01)
        # simulate trypins digest and calculate peptide masses
        col_pep_df = run_cleave_mass(value, key)
        # save to correct folder
        species_name = key.species.replace(" ", "_")
        output_name = f"{species_name}_col1_peptides.csv"
        output_filename = output_folder / output_name
        col_pep_df.to_csv(output_filename)

if __name__ == "__main__":
    sys.exit()
