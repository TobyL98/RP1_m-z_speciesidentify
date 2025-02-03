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

import pandas as pd
from tqdm import tqdm

from casi.theoretical_peptides.generate_peptides import cleave_mass


def run_cleave_mass(col1_row):
    """In silico trypsin digest of the COL1 sequence
    and mass calculation. Includes common
    post-translational modifications."""

    # generate peptides from trypsin digest
    collagen_seq = col1_row["sequence"]
    # cleave and calculate mass
    collagen_pep_df = cleave_mass.cleave_and_mass(collagen_seq, "trypsin", 1)
    # add species and genus information
    collagen_pep_df["GENUS"] = col1_row["GENUS"]
    collagen_pep_df["SPECIES"] = col1_row["SPECIES"]
    return collagen_pep_df


def collagen_peptide_mass(output_folder):
    """Runs throughs each input csv. Extracts the sequence.
    Sequence is then in silico cleaved and masses for
    each peptide is calculated

    Input: dataframe with collagen sequence and species information
    Output: csv file with generated collagen peptides for each peptide"""

    print("""Generating possible peptides from trypsin digest and
their possible masses. This code will take longer to run.""")
    # find correct input path and read csv
    input_path = output_folder / "sequences_taxon_NCBI.csv"
    col_names = ["sequence", "GENUS", "SPECIES"]
    input_df = pd.read_csv(input_path, sep=",", usecols=col_names)

    output_folder = output_folder / "unfiltered_peptides"
    output_folder.mkdir(exist_ok=True)
    total_iterations = len(input_df)
    for index, row in tqdm(
        input_df.iterrows(),
        total=total_iterations,
        desc="Generating Theoretical Peptides",
    ):
        time.sleep(0.01)
        # simulate trypins digest and calculate peptide masses
        col_pep_df = run_cleave_mass(row)
        # save to correct folder
        species_name = row["SPECIES"].replace(" ", "_")
        output_name = f"{species_name}_col1_peptides.csv"
        output_filename = output_folder / output_name
        col_pep_df.to_csv(output_filename)


if __name__ == "__main__":
    sys.exit()
