"""mass_calculator.py

This code will split a COL1 sequence input into
it respetive peptides allowing up to 1 cleavage.
Also will calculate the mass of the peptide and
the mass of all the combinations of hydroxylations
and deamidations

Input: csv with species name and COL sequence
Output: Dataframe with peptide sequence, mass,
num of hydroxylations, num of deamidations."""

from pathlib import Path

import pandas as pd
import numpy as np
from pyteomics import parser, mass

def peptide_creator(seq):
    """Creates the peptides that would
    be formed from trypsin digest using
    expasy rules.
    
    Input: Full sequence (str)
    Output: set of peptides"""

    peptides_set = parser.cleave(seq, parser.expasy_rules['trypsin'], 1)

    return peptides_set


def peptides_mass_calculator(peptides):
    """calculates the masses for
    a set of peptides"""

    for peptide in peptides:
        print(peptide)

def main():
    """Runs the main part of
    the program"""

    # find correct input path and read csv
    folder_path = Path(__file__).parents[1]
    input_path = folder_path / "Sequences/sequences_taxon_NCBI.csv"
    col_names = ["sequence", "GENUS", "SPECIES"]
    input_df = pd.read_csv(input_path, sep=",", usecols=col_names)
    
    for index, row in input_df.iterrows():

        # generate peptides from trypsin digest
        peptides = peptide_creator(row["sequence"])
        for peptide in peptides:
            # calculate mass
        
        
    

if __name__ == "__main__":
    main()