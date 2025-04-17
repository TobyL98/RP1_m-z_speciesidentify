#####################
## fasta_to_csv_NCBI.py
#####################

# reads in the fasta file (in this case combined COL1A1 and COL1A2 file for NCBI data)
# outputs a csv with columns for the taxonomic information and the sequence

import sys
from pathlib import Path

import pandas as pd

def FastaToCSV(output_folder):
    file_name = output_folder / 'COL1A1A2_combined_seqs_NCBI.fasta'
    
    file_obj = open(file_name, 'r', encoding="utf-8")
    known_taxons = []       #  ... and a list to contain their names
    sequences = {"sequence": []} # creates sequence dictionary
    sequence_count = 0

    for line in file_obj:
        if line.startswith('>'):   # Ah ha ! a new sequence?

            if sequence_count > 0: #ensures first sequence created before append
                sequences["sequence"].append(sequence) #appends sequences to dictionary
            sequence = "" # resets the sequence when we reach new entry
            sequence_count = 1

            header = line[1:].rstrip('\n')  # remove newline, ignore '>'

            ## gets taxonomic information format: 
            # CLASS=Mammalia|SUBCLASS=Theria|INFRACLASS=Prototheria|ORDER=Monotremata
            # assign each taxonomic level to list in dictionary
            taxons = header.split("|")
            for info in taxons:
                taxon_list = info.split("=")
                
                if taxon in known_taxons: # if taxon already a key appends new name to list
                    sequences[taxon].append(name)
                else: # otherwise adds the name to list of known taxon key
                    sequences[taxon] = [name]
                    known_taxons.append(taxon)

        else:
            # adds lines to create sequence
            sequence += line.rstrip('\n')

    # ensure last sequence is appended
    sequences["sequence"].append(sequence)
    
    # convert dictionary to pandas dataframe
    df = pd.DataFrame(data = sequences)
    
    file_obj.close()
    output_file = output_folder / "sequences_taxon_NCBI.csv"
    df.to_csv(output_file)
    print("Fasta file converted to CSV")
    print("Output is sequences_taxon_NCBI.csv")
    print("######################################")

if __name__ == "__main__":
    sys.exit()



