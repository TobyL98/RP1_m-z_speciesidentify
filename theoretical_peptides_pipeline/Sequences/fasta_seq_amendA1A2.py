#####################
## fasta_seq_amendA1A2.py
#####################

# reads in the fasta file taken from e.g, ensembl or NCBI entrez
# does this for both COL1A1 and COL1A2
# filters by class 'Mammalia'
# creates dictionaries where key is species with taxonomic information
# and sequence is COL1A1 or COL1A2 sequence
# alters taxonomic information to a standard format
# combines the COL1A1 and COL1A2 dictionary values (sequences) by keys
# outputs as a fasta file

import pandas as pd
import re

# function reads the fasta file and outputs as dictionary
# also standardises taxonomic information
def readCOLFastaFile(fileName):
    fileObj = open(fileName, 'r')
    sequences = {}   #  a dict, to contain all our sequences ...
    names = []       #  ... and a list to contain their names 
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            header = line[1:].rstrip('\n')  # remove newline, ignore '>
            class_name = re.split(r'[-_]', header)[0]
            if class_name == "Mammalia":
                taxon = re.split(r'-[XN]P_', header)[0] # obtain the taxonomic information
                #splits to seperate different taxonomic information
                taxon_names = re.split(r'[-_]', taxon)
                

                ## Organising all taxnomic names into standard format (quite confusing)
                # adding the class name (Mammalia)
                name = "CLASS=" + taxon_names[0]

                # When no theria subclass given EXAMPLES:
                # Mammalia-Prototheria-Monotremata
                if (taxon_names[1] == "Eutheria" or 
                taxon_names[1] == "Prototheria"):
                    # adds in theria class when missing
                    name += "|SUBCLASS=Theria"
                    name += "|INFRACLASS=" + taxon_names[1]

                    # Ignoring extra taxanomic information before ORDER. Examples:
                    # Mammalia_Eutheria_Euarchontoglires_Glires_Rodentia
                    # Mammalia_Eutheria_Laurasiatheria_Artiodactyla
                    if (taxon_names[2] == "Laurasiatheria" or
                        taxon_names[2] == "Euarchontoglires"):
                        if taxon_names[3] == "Glires":
                            name += "|ORDER=" + taxon_names[4]
                        else:
                            name += "|ORDER=" + taxon_names[3]
                    else:
                        name += "|ORDER=" + taxon_names[2]
                # Most examples:
                # Mammalia-Theria-Eutheria-Carnivora
                # ORDER = Carnivora
                else:
                    name += "|SUBCLASS=" + taxon_names[1] 
                    name += "|INFRACLASS=" + taxon_names[2]
                    name += "|ORDER=" + taxon_names[3]

                # finding family name and then adding
                # all family names end 'idae'
                for family_name in taxon_names:
                    if family_name.endswith("idae"):
                        name += "|FAMILY=" + family_name
                #finding genus and species name
                genus = taxon_names[-1].split(" ")[0]
                name += "|GENUS=" + genus
                name += "|SPECIES=" + taxon_names[-1]

                # will only create new dictionary item if
                # name is not already in names list
                if name in names:
                    name = name
                names.append(name)
                sequences[name] = ''
        # adds the sequences lines as the dictionary value
        elif class_name == "Mammalia":
            sequences[name] += line.rstrip('\n')  # concat to growing dict value
        else:
            continue
    fileObj.close()
    return (sequences)

# function combines the COL1A1 and COL1A2 sequences to one sequence
def COLA1A2combine():
    # formats the A1 sequences
    file = 'Sequences/COL1A1_seqs_clean.fasta'
    output_COL1A1_dict = readCOLFastaFile(file)

    #print(len(output_COL1A1_dict))

    #formats COL1A2 sequence
    file2 = 'Sequences/COL1A2_seqs_clean.fasta'
    output_COL1A2_dict = readCOLFastaFile(file2)

    #print(len(output_COL1A2_dict))

    ### Combining COL1A1 and COL1A2 dictionaries
    # now one sequence with R (Arginine) in between
    # combine if they have same key (i.e., same species)
    A1A2_combined_dict = {}
    for key in output_COL1A1_dict:
        if key in output_COL1A2_dict:
            A1A2_combined_dict[key] = output_COL1A1_dict[key] + "R" + output_COL1A2_dict[key]

    # Converting back to fasta format
    fh = open("Sequences/COL1A1A2_combined_seqs.fasta", "w")
    for key in A1A2_combined_dict:
        # first line is key
        print(">{0}".format(key), file = fh)
        # add lines for sequence
        print(A1A2_combined_dict[key], file = fh)

    fh.close()
    print("COL1A1 and COL1A2 sequences have been combined and taxonomic information added.")
    print("Number of complete COL1A1 and COL1A2 Sequences = {0}".format(len(A1A2_combined_dict)))
    print("Output is COL1A1A2_combined_seqs.fasta")
    print("######################################")









