#####################
## fasta_seq_amendA1A2.py
#####################

# reads in the fasta file taken from e.g, ensembl or NCBI entrez
# filters by class 'Mammalia'
# creates dictionaries where key is species with taxonomic information
# and sequence is COL1A1 or COL1A2 sequence
# alters taxonomic information to a standard format
# filters COL1A1 data as not already cleaned

import pandas as pd
import re

def readFastaFile(fileName):
    fileObj = open(fileName, 'r')
    sequences = {}   #  a dict, to contain all our sequences ...
    names = []       #  ... and a list to contain their names 
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            header = line[1:].rstrip('\n')  # remove newline, ignore '>
            class_name = re.split(r'[-_]', header)[0]
            if class_name == "Mammalia":
                taxon = re.split(r'-[XN]P_', header)[0] # obtain the taxonomic information
                taxon_names = re.split(r'[-_]', taxon)
                #print(taxon_names)

                ## Organising all taxnomic names into standard format (quite confusing)
                # adding the class name (Mammalia)
                name = "CLASS=" + taxon_names[0]

                # When no theria subclass given EXAMPLES:
                # Mammalia-Prototheria-Monotremata
                if (taxon_names[1] == "Eutheria" or 
                taxon_names[1] == "Prototheria"):
                    name += "|SUBCLASS=Theria"
                    name += "|INFRACLASS=" + taxon_names[1]

                    # Ignoring extraq taxnomic information before ORDER. Examples:
                    # Mammalia_Eutheria_Euarchontoglires_Glires_Rodentia
                    # Mammalia_Eutheria_Laurasiatheria_Artiodactyla
                    if (taxon_names[2] == "Laurasiatheria" or
                        taxon_names[2] == "Euarchontoglires"):
                        if taxon_names[3] == "Glires":
                            name += "|ORDER" + taxon_names[4]
                        else:
                            name += "|ORDER" + taxon_names[3]
                    else:
                        name += "|ORDER" + taxon_names[2]
                # Most examples:
                # Mammalia-Theria-Eutheria-Carnivora
                # ORDER = Carnivora
                else:
                    name += "|SUBCLASS=" + taxon_names[1] 
                    name += "|INFRACLASS=" + taxon_names[2]
                    name += "|ORDER=" + taxon_names[3]

                for family_name in taxon_names:
                    if family_name.endswith("idae"):
                        name += "|FAMILY=" + family_name
                genus = taxon_names[-1].split(" ")[0]
                name += "|GENUS=" + genus
                name += "|SPECIES=" + taxon_names[-1]

                

                if name in names:
                    name = name
                names.append(name)
                sequences[name] = ''
        elif class_name == "Mammalia":
            sequences[name] += line.rstrip('\n')  # concat to growing dict value
        else:
            continue
    fileObj.close()
    return (sequences)

file = 'COL1A1_seqs.fasta'
output_COL1A1_dict = readFastaFile(file)
#print(output_seq)

output_COL1A1_final_dict = {}
for key, value in output_COL1A1_dict.items():
    if len(value) >= 1055 and len(value) <= 1058:
        output_COL1A1_final_dict[key] = value

#print(output_COL1A1_final_dict)
print(len(output_COL1A1_final_dict))

file2 = 'COL1A2_seqs.fasta'
output_COL1A2_dict = readFastaFile(file2)

print(output_COL1A1_final_dict)
print(len(output_COL1A1_final_dict))





