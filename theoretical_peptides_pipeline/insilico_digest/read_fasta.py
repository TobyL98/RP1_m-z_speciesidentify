#####################
## read_fasta.py
#####################

# reads in the fasta file taken from e.g, ensembl or NCBI entrez
# example in github is COL1A2 Final_filtered
# and currently outputs just the sequence

## Toby note
# will need to amend to add in to turn into csv that cane be processed
# by run_parseseq

import pandas as pd

def readFastaFile(fileName):
    fileObj = open(fileName, 'r')
    sequences = {}   #  a dict, to contain all our sequences ...
    names = []       #  ... and a list to contain their names
    count = 0 
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            header = line[1:].rstrip('\n')  # remove newline, ignore '>'
            name_list = header.split()[-2:]       # take first and second of line, name of organism
            name = ' '.join(name_list)
            name = name.strip("[]")

            if name in names:
                count +=1 
                name = name + " " + str(count)
            names.append(name)
            sequences[name] = ''
        else:
            sequences[name] += line.rstrip('\n')  # concat to growing dict value
    fileObj.close()
    return (sequences)

file = 'sequence.fasta'
output_seq = readFastaFile(file)
print(output_seq)
