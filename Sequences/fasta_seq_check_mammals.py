#####################
## fasta_seq_check_mammals.py
#####################

# reads in the fasta file taken from e.g, ensembl or NCBI entrez
# example in github is COL1A2 Final_filtered
# converts to dictionary where key is identifier and value is sequence
# looks at start and end of sequence and length of sequence
# to make potential rules for filtering



import pandas as pd
import re

def readFastaFile(fileName):
    fileObj = open(fileName, 'r')
    sequences = {}   #  a dict, to contain all our sequences ...
    names = []       #  ... and a list to contain their names
    count = 0 
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            header = line[1:].rstrip('\n')  # remove newline, ignore '>
            class_name = re.split(r'[-_]', header)[0]
            if class_name == "Mammalia":
                name = header.split('[')[-1] # obtain the species name
                name = name.strip(']') #remove unwanted second ]

                if name in names:
                    count +=1 
                    name = name + " " + str(count)
                names.append(name)
                sequences[name] = ''
        elif class_name == "Mammalia":
            sequences[name] += line.rstrip('\n')  # concat to growing dict value
        else:
            continue
    fileObj.close()
    return (sequences)

file = 'COL1A2_seqs.fasta'
output_seq = readFastaFile(file)
#print(output_seq)


# start of the dictionaries ofr the rules
seq_three_dict = {}
seq_three_list = []

seq_five_dict = {}
seq_five_list = []

seq_last_three_dict = {}
seq_last_three_list = []

seq_last_five_dict = {}
seq_last_five_list = []

seq_length_dict = {}
seq_length_list = []

## looking at rules I can make
for seq in output_seq.values():
    seq_three = seq[:3]
    if seq_three in seq_three_list:
        seq_three_dict[seq_three] += 1
    else:
        seq_three_dict[seq_three] = 1
        seq_three_list.append(seq_three)

    seq_five = seq[:5]
    if seq_five in seq_five_list:
        seq_five_dict[seq_five] += 1
    else:
        seq_five_dict[seq_five] = 1
        seq_five_list.append(seq_five)

    seq_last_three = seq[-3:]
    if seq_last_three in seq_last_three_list:
        seq_last_three_dict[seq_last_three] += 1
    else:
        seq_last_three_dict[seq_last_three] = 1
        seq_last_three_list.append(seq_last_three)

    seq_last_five = seq[-5:]
    if seq_last_five in seq_last_five_list:
        seq_last_five_dict[seq_last_five] += 1
    else:
        seq_last_five_dict[seq_last_five] = 1
        seq_last_five_list.append(seq_last_five)

    seq_length = len(seq)
    if seq_length in seq_length_list:
        seq_length_dict[seq_length] += 1
    else:
        seq_length_dict[seq_length] = 1
        seq_length_list.append(seq_length)

print(seq_three_dict)
print(seq_five_dict)
print(seq_last_three_dict)
print(seq_last_five_dict)
print(sorted(seq_length_dict.items(), key= lambda item: item[1]))


