#####################
## fasta_to_csv_NCBI.py
#####################

# reads in the fasta file (in this case combined COL1A1 and COL1A2 file for NCBI data)
# outputs a csv with columns for the taxonomic information and the sequence



import pandas as pd

def Fasta_Database():
    fileName = 'COL1A1A2_combined_seqs_NCBI.fasta'
    
    fileObj = open(fileName, 'r')

    # create 
    new_file = open("database.fasta", "w")
    count = 0
    for line in fileObj:
        if line.startswith('>'):   # Ah ha ! a new sequence?
            count += 1

            header = line[1:].rstrip('\n')  # remove newline, ignore '>'

            ## gets taxonomic information format: 
            # CLASS=Mammalia|SUBCLASS=Theria|INFRACLASS=Prototheria|ORDER=Monotremata
            # assign each taxonomic level to list in dictionary
            taxons = header.split("|")
            info = taxons[1]
            name = info.split("=")[1]
            genus, species = name.split(" ")[0:2]
            letter1 = genus[0]
            letter2 = species[0].upper()
            initials = letter1 + letter2 + "COL1"
            number = "00000" + str(count)
            number = number[-5:]
            ID_num = "M" + number
            
            print(">cu|{0}|{1} {2} Collagen 1".format(ID_num, initials, genus), file= new_file)


            


        else:
            # adds lines to create sequence
            print(line.rstrip('\n'), file= new_file)
    fileObj.close()

Fasta_Database()

