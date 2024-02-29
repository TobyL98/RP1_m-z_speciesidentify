filename = 'RatSeqeunce.txt'

fileObj = open(filename, 'r')

sequence = ""
for line in fileObj:
    line = line.strip('\n')
    sequence += line

fileObj = open(filename, 'w')
fileObj.write(sequence)
