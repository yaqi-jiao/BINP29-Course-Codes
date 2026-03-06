#!/usr/bin/python3                                                             
                                                                      
# A quick and dirty program              
# sys.argv[1] is input fasta file name
# sys.argv[2] is threshold GC content as integer
# sys.argv[3] is the output file
# sys.argv[4] is the minimum length for scaffolds to keep

# Run like: removeScaffolds.py Haemoproteus_tartakovskyi.genome 35 Ht.genome 3000
import sys

seq = {} # A hash with ids as keys and their sequences as values
ids = [] # A list of ids to keep the original order

with open(sys.argv[1], 'r') as fin:
    for line in fin:
        line = line.rstrip()
        if line[0] == '>':
            id = line.split(' ')[0] # Keep id up to the first space
            id = id[1:] # Remove the >
            seq[id] = '' 
            ids.append(id)
        else:
            seq[id] += line.upper()

with open(sys.argv[3], 'w') as fout:
    for id in ids:
        if len(seq[id]) < int(sys.argv[4]):
            continue
        gcCount = seq[id].count('G') + seq[id].count('C')
        atCount = seq[id].count('A') + seq[id].count('T')
        gc = gcCount / (gcCount + atCount)
        if gc <= float(sys.argv[2]) / 100:
            outGC = round(gc, 2)
            print('>' + id + ' GC=' + str(outGC) + ' Length=' + str(len(seq[id])), file=fout)
            print(seq[id], file=fout)
