#! /usr/bin/env python

import sys
import csv
import os

input_file = "genes.fna"

output_file_name = input_file.strip(".fna")+"_nucl_counts.csv"

NucDict = {
    "A": 0,  "C": 0, "G": 0, "T": 0}

colnames = ["ID"]
colnames.extend(NucDict.keys())
# Create output file, write headers
output_file = open(output_file_name, "w")
wr = csv.writer(output_file, quoting=csv.QUOTE_ALL)
wr.writerow(colnames)

output_line = list()
seq = str()
ID = str()
row = 0

with open(input_file, 'r') as list_paths:
    # For each line, identify if it is a label or string of AA and parse accordingly
    for line in list_paths:
        if line.startswith(">"):
            # If it's the very first line, go straight to parse
            if row == 0:
                next
            else:
                for i in range(0, len(seq)):
                    nuc = seq[i]
                    NucDict[nuc] = NucDict[nuc] + 1
                output_line.append(ID)
                output_line.extend(list(NucDict.values()))
                # Marks the end of a sequeence and the start of a new one
                # Append the last sequence to the DF
                wr.writerow(output_line)
                # Re-initialize dictionary so everything is back to zero
                seq = str()
                NucDict = dict.fromkeys(NucDict, 0)
                output_line = list()
            row += 1
            ID = line.strip("\n").replace(" # ", ";").split(';')[0]
            ID = ID.replace('>','')
        else:
            seq = seq+line.strip("\n")
    for i in range(0, len(seq)):
        nuc = seq[i]
        NucDict[nuc] = NucDict[nuc] + 1
    output_line.append(ID)
    output_line.extend(list(NucDict.values()))
    wr.writerow(output_line)

output_file.close()
