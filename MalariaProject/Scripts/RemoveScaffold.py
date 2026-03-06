#!/usr/bin/env python3

"""
Author: Yaqi Jiao
Date: 27th February, 2026

RemoveScaffold.py
-----------------

Description: this script takes genome data as input, and removes sequence less than 3000 bp, and GC content larger than 35%

Input: genome 
Output: filtered genome data

Usage Example:
RemoveScaffold.py genomename outputname
"""

# package preparation
import sys


def fasta_reader(filename):
    """generator: read one by one, and return (header, sequence)"""
    header, seqs = None, []
    with open(filename, 'r') as f:
        for line in f:
            line = line.strip()
            if not line: continue
            if line.startswith(">"):
                if header:
                    yield header, "".join(seqs)
                header, seqs = line, []
            else:
                seqs.append(line)
        if header:
            yield header, "".join(seqs)

def calc_gc(seq):
    """calculate GC content"""
    if not seq: return 0
    # iterate sequence for once
    gc_count = sum(1 for base in seq if base in "GCgc")
    return (gc_count / len(seq)) * 100

def main():
    if len(sys.argv) < 3:
        print("Usage: RemoveScaffold.py <input_fasta> <output_fasta>")
        sys.exit(1)

    input_file, output_file = sys.argv[1], sys.argv[2]

    with open(output_file, "w") as out:
        for header, seq in fasta_reader(input_file):
            # calculate length directly
            length = len(seq)
            
            if length >= 3000:
                gc_content = calc_gc(seq)
                if gc_content <= 35:
                    gc_content = round(gc_content,2)
                    # write Header
                    out.write(f"{header}\tgc_content={gc_content}\n")
                    # 60 bases each line
                    for i in range(0, length, 60):
                        out.write(f"{seq[i:i+60]}\n")

if __name__ == "__main__":
    main()
