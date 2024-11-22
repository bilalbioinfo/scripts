#!/bin/usr/env python3
import sys
import gzip

"""
Author: Bilal Sharif <bilal.bioinfo@gmail.com>
"""

if len(sys.argv) != 3:
    print(f"Usage: python3 {sys.argv[0]} <filename.haplo.gz> <filename.fasta>")
    sys.exit(1)

haplo_file = sys.argv[1]
fasta_file = sys.argv[2]

def haplo_to_fasta(haplo_file):
    """
    Read a haplotype file and return a dictionary with the sequences for each individual.
    """
    sequences = {}
    with gzip.open(haplo_file, "rt") as haplo_file_handle:
        ## read header and get the individuals
        header = haplo_file_handle.readline().strip().split("\t")
        inds = header[3:]
        for ind in inds: # Skip the first three columns
            sequences[ind] = []

        ## read data and get the sequences
        for line in haplo_file_handle:
            cols = line.strip().split("\t")
            genotypes = cols[3:]
            for ind, base in zip(inds, genotypes):
                sequences[ind].append(base)
    return sequences

sequences = haplo_to_fasta(haplo_file)

# Write fasta file
with open(fasta_file, "w") as fasta_file_handle:
    for ind, seq in sequences.items():
        fasta_file_handle.write(f">{ind}\n")  # Write the header
        for i in range(0, len(seq), 75):
            fasta_file_handle.write("".join(seq[i:i+75]) + "\n")
