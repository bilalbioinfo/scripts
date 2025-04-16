#!/usr/bin/env python3
import sys
"""
Author = "Bilal Sharif: bilal.bioinfo@gmail.com"
This script is used to create a feature table for NCBI annotation upload from a tsv file. The tsv file should have the following columns:
1. row number/ID - by default Genious exports Type as row ID - this column is not used in this script!!
2. Sample name
3. Start position
4. End position
5. NCBI feature
6. Gene
7. Product
8. Translation table
9. Translation except
10. Codon start
11. Direction
"""
## check if tsv file is provided
if len(sys.argv) != 2:
    print("Usage: python3 create_ncbi_feature_table.py <tsv_file> > feature_table.txt")
    sys.exit(1)
## function to parse tsv file and print feature table
def parse_tsv(tsv_file):
    last_sample = ''
    with open(tsv_file, 'r') as f1:
        f1.readline()  ## skip header
        for line in f1:
            fields = line.strip().split('\t')
            ## extract fields
            sample = fields[1]
            min_pos = int(fields[2])
            max_pos = int(fields[3])
            ncbi_feature = fields[4]
            gene = fields[5]
            product = fields[6]
            transl_table = fields[7]
            transl_except = fields[8]
            codon_start = fields[9]
            direction = fields[10]
            ## skip if ncbi_feature is empty as these are duplicates
            if not ncbi_feature:
                continue
            ## print sample header
            if sample != last_sample:
                print(f">Feature {sample}")
            last_sample = sample
            # print feature details based on direction
            if direction == 'forward':
                print(f"{min_pos}\t{max_pos}\t{ncbi_feature}")
            else:
                print(f"{max_pos}\t{min_pos}\t{ncbi_feature}")
            # print additional feature details if available
            if gene:
                print(f"\t\t\tgene\t{gene}")
            if product:
                print(f"\t\t\tproduct\t{product}")
            if transl_table:
                print(f"\t\t\ttransl_table\t{transl_table}")
            if transl_except:
                print(f"\t\t\ttransl_except\t{transl_except}")
            if codon_start:
                print(f"\t\t\tcodon_start\t{codon_start}")
## run script
tsv_file = sys.argv[1]
parse_tsv(tsv_file)
