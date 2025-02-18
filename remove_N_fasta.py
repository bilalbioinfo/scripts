from Bio import SeqIO
from sys import argv
import os

'''
Author: Bilal Sharif <bilal.bioinfo@gmail.com>
Date: 2025-02-13
Description: This script removes 'N' bases from the aligned full-genome fasta files.
It works well with big fasta files. it process the fasta files in chunks of 100M bases. IMPORTANT: modify the chunk size as needed in the script.
The script reads the input FASTA file and checks if all sequences are of equal length.
The script writes the filtered sequences to temporary files and then merges the temporary files to create the final output FASTA file.
The script also takes an argument for the number of 'N' bases allowed in the sequences.
'''

if len(argv) != 4:
    raise ValueError("Usage: python remove_N_fasta.py <input_fasta> <output_fasta> <missing_allowed_number>")
    exit()

## Input arguments
fasta_file_in = argv[1]
fasta_file_out = argv[2]
missing_allowed = int(argv[3])

## Create a directory to store temporary files
os.makedirs("tmp_files", exist_ok=True)

## Initialize variables
start_pos = 0
chunk_size = 100_000_000 ## modify this value to change the chunk size as needed
chunk_num = 1

## Sanity check: check if all sequences are of equal length, get the length of the sequences and sequence IDs
print(f"Checking if all sequences are of equal length in {fasta_file_in}.....................:", flush=True)
with open(fasta_file_in, "r") as f1:
    seq_ids = []
    seq_length = None
    for record in SeqIO.parse(f1, "fasta"):
        seq_ids.append(record.id)
        if seq_length is None:
            seq_length = len(record.seq)
        elif len(record.seq) != seq_length:
            raise ValueError(f"Error: Sequences are not of equal length! Check {record.id}")
    if seq_length is None:
        raise ValueError("Error: No sequences found in the input FASTA file")
    if len(seq_ids) < missing_allowed:
        raise ValueError(f"Error: Missing allowed number {missing_allowed} is greater than the number of sequences {len(seq_ids)}")
print(f"All sequences are of equal length: {seq_length}. Proceeding further", flush=True)



def remove_bases(fasta_file_in, start_pos, end_pos, chunk_num, missing_allowed=0):
    rm_pos = {}
    # Step 1: Read and collect positions of 'N' in the chunk
    print (f"\nCollecting 'N' positions in samples: ", flush=True, end="")
    with open(fasta_file_in, "r") as f1:
        for record in SeqIO.parse(f1, "fasta"):
            seq = str(record.seq)
            seq_chunk = seq[start_pos:end_pos]
            for i, base in enumerate(seq_chunk):
                if base == 'N':
                    rm_pos[i] = rm_pos.get(i, 0) + 1
            print(f"{record.id} ", flush=True, end="")
    print(f"\nFound all 'N' positions in the sequences", flush=True)

    # Step 2: Write filtered sequences to the output file for the current chunk
    print(f"\nWriting filtered sequences: " , flush=True, end="")
    with open(fasta_file_in, "r") as f1:
        for record in SeqIO.parse(f1, "fasta"):
            filtered_seq = ""
            seq = str(record.seq)
            seq_chunk = seq[start_pos:end_pos]
            for i, base in enumerate(seq_chunk):
                if rm_pos.get(i, 0) <= missing_allowed:
                    filtered_seq += base
            fasta_file_out = f"tmp_files/{record.id}_chunk_{chunk_num}.fasta"
            with open(fasta_file_out, "w") as f2:
                f2.write(f">{record.id}\n")
                for i in range(0, len(filtered_seq), 50):
                    f2.write(filtered_seq[i:i+50] + "\n")
                print(f"{record.id} ", flush=True, end="")
    print(f"\nDone writing chunk {chunk_num} to {fasta_file_out}", flush=True)



while start_pos < seq_length:
    end_pos = start_pos + chunk_size
    if end_pos > seq_length:
        end_pos = seq_length  # Don't exceed the total number of bases
    print(f"\nPROCESSING chunk {chunk_num} from {start_pos}bp to {end_pos}bp.....................................:", flush=True)
    remove_bases(fasta_file_in, start_pos, end_pos, chunk_num, missing_allowed)

    start_pos = end_pos
    chunk_num += 1

## Merge all chunks per sample
for sample_id in seq_ids:
    with open(f"tmp_files/{sample_id}_all_chunks.fasta", "w") as f3:
        f3.write(f">{sample_id}\n")
        for chunk in range(1, chunk_num):
            with open(f"tmp_files/{sample_id}_chunk_{chunk}.fasta", "r") as f4:
                for record in SeqIO.parse(f4, "fasta"):
                    seq = str(record.seq)
                    for i in range(0, len(seq), 50):
                        f3.write(seq[i:i+50] + "\n")

## Merge all samples into the final output FASTA
with open(fasta_file_out, "w") as final_out:
    for sample_id in seq_ids:
        with open(f"tmp_files/{sample_id}_all_chunks.fasta", "r") as f5:
            for record in SeqIO.parse(f5, "fasta"):
                seq = str(record.seq)
                final_out.write(f">{record.id}\n")
                for i in range(0, len(seq), 50):
                    final_out.write(seq[i:i+50] + "\n")


print(f"\n\nSuccess Run: Final merged FASTA file saved as {fasta_file_out}", flush=True)
