### This script is used to convert a list of FASTA sequences to a pairwise ensemble of interactions

import sys

# Parse input variables

if len(sys.argv) != 3:
    print("Usage: {} <input_fasta_file> <output_fasta_file>".format(sys.argv[0]))
    sys.exit(1)

input_fasta_file = sys.argv[1]
output_fasta_file = sys.argv[2]

try:
    with open(input_fasta_file, 'r') as file:
        lines = file.readlines()
except FileNotFoundError:
    print("Error: File not found: {}".format(input_fasta_file))
    sys.exit(1)

ids_and_seqs = {}

# Read the input FASTA file
# Iterates through each line in the file and stores the identifier and sequence in a dictionary with fasta as the value
for line in lines:
    if line.startswith(">"):
        if "|" in line:
            prev_identifier = line.split('|')[1].strip()
        else:
            prev_identifier = line.replace(">", "").strip()
    elif not line.startswith(">"):
        ids_and_seqs[prev_identifier] = line.strip()

# Create the pairwise combinations and write to the output file
with open(output_fasta_file, 'w') as output_file:
    for id1, seq1 in ids_and_seqs.items():
        for id2, seq2 in ids_and_seqs.items():
            if len(seq1) + len(seq2) > 4000:
                print(f"Skipping {id1}_{id2} due to length")
            else:
                output_file.write(f">{id1}_{id2}\n{seq1}:{seq2}\n")
