### This script is for making the "chains" file in all_vs_all mode 
### Only used when doing analysis of complexes (with multimer_colabfold_analysis.py)

import sys

# Parse command line options

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
# Iterates through each line in the file and store the identifier and sequence in a dictionary with fasta as the value
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
                count1 = seq1.count(":") + 1
                count2 = seq2.count(":") + 1
                output_file.write(f"{id1}_{id2},{count1},{count2}\n")
