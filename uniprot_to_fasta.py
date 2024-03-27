### SCRIPT TO CONVERT A LIST OF UNIPROT IDS TO THE CORRESPONDING FASTA SEQUENCES VIA REQUESTING URL ON UNIPROT

import argparse
import os
import requests

# Get protein sequences from the accession number just using URL requests
def get_protein_sequence(accession):
    url = f"https://www.uniprot.org/uniprot/{accession}.fasta"
    response = requests.get(url)

    if response.status_code == 200:
        return response.text.split('\n', 1)[1].replace('\n', '')  # Extract sequence excluding the header
    else:
        print(f"Failed to retrieve sequence for UniProt accession {accession}")
        return None

def main():
    # Argument parser
    parser = argparse.ArgumentParser(description="Convert UniProt accession numbers to a single FASTA file with concatenated protein sequences.")
    parser.add_argument("-i", "--input", help="Path to a text file containing UniProt accession numbers, one per line", required=True)
    parser.add_argument("-o", "--output", help="Output file in FASTA format", required=True)
    args = parser.parse_args()

    # Read UniProt accession numbers from input file
    with open(args.input, "r") as input_file:
        accession_list = [line.strip() for line in input_file]

    # Concatenate protein sequences (allows processing as list)
    concatenated_sequence = ""
    for accession in accession_list:
        protein_sequence = get_protein_sequence(accession)
        if protein_sequence is not None:
            concatenated_sequence += f">UniProt|{accession}\n{protein_sequence}\n"

    # Write to the output file
    with open(args.output, "w") as output_file:
        output_file.write(concatenated_sequence)

    print(f"All protein sequences saved to {args.output}")

if __name__ == "__main__":
    main()
