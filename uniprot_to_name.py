import argparse
from Bio import ExPASy
from Bio import SwissProt

def get_protein_names(uniprot_ids_file):
    protein_names = []

    with open(uniprot_ids_file, "r") as file:
        uniprot_ids = [line.strip() for line in file]

    for uniprot_id in uniprot_ids:
        try:
            handle = ExPASy.get_sprot_raw(uniprot_id)
            record = SwissProt.read(handle)
            protein_names.append(record.entry_name)
        except Exception as e:
            print(f"Error fetching data for UniProt ID {uniprot_id}: {str(e)}")

    return protein_names

def main():
    parser = argparse.ArgumentParser(description="Get protein names from UniProt IDs")
    parser.add_argument("uniprot_ids_file", help="Path to a text file containing line-separated UniProt IDs")
    parser.add_argument("output", help="Output file in text format")
    args = parser.parse_args()

    protein_names = get_protein_names(args.uniprot_ids_file)
    output_file = args.output

    with open(output_file, "w") as output:
        for protein_name in protein_names:
            output.write(f"{protein_name}\n")

    print(f"Protein names saved to {output_file}")

if __name__ == "__main__":
    main()
