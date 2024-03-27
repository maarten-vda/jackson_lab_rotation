#!/bin/bash
### The file goa_human.gaf is found at https://current.geneontology.org/annotations/goa_human.gaf.gz ###

# Check if the required arguments are provided
if [ "$#" -ne 2 ]; then
  echo "Usage: $0 <GO_term> <output_filepath>"
  exit 1
fi

# Extract input variables
go_term="$1"
output_filepath="$2"

# Check if the GO term starts with "GO:"
if [[ $go_term == GO:* ]]; then
  # Remove "GO:" from the beginning if present
  go_term="${go_term#GO:}"
fi

# Check if the GO term is a valid string of numbers
if ! [[ $go_term =~ ^[0-9]+$ ]]; then
  echo "Invalid GO term. Please provide a valid GO term as a string of numbers or GO:##### format."
  exit 1
fi

# Execute the grep and awk command with the provided GO term, excluding lines with "NOT" in the 4th column
grep "GO:$go_term" goa_human.gaf | awk '$4 !~ /NOT/ {print $2}' | sort -u > "$output_filepath"

echo "Output saved to: $output_filepath"

