### This script is for making the "chains" file in pulldown mode, only use for analysing complexes using multimer_colabfold_analysis.py

#!/bin/bash


# Function for processing names for interactions 
process_first_line_seq() {
    local input_variable="$1"

    # Check if the variable starts with ">"
    if [[ $input_variable == ">"* ]]; then
        input_variable="${input_variable:1}"
    fi

    # Check if the variable contains "|", this is assuming its a UniProt ID
    if [[ $input_variable == *"|"* ]]; then
        # Split the string at "|"
        IFS='|' read -ra fields <<< "$input_variable"

        # Get the last field and assign it to the variable
        input_variable="${fields[@]: -1}"
    fi

    echo "$input_variable"
}

if [ "$#" -ne 3 ]; then
    echo "Usage: $0 <single_sequence.fasta> <multi_line.fasta> <output.fasta>"
    exit 1
fi

# Read the first line of the single sequence file
seq_file="$1"
first_line_seq=$(head -n 1 "$seq_file")
first_line_seq=$(process_first_line_seq "$first_line_seq")
second_line_seq=$(head -n 2 "$seq_file" | tail -n 1)
second_line_count=$(( $(echo "$second_line_seq" | grep -o ":" | wc -l) + 1 ))

echo "$( echo "${first_line_seq}_${first_line_seq},${second_line_count},${second_line_count}" | tr -d '[:space:]') " > $3

counter=0

# Iterate through the lines of the multi-line fasta file
while IFS= read -r line; do
    ((counter++))
    if ((counter % 2 != 0)); then
        line_id=$(process_first_line_seq "$line")
    elif ((counter % 2 == 0)); then
        if (( ${#second_line_seq} + ${#line} < 4000 )); then
            line_count=$(( $(echo "$line" | grep -o ":" | wc -l) + 1 ))
            echo "$( echo "${first_line_seq}_${line_id},${second_line_count},${line_count}" | tr -d '[:space:]') " >> $3
        else
            echo "${first_line_seq} and ${line_id} are over 4000 amino acids, skipping"
        fi
    fi
done < "$2"

