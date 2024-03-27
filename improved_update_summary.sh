#!/bin/bash

export PATH="/exports/applications/gridengine/ge-8.6.5/bin/lx-amd64:/exports/applications/apps/SL7/modules/5.2.0/bin:/usr/local/bin:/usr/bin:/usr/local/sbin:/usr/sbin:/home/s2037423/.local/bin:/home/s2037423/bin"
source /etc/profile.d/modules.sh
module load igmm/apps/openssl/3.0.5
module load anaconda/5.0.1
conda activate /gpfs/igmmfs01/eddie/ajackson-wrkgrp/maarten/software/prodigy/prodigy_env


# Check if the correct number of arguments are provided
if [ "$#" -ne 1 ]; then
    echo "Usage: $0 output_folder"
    exit 1
fi

output_dir="$PWD/$1"
exec_dict="$PWD"

dr_sasa="/exports/igmm/eddie/ajackson-wrkgrp/maarten/software/drsasa2/build/dr_sasa"

mkdir $output_dir/drsasa

# Copy the header from summary.csv to the temporary file
head -n 1 "$output_dir/colabfold_predictions_analysis/summary.csv" | sed 's/$/,pTM,ipTM,ranking_confidence,interface_area,total_surface_area,affinity,kd25/' > "$output_dir/updated_summary.csv"

# Iterate through the rows of summary.csv (excluding the header)
tail -n +2 "$output_dir/colabfold_predictions_analysis/summary.csv" | while IFS=, read -r complex_name avg_n_models max_n_models num_contacts_with_max_n_models num_unique_contacts best_model_num best_pdockq best_plddt_avg best_pae_avg; do
    # Find the corresponding row in combined_metrics.csv
    result=$(awk -F, -v complex_name="$complex_name" -v best_model_num="$best_model_num" 'BEGIN {OFS=",";} \
    $1 == complex_name && $2 == best_model_num {print $(NF-2), $(NF-1), $NF}' "$output_dir/combined_metrics.csv")
    # Extract pTM, ipTM, and ranking_confidence from the result
    pTM=$(echo "$result" | cut -d',' -f1)
    ipTM=$(echo "$result" | cut -d',' -f2)
    ranking_confidence=$(echo "$result" | cut -d',' -f3)
    mkdir -p ${output_dir}/drsasa_output
    cd ${output_dir}/drsasa_output
    # Use find to dynamically locate the file matching the pattern
    pdb_file=$(find "${output_dir}/colabfold_predictions/" -name "${complex_name}_*model_${best_model_num}*.pdb" -print -quit)
    # Check if exactly one file is found
    if [ -f "$pdb_file" ]; then
        # Run dr_sasa only if the file is found
        $dr_sasa -m 1 -i "$pdb_file" -chain A -chain B > ${output_dir}/drsasa/areas_${complex_name}.txt
    else
        echo "Error: No or multiple matching PDB files found for ${complex_name}_*rank_00${best_model_num}*"
    fi
    echo "Saved interface areas for ${pdb_file} to ${output_dir}/drsasa/areas_${complex_name}.txt"
    interface_area=$(tail -n 1 ${output_dir}/drsasa/areas_${complex_name}.txt | awk '{ print $NF }')
    surface_area=$(head -n 3 ${output_dir}/drsasa/areas_${complex_name}.txt | tail -n 1 | awk '{ print $NF }')
    [[ -z "${interface_area}" ]] && interface_area=0
    [[ -z "${surface_area}" ]] && surface_area=0
    cd $output_dir
    rm -r drsasa_output
    prodigy $pdb_file > prodigy_output.txt
    affinity=$(tail -n 2 prodigy_output.txt | head -n 1 | awk '{print $NF}')
    kd25=$(tail -n 1 prodigy_output.txt | awk '{print $NF}')
    echo "$complex_name,$avg_n_models,$max_n_models,$num_contacts_with_max_n_models,$num_unique_contacts,$best_model_num,$best_pdockq,$best_plddt_avg,$best_pae_avg,$pTM,$ipTM,$ranking_confidence,$interface_area,$surface_area,$affinity,$kd25" >> $output_dir/updated_summary.csv
    cd $exec_dict
done

echo "Updated summary.csv with pTM, ipTM, and complex SASAs"
