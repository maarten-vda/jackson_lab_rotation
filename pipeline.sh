#!/bin/bash

### THIS SCRIPT COMBINES ALL ASPECTS OF THE ANALYSIS PIPELINE, GOING FROM UNIPROT IDS OR GO ANNOTATIONS TO A RANKED LIST OF POTENTIAL INTERACTORS

##Eddie qsub flags
#$ -N protein_interactors_set       # {JOB_NAME}
#$ -V                 # environmental variables retain their values
#$ -R y               # reserve nodes as they become available
#$ -cwd               # current working directory
#$ -q gpu             # request a GPU
#$ -pe gpu-a100 2     # GPU resource
#$ -l h_vmem=300G     # memory resource per core
#$ -l h_rt=96:00:00   # wallclock time
#$ -m ae
#$ -M s2037423@ed.ac.uk

# dependencies
source /exports/applications/support/set_cuda_visible_devices.sh
source /etc/profile.d/modules.sh
module load phys/compilers/gcc/11.2.0
module load igmm/apps/openssl/3.0.5
module load igmm/apps/python/3.10.6
module load cuda/12.1.1
export PATH="/exports/igmm/eddie/ajackson-wrkgrp/maarten/software/localcolabfold/colabfold-conda/bin:$PATH"

# Initialize default variables
scripts_dict=$PWD
target=""
interactors="../data/uniprot_ids/protein_interactions_stratified_set.txt"
go_term=""
output_folder="./testing"
num_models="3"
num_recycle="5"
mode="all_vs_all"
msa_mode="mmseqs2_uniref"
model_type="alphafold2_multimer_v3"
templates="On"
model_random_seed="42"

# Function to display usage information
usage() {
    echo "Usage: pipeline.sh --target <Target Protein UniProt ID> \
        (--interactors <Interactors UniProt ID> XOR --go_term <GO Term>) \
        --output_folder <Output Folder> [--mode <Mode (pulldown, all_vs_all, custom)> \
        --num-models <Number of models per interaction to fold>] \
        [--num-recycle <Number of recycles per model>] \
        [--msa-mode <MSA mode (default: mmseqs2_uniref)>] \
        [--model-type <Model Type (default: alphafold2_multimer_v3)>] \
        [--templates <On or Off (default: On)>] \
        [--model-random-seed <Model Random Seed (default: 42)>]"
    echo "Display help information: pipeline.sh --help"
    exit 1
}

help() {
    echo "Usage: $0 --target <Target Protein UniProt ID> \
        (--interactors <Interactors UniProt ID> XOR --go_term <GO Term>) \
        --output_folder <Output Folder> [--mode <Mode (pulldown, all_vs_all, custom)>] \
        [--num-models <Number of models per interaction to fold>] \
        [--num-recycle <Number of recycles per model>] \
        [--msa-mode <MSA mode (default: mmseqs2_uniref)>] \
        [--model-type <Model Type (default: alphafold2_multimer_v3)>] \
        [--templates <Templates (default: On)>] \
        [--model-random-seed <Model Random Seed (default: 42)>]"
    echo "Display help information: $0 --help"
    echo "Description:"
    echo "  This pipeline is for predicting complexes using AlphaFold-Multimer, specifically the Colabfold implementation."
    echo "  The following options are available:"
    echo "	--target"
    echo "	A text file containing the UniProt ID of the target protein. Required in pulldown and all_vs_all mode"
    echo "      --interactors"
    echo "	A text file containing the line separated UniProt IDs of potential interactors to screen through. Required if not using Gene Ontology sets"
    echo "      --go_term"
    echo "	The Gene Ontology ID to screen through. Either a list of potential interactors or a Gene Ontology ID is required. Can either be in the format \"GO:000xxxx\" or \"000xxxx\" (Including quotation marks)"
    echo "      --output_folder"
    echo "	A filepath to the desired directory for output files. If no directory is provided then the current working directory is used. If a directory is supplied that does not exist it will be created."
    echo "      --mode"
    echo "	The mode in which the complex structure prediction will be run (Default: pulldown). Pulldown mode compares every interacting protein to the target protein. All_vs_all mode compares every protein in the target and interactors list with every other protein in the target and interactors list. Custom mode does not provide any comparison, instead it requires the input specified by --interactors to be a multi-line FASTA file where complexes are represented by : separated fasta files on the same line (eg. MYAAPQG:MYAAPQG:MYAAPQG)."
    echo "      --num-models"
    echo "	The number of models which AlphaFold should generate (Default: 3)"
    echo "      --num-recycle"
    echo "	The number of recycles per AlphaFold model (Default: 5)"
    echo "      --msa-mode"
    echo "	The MSA mode AlphaFold should use (Default: mmseqs2_uniref)"
    echo "      --model-type"
    echo "	The AlphaFold model to be used (Default: alphafold2_multimer_v3)"
    echo "      --templates"
    echo "	Whether AlphaFold should use Templates. Can be either On or Off (Default: Templates on)"
    echo "      --model-random-seed"
    echo "	The model random seed, which is an arbitrary number. Keep the same for consistency across runs. (Default: 42)"
    echo "	--help"
    echo "	Display this message."
    exit 1
}

# Check for --help option
if [[ "$#" -eq 1 && "$1" == "--help" ]]; then
    help
fi

# Parse command line options
while [[ "$#" -gt 0 ]]; do
    case $1 in
        --target)
            target="$2"
            shift 2
            ;;
        --interactors)
            interactors="$2"
            shift 2
            ;;
        --go_term)
            go_term="$2"
            shift 2
            ;;
        --output_folder)
            output_folder="$2"
            shift 2
            ;;
	--mode)
	    mode="$2"
	    shift 2
	    ;;
        --num-models)
            num_models="$2"
            shift 2
            ;;
        --num-recycle)
            num_recycle="$2"
            shift 2
            ;;
        --msa-mode)
            msa_mode="$2"
            shift 2
            ;;
        --model-type)
            model_type="$2"
            shift 2
            ;;
        --templates)
            templates="$2"
            shift 2
            ;;
        --model-random-seed)
            model_random_seed="$2"
            shift 2
            ;;
        *)
            echo "Unknown option: $1"
            usage
            ;;
    esac
done

# Check that either --interactors or --go_term is provided, but not both or neither
if [[ -z "$interactors" && -z "$go_term" ]] || [[ -n "$interactors" && -n "$go_term" ]]; then
    echo "Error: Please provide either --interactors <Interactors UniProt ID> or --go_term <GO Term UniProt ID>, but not both or neither."
    usage
fi

# Make --target essential only if --mode is not 'custom'
if [[ "$mode" == "pulldown" && -z "$target" ]]; then
    echo "Error: --target is essential when --mode is 'pulldown'."
    usage
fi

# Make the output folder an absolute filepath and not relative
output_folder="${PWD}/${output_folder}"
# Check if the output directory exists, and create it if not
if [ ! -d "$output_folder" ]; then
    mkdir -p "$output_folder"
fi


# Construct the name for this run, which is just the target name followed by either the interactors or gene ontology filename
# This part converts interactors or the gene ontology part to FASTA
if [[ "$target" == *.* ]]; then
    targetbasename=$(basename "$target")
    targetname="${targetbasename%.*}"
else 
    targetname="$target"
fi
if [ -n "$go_term" ]; then
    interactor_name=$go_term
    # Check if colon is present in the variable
    if [[ $interactor_name == *:* ]]; then
        interactor_name="${interactor_name#*:}"  # Extract part after the colon
    fi
    echo "Converting GO:${interactor_name} to fasta..."
    filename="${targetname}_GO${interactor_name}"
    bash GO_terms.sh $go_term ${output_folder}/${filename}_go_uniprot_ids.txt
    python uniprot_to_fasta.py -i ${output_folder}/${filename}_go_uniprot_ids.txt -o ${output_folder}/${filename}_go_fastas.fasta
    binding_partners=${output_folder}/${filename}_go_fastas.fasta
elif [ -n "$interactors" ]; then
    if [[ "$interactors" == *.* ]]; then
        interactorbasename=$(basename "$interactors")
        interactor_name="${interactorbasename%.*}"
    else
        echo "Error: interactors must be a file or GO annotation"
        exit 1
    fi
    if [[ -z "$targetname" ]]; then
        filename="$interactor_name"
    else
        filename="${targetname}_${interactor_name}"
    fi
    if [[ ! "$interactors" =~ \.fasta$ ]] && [[ "$interactors" == *.* ]]; then
        echo "Converting ${interactors} to fasta..."
        python uniprot_to_fasta.py -i ${interactors} -o ${output_folder}/${filename}_interactors.fasta
        binding_partners=${output_folder}/${filename}_interactors.fasta
    elif [[ "$interactors" == *.fasta ]]; then
        binding_partners=$interactors
    fi
else
    echo "Error: Either go_term or interactors must have a value."
    exit 1
fi



#Convert target to fasta, and generate the different combinations to enter into Colabfold
# Check if $target does not have the .fasta extension but does have an extension, which implies it is a text file with a uniprot ID, and convert to fasta
# This part processes the target variable depending on if its a fasta, uniprot TXT file, or uniprot string
if [[ ! "$target" =~ \.fasta$ ]] && [[ "$target" == *.* ]]; then
    # Execute the following commands
    python uniprot_to_fasta.py -i "$target" -o "${output_folder}/${filename}_target.fasta"
    target_fasta="${output_folder}/${filename}_target.fasta"
# If target is in the fasta format then just continue
elif [[ "$target" == *.fasta ]]; then
    # If $target already has the .fasta extension, use it as is
    target_fasta="$target"
# If target has no filetype then it must be the uniprot ID as text
elif [[ "$mode" != "custom" && "$target" != *.* ]]; then
    echo "Converting ${target} to fasta..."
    echo $target > "${output_folder}/${filename}_uniprot.txt"
    python uniprot_to_fasta.py -i "${output_folder}/${filename}_uniprot.txt" -o "${output_folder}/${filename}_target.fasta"
    target_fasta="${output_folder}/${filename}_target.fasta"
fi




# Generate the AF-M inputs depending on the mode
echo "Generating ${mode} mode multimer representations..."
if [[ "$mode" == "pulldown" ]]; then
    bash pulldown.sh $target_fasta $binding_partners ${output_folder}/${filename}_pulldown.fasta
    colabfold_input="${output_folder}/${filename}_pulldown.fasta"
elif [[ "$mode" == "all_vs_all" ]]; then
    cat $target_fasta > ${output_folder}/${filename}_combined.fasta
    cat $binding_partners >> ${output_folder}/${filename}_combined.fasta
    python list_to_pairwise.py ${output_folder}/${filename}_combined.fasta ${output_folder}/${filename}_all_vs_all.fasta
    colabfold_input="${output_folder}/${filename}_all_vs_all.fasta"
elif [[	"$mode"	== "custom" ]]; then
    echo "Initialising custom mode"
fi



# run localcolabfold
echo -e "Running localcolabfold with the following options: \n --num-models $num_models \n --num-recycle $num_recycle \
    \n --random-seed $model_random_seed \n --msa-mode $msa_mode \n --model-type $model_type \
    \n input $colabfold_input \n output folder ${output_folder}/colabfold_predictions"

if [[ "$mode" == "pulldown" ]]; then
    colabfold_batch --num-models $num_models --num-recycle $num_recycle \
    --random-seed $model_random_seed --msa-mode $msa_mode --model-type $model_type \
    $colabfold_input ${output_folder}/colabfold_predictions
elif [[ "$mode" == "all_vs_all" ]]; then
    colabfold_batch --num-models $num_models --num-recycle $num_recycle \
    --random-seed $model_random_seed --msa-mode $msa_mode --model-type $model_type \
    $colabfold_input ${output_folder}/colabfold_predictions
elif [[ "$mode" == "custom" ]]; then
    colabfold_batch --num-models $num_models --num-recycle $num_recycle \
    --random-seed $model_random_seed --msa-mode $msa_mode --model-type $model_type \
    $interactors ${output_folder}/colabfold_predictions
fi



# analyse results using DONSON paper colabfold_analysis.py script, then add ipTM and pTM using combine_metrics.py, then score interactions and seperate interactions into those that pass threshold and those that dont
cd ${output_folder}
python ${scripts_dict}/colabfold_analysis.py ${output_folder}/colabfold_predictions
python ${scripts_dict}/combine_metrics.py colabfold_predictions_analysis/interfaces.csv ${output_folder}/colabfold_predictions combined_metrics.csv
bash ${scripts_dict}/improved_update_summary.sh ${output_folder}
# Here is where the script to graph the interactions (graph_interactions.py) should go
cd $scripts_dict
