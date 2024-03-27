#!/bin/bash
#$ -N error_pipeline_handling       # {JOB_NAME}
#$ -V                 # environmental variables retain their values
#$ -R y               # reserve nodes as they become available
#$ -cwd               # current working directory
#$ -q gpu             # request a GPU
#$ -pe gpu-a100 2     # GPU resource
#$ -l h_vmem=200G     # memory resource per core
#$ -l h_rt=24:00:00   # wallclock time
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

# run localcolabfold
colabfold_batch --num-models 3 \
                --num-recycle 5 \
                --random-seed 42 \
                --templates \
                --msa-mode mmseqs2_uniref \
                --model-type alphafold2_multimer_v3 \
                ./protein_interactions_set/protein_chunk_004.fasta \
                ./protein_interactions_set/protein_chunk_004

# end
