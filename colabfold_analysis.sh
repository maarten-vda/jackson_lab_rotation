#!/bin/bash
#$ -N analyse_results       # {JOB_NAME}
#$ -V                 # environmental variables retain their values
#$ -R y               # reserve nodes as they become available
#$ -cwd               # current working directory
#$ -l h_vmem=100G     # memory resource per core
#$ -l h_rt=48:00:00   # wallclock time
#$ -m ae
#$ -M s2037423@ed.ac.uk

# dependencies
source /etc/profile.d/modules.sh
module load igmm/apps/python/3.10.6


# run colabfold_analysis.py
python colabfold_analysis.py ../data/output/interactions_test


# end
