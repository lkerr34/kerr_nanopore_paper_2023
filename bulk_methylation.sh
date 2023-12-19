#!/bin/bash

# Specifications for eddie

#$ -N summarise_meth_job			# Name the job
#$ -pe sharedmem 1		           	# Request 4 cpus
#$ -R y                             # Start reserving requested nodes
#$ -l h_vmem=128G 					# Request 128G memory per cpu
#$ -cwd
#$ -l h_rt=12:00:00					# Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh	# Initialise the environment modules

# Arguments:
# $1: path to Nanopolish methylation tsv
# $2: path to bulk output directory
# $3: path to anaconda directory
# $4: name of anaconda environment
# $5: path to calculate_methylation_frequency.py script
# $6: bed file containing genomic regions to exclude from analysis 
# $7: log-likelihood threshold to use for methylation calls

# Load modules

module load anaconda

# Activate conda environment

source activate $4

LD_LIBRARY_PATH=$3/envs/$4/lib:$LD_LIBRARY_PATH #(add a folder to path which stores a file needed)

# Run script to summarise methylation 
# Remove non-canonical chromosomes from output

$5 $1 -c $7 -s| grep -v -P 'chr.+_' > $2/summary_all_meth_calls.tsv



# Convert summary CpG tsv file to bed format

awk 'OFS="\t" {print $1,$2,$3+2,$5,$6,$7,$8,$9}' $2/summary_all_meth_calls.tsv > $2/summary_all_meth_calls.bed

# Exclude CpGs where coverage is less than 10
# Exclude regions that overlap with reference sequence gaps, centromeres and generally poorly sequenced parts of the genome
# Exclude chr X, Y and M chromosomes

awk 'OFS="\t" {if($4>=10) {print $1,$2,$3,$4,$5,$6}}' $2/summary_all_meth_calls.bed | bedtools subtract -a - -b $6 | grep -v -E 'chrX|chrY|chrM' > $2/summary_meth_calls_excl.bed
