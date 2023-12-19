#!/bin/bash

# Specifications for eddie

#$ -N individual_molecule_dist_corr_job                     # Name the job
#$ -pe sharedmem 1                              # Request 1 cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=128G                                 # Request 128G memory per cpu
#$ -cwd
#$ -l h_rt=240:00:00                              # Allocate a maximum of 8 hours to the job
source /etc/profile.d/modules.sh                # Initialise the environment modules

# Arguments:
# $1: path to file containing file set names
# $2: path to folder containing the bed files with distance info for reads containing at least 100 CpGs (gzipped)
# $3: path to folder containing bed file of genomic coords of reads
# $4: output directory for distance info
# $5: path to all_vs_all_individual_molecule_distance_dependent_correlations.R script
# $6: LLR threshold used to call methylation
# $7: max distance to consider correlations for
# $8: output path for correlations

# Input the list to be used as arguments 

F=`sed -n ${SGE_TASK_ID}p < $1`

# Load the required modules

module load roslin/R/4.1.0
module load igmm/apps/BEDTools/2.30.0

# Extract distance info for reads of interest

bedtools intersect -wa -wb -a $2/read_dist_info_$F.bed.gz -b $3 | awk 'OFS="\t" {if ($4==$11){print $1,$2,$3,$4,$5,$6,$7}}' | gzip > $4/read_dist_info_$F.bed.gz

# Make a directory to store the correlations

mkdir $8/$F

# Calculate distance-dependent correlations for each read in the individual reads

Rscript $5 $4/read_dist_info_$F.bed.gz $8/$F/ $6 $7

