#!/bin/bash

# Specifications for eddie

#$ -N single_read_stats_job                     # Name the job
#$ -pe sharedmem 1                              # Request 1 cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=2G                                 # Request 2G memory per cpu
#$ -cwd
#$ -l h_rt=2:00:00                              # Allocate a maximum of 2 hours to the job
source /etc/profile.d/modules.sh                # Initialise the environment modules

# Arguments:
# $1: path to file containing file set names
# $2: path to single_read_statistics.R script
# $3: path to the folder containing methylation bed files
# $4: path to desired output folder for read stats
# $5: path to genomic window bed file
# $6: path to output folder for window single-read info
# $7: log-likelihood ratio threshold used to call methylation

F=`sed -n ${SGE_TASK_ID}p < $1`

# Load the required modules

module load roslin/R/4.1.0
module load igmm/apps/BEDTools/2.30.0

# Install the packages required by R to read in the data

#R
#install.packages("readr",repos="http://cran.us.r-project.org")
#q()

# Run Rscript

Rscript $2 $3/meth_calls_$F.bed.gz $4/stats_$F.bed $7

# Use bedtools to identify the reads that are in each genomic window
# Save read in a bed file that stores the name of the genomic window

gzip -dc $4/stats_$F.bed.gz | sed 1d - | bedtools intersect -bed -wa -wb -a - -b $5 -f 1 | awk 'OFS="\t" {print $1,$2,$3,$4,$8,$9,$10,$11,$12,$13,$14,$15,$19}' | gzip - > $6/stats_$F.bed.gz
