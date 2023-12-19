#!/bin/bash

# Specifications for eddie

#$ -N window_pop_stats_job                     # Name the job
#$ -pe sharedmem 1                              # Request 1 cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=128G                                 # Request 128G memory per cpu
#$ -cwd
#$ -l h_rt=48:00:00                              # Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh                # Initialise the environment modules

# Arguments:
# $1: path to chromosome size file for reference genome
# $2: genomic window size
# $3: main window directory
# $4: path to bed file containing bulk CpG methylation levels (as outputted by bulk_methylation.sh script)
# $5: path to window_population_statistics.R script

# Load the required modules

module load igmm/apps/BEDOPS/2.4.26
module load igmm/apps/BEDTools/2.30.0
module load roslin/R/4.1.0

# Create a bed file containing hg38 autosomes in {$2}kb windows

awk 'OFS="\t" {print $1, 0, $2}' $1 | grep -v -P 'chr.+_' | grep -v -E 'chrX|chrY|chrM' | bedops --chop $(($2*1000)) - | awk 'OFS="\t" {print $0, "window_"NR}' | sort -k1,1 -k2,2n > $3/${2}kb_autosome_windows.bed

# Create a file containing the window names

awk '{print $4}' $3/${2}kb_autosome_windows.bed > $3/windows.txt

# Create a folder to store bulk results

mkdir $3/population_level

# Extract CpGs that lie entirely within a single window

bedtools intersect -wa -wb -a $4 -b $3/${2}kb_autosome_windows.bed -f 1 | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$10}' > $3/population_level/population_CpG_summary.bed

Rscript $5 $3/population_level/population_CpG_summary.bed $3/population_level/window_population_stats.tsv