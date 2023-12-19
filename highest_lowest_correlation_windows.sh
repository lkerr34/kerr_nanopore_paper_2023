#!/bin/bash

# Specifications for eddie

#$ -N window_mean_stats_job                     # Name the job
#$ -pe sharedmem 1                              # Request 1 cpu
#$ -R y                                         # Start reserving requested nodes
#$ -l h_vmem=32G                                 # Request 32G memory per cpu
#$ -cwd
#$ -l h_rt=48:00:00                              # Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh                # Initialise the environment modules

# Arguments:
# $1: path to calculate_mean_window_stats.R script
# $2: path to the bed file with read-level stats for each window
# $3: path to single-molecule window directory
# $4: path to bed file with genomic windows
# $5: path to bed file containing bulk methylation levels

# Load the required modules

module load roslin/R/4.1.0
module load igmm/apps/BEDTools/2.30.0

Rscript $1 $2 $3/mean_window_stats.tsv

# Extract locations of all windows for which there was sufficient data to perform the single-molecule analysis
# Also save the number of these windows as a variable

gzip -dc $3/mean_window_stats.tsv.gz | awk '{print $1}' | grep -wf - $4 > $3/windows_analysed.bed
num_windows=$(wc -l $3/windows_analysed.bed | awk '{print $1}')

# Create bedgraphs with the mean methylation and correlation for each window

gzip -dc $3/mean_window_stats.tsv.gz | sort -k1,1 | join -1 1 -2 4 - <(sort -k4,4 $3/windows_analysed.bed) | awk 'OFS="\t" {print $8,$9,$10,$2}' > $3/windows_mean_meth.bedGraph
gzip -dc $3/mean_window_stats.tsv.gz | sort -k1,1 | join -1 1 -2 4 - <(sort -k4,4 $3/windows_analysed.bed) | awk 'OFS="\t" {print $8,$9,$10,$5}' > $3/windows_correlation.bedGraph


# Extract the 10% of windows with highest correlation and the 10% of windows with the lowest correlation
# Extract locations of these windows and save as bed files

gzip -dc $3/mean_window_stats.tsv.gz | awk 'NR==1; NR > 1 {print $0 | "sort -nr -k5,5"}' | head -$(((${num_windows}+10-1)/10)) | gzip - > $3/top_10_percent_correlation.tsv.gz
gzip -dc $3/mean_window_stats.tsv.gz | awk 'NR==1; NR > 1 {print $0 | "sort -n -k5,5"}' | head -$(((${num_windows}+10-1)/10)) | gzip - > $3/bottom_10_percent_correlation.tsv.gz

gzip -dc $3/top_10_percent_correlation.tsv.gz | awk '{print $1}' | grep -wf - $4 > $3/top_10_percent_correlation_windows.bed
gzip -dc $3/bottom_10_percent_correlation.tsv.gz | awk '{print $1}' | grep -wf - $4 > $3/bottom_10_percent_correlation_windows.bed

# Identify the methylation levels associated with the 10% of windows with the highest correlation and the 10% of windows with the lowest correlation

bedtools intersect -wa -a $5 -b $3/top_10_percent_correlation_windows.bed > $3/least_disordered_correlation_CpG_summary.bed
bedtools intersect -wa -a $5 -b $3/bottom_10_percent_correlation_windows.bed > $3/most_disordered_correlation_CpG_summary.bed

# Identify the methylation levels associated with all CpGs except those in the 10% of windows with the highest correlation 
# Similarly for the 10% of windows with the lowest correlation

bedtools subtract -a $5 -b $3/top_10_percent_correlation_windows.bed > $3/pop_CpG_summary_excl_least_disordered.bed
bedtools subtract -a $5 -b $3/bottom_10_percent_correlation_windows.bed > $3/pop_CpG_summary_excl_most_disordered.bed




