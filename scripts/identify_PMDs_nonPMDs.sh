#!/bin/bash

# Specifications for eddie

#$ -N pmd_detection_job			# Name the job
#$ -pe sharedmem 4		           	# Request 4 cpus
#$ -R y                             # Start reserving requested nodes
#$ -l h_vmem=8G 					# Request 8G memory per cpu
#$ -cwd
#$ -l h_rt=8:00:00					# Allocate a maximum of 8 hours to the job
source /etc/profile.d/modules.sh	# Initialise the environment modules

# Arguments:
# $1: path to directory where PMD and non-PMD folders will be created
# $2: path to methylation tsv file
# $3: path to chromosome size file
# $4: path to bed file containing regions to exclude

# Load required modules

module load igmm/apps/methpipe/5.0.0
module load igmm/apps/BEDTools/2.30.0 

# Make directories to work in

mkdir $1/PMD_info
mkdir $1/nonPMD_info

# Use methylation tsv to create a file in format required my methpipe

awk -F '\t' 'OFS="\t" {if($7!="NA"){print $1, $2, "+", "CpG", $7, $5}}' $2 > $1/PMD_info/methpipe_input.meth

# Run the pmd detector

pmd -i 1000 -o $1/PMD_info/locations.pmd $1/PMD_info/methpipe_input.meth

# The output from methpipe does not include the G within the final CpG as part of the PMD
# We create a bed file to rectify this

awk 'OFS="\t" {print $1,$2,$3+1,$4}' $1/PMD_info/locations.pmd > $1/PMD_info/locations.pmd.bed



# Identify non-PMDs

bedtools complement -i $1/PMD_info/locations.pmd.bed -g $3 > $1/nonPMD_info/locations.nonpmd.bed

# Exclude regions that overlap with reference sequence gaps, centromeres and generally poorly sequenced parts of the genome
# Exclude chr X, Y and M
# Exclude PMDs/nonPMDs that are smaller than 200kb (note all PMDs are >500kb anyway)
# Add name to PMDs/non-PMDs

bedtools subtract -a $1/PMD_info/locations.pmd.bed -b $4 | grep -v -E 'chrX|chrY|chrM' | awk '{if($3-$2>=200000) {print $0}}' | awk 'OFS="\t" {print $1,$2,$3,"PMD:"NR-1}' > $1/PMD_info/locations.pmd.excl.bed
bedtools subtract -a $1/nonPMD_info/locations.nonpmd.bed -b $4 | grep -v -E 'chrX|chrY|chrM' | awk '{if($3-$2>=200000) {print $0}}' | awk 'OFS="\t" {print $1,$2,$3,"nonPMD:"NR-1}' > $1/nonPMD_info/locations.nonpmd.excl.bed

