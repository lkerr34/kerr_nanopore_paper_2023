#!/bin/bash

# Specifications for eddie

#$ -N read_distance_info_job		# Name the job
#$ -pe sharedmem 1		           	# Request 1 cpus
#$ -R y                             # Start reserving requested nodes
#$ -l h_vmem=4G 					# Request 4G memory per cpu
#$ -cwd
#$ -l h_rt=1:00:00					# Allocate a maximum of 1 hour to the job
source /etc/profile.d/modules.sh	# Initialise the environment modules

# Calculate the distance between CpGs in individual reads

# Arguments:
# $1: path to file containing file set names
# $2: path to folder containing the methylation bed files
# $3: path to store output relating to all reads
# $4: path to bed file containing the genomic locations and read names of all reads with minimum 100 CpGs (gzip compressed)
# $5: path to store output relating to reads with minimum 100 CpGs

# Input the list to be used as arguments 

F=`sed -n ${SGE_TASK_ID}p < $1`

# Extract the distances between CpGs from the bed files (only keep canonical autosomal chromosomes)

gzip -dc $2/meth_calls_$F.bed.gz | grep -v -E 'chrX|chrY|chrM|chr.+_' | awk '{print $5}' | sed 's/m/,/g' - | sed 's/u/,/g' | sed 's/x/,/g'| sed 's/.$//' | sed 's/^.//' | sed 's/^.//' | sed 's/^$/NA/g' > $3/distances_$F.txt

# Extract all items except last from log likelihood ratios in bed files

gzip -dc $2/meth_calls_$F.bed.gz | grep -v -E 'chrX|chrY|chrM|chr.+_' | awk -F '\t' '{print $6}' | awk -F ',' 'OFS="," {if (NF>1) {$NF=""; print $0} else {print "NA"}}' | sed 's/,$//'  > $3/log_lik_1_$F.txt

# Extract all items except last from log likelihood ratios in bed files

gzip -dc $2/meth_calls_$F.bed.gz | grep -v -E 'chrX|chrY|chrM|chr.+_' | awk -F '\t' '{print $6}' | awk -F ',' 'OFS="," {if (NF>1) {$1=""; print $0} else {print "NA"}}' | sed 's/^,//' > $3/log_lik_2_$F.txt



# Create a new data frame containing the chromosome, start/end position, read name, log likelihood info and context info for neighbouring CpGs extracted above and the distances between CpGs in the read

gzip -dc $2/meth_calls_$F.bed.gz | awk 'OFS="\t" {print $1, $2, $3, $4}' | paste -d '\t' - $3/log_lik_1_$F.txt $3/log_lik_2_$F.txt $3/distances_$F.txt | awk -F '\t' '{num_distances=split($5,a,","); num_llr=split($7,a,","); if(num_distances==num_llr) {print $0}}' | gzip > $3/read_dist_info_$F.bed.gz


# Remove intermediate files

rm $3/distances_$F.txt $3/log_lik_1_$F.txt $3/log_lik_2_$F.txt



# Extract information for all reads containing at least 100 CpGs

join -t $'\t' -1 4 -2 1 <(gzip -dc $3/read_dist_info_$F.bed.gz | sort -k4,4) <(gzip -dc $4 | awk '{print $4}' | sort -k1,1) | awk 'OFS="\t" {print $2,$3,$4,$1,$5,$6,$7}' | gzip > $5/read_dist_info_$F.bed.gz

