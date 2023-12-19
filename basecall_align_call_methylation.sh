#!/bin/bash

# Specifications for eddie

#$ -N base_align_call_meth			# Name the job
#$ -pe gpu-a100 1	           	# Request 1 gpu
#$ -R y                             # Reserve requested nodes as they become available
#$ -l h_vmem=64G					# Request 64G memory
#$ -cwd
#$ -l h_rt=48:00:00					# Allocate a maximum of 48 hours to the job
source /etc/profile.d/modules.sh	# Initialise the environment modules

# This script takes in nanopore fast5 files and outputs methylation tsv files

# Arguments:

# $1: File containing the a list of folder names where Nanopore fast5 files are stored
# $2: Path to the main folder containing the subfolders above
# $3: Path to reference genome
# $4: Name of conda environment where required software is installed
# $5: Path to main anaconda folder
# $6: Guppy configuration required
# $7: Output path for basecalled reads
# $8: Output path for aligned reads
# $9: Output path for methylation tsv files
# $10: Path to mtsv2bedGraph script
# $11: Log likelihood ratio threshold used to call methylation

# Input the list to be used as arguments 

F=`sed -n ${SGE_TASK_ID}p < $1`

# Load modules

module load anaconda
module load cuda/11.0.2
module load roslin/guppy/5.0.11-gpu

# Activate conda environment

source activate $4

LD_LIBRARY_PATH=$5/envs/$4/lib:$LD_LIBRARY_PATH 



# Perform basecall with guppy

guppy_basecaller --input_path $2/$F --save_path $7/$F --config $6 --num_callers 4 --device auto --compress_fastq

# Combine output files

cat $7/$F/pass/*.fastq.gz > $7/all_guppy_$F.fastq.gz



# Perform alignment of basecalled reads using Minimap2

# Create an index file that links read ids with their signal-level fast5 data.

nanopolish index -d $2/$F $7/all_guppy_$F.fastq.gz

# Next align to the reference genome (-a option produces a sam file, -x means we choose map-ont as the performance tuner).
# We pip the alignment output to samtools to create a sorted bam file. -T stores temporary files as tmp.$F.0000.bam and -o specifies output.

minimap2 -ax map-ont $3 $7/all_guppy_$F.fastq.gz --MD | samtools sort -T tmp.$F -o $8/alignments_$F.bam

# Filter out secondary alignments and supplementary alignments.
# The first command filters out supplementary alignments and the second filters out secondary alignments.

samtools view -h -F 0X800 $8/alignments_$F.bam | samtools view -h -F 0X100 | samtools view -bh -@ 8 - > $8/primary_$F.bam

samtools index $8/primary_$F.bam



# Call methylation using Nanopolish

nanopolish call-methylation -t 8 -r $7/all_guppy_$F.fastq.gz -b $8/primary_$F.bam -g $3 --progress > $9/meth_calls_$F.tsv

# Convert Nanopolish tsv to bed file

# Run python script to convert tsv file to bed file
# Remove non-canonical bases and sort
# (Then index the resulting bed file)

mkdir $9/bed_files
${10} -i $9/meth_calls_$F.tsv -g $3 -c ${11} -w 3 -v | grep -v -P 'chr.+_' | sort -k1,1 -k2,2n | bgzip > $9/bed_files/meth_calls_$F.bed.gz 
tabix -p bed $9/bed_files/meth_calls_$F.bed.gz 

