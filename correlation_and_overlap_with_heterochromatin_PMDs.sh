#!/bin/bash

# Arguments:

# $1: path to chromHMM bed file
# $2: path to PMD location bed file
# $3: path to output folder for single-read info
# $4: path to output folder for window info
# $5: bed file containing genomic coordinates of autosomal reads with at least 100 CpGs
# $6: tsv file containing stats for autosomal reads with at least 100 CpGs (gzip compressed)
# $7: path to bed file for all windows containing correlation binned group info
# $8: path to tsv file containing mean stats for each window (gzip compressed)

# Load required modules

module load igmm/apps/BEDTools/2.30.0

# EXTRACT THE INFORMATION REQUIRED TO CALCULATE THE CORRELATION BETWEEN CORRELATION AND OVERLAP WITH HETEROCHROMATIN/PMD AT THE WINDOW/READ LEVEL


# HETEROCHROMATIN READ LEVEL

# Extract the overlap and proportion overlap between all reads with min 100 CpGs and heterochromatin
# Combine this with correlation info

grep -w "heterochromatin" $1 | bedtools intersect -wao -a $5 -b - | awk 'OFS="\t" {print $4,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | join -1 1 -2 1 - <(gzip -dc $6 | sort -k1,1) | awk 'OFS="\t" {print $1,$2,$6}' | gzip > $3/min_100_CpGs_read_hetero_overlap_correlation.tsv.gz

# HETEROCHROMATIN WINDOW LEVEL

# Extract the proportion of each window that overlaps with heterochromatin 

grep 'heterochromatin' $1 | bedtools intersect -wao -a $7 -b - | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | bedtools groupby -i - -g 1,2,3,4,5 -c 10 -o sum | awk 'OFS="\t" {print $4,$5,$6/($3-$2)}' > $4/heterochromatin_proportions.tsv

# Combine the heterochromatin proportion info with the correlation window info

gzip -dc $8 | sort -k1,1 | join -1 1 -2 1 <(sort -k1,1 $4/heterochromatin_proportions.tsv) - | awk 'OFS="\t" {print $1,$3,$7}' > $4/windows_heterochromatin_correlation.tsv



# PMD READ LEVEL

# Extract the overlap and proportion overlap between all reads with min 100 CpGs and PMDs
# Combine this with correlation info

bedtools intersect -wao -a $5 -b $2 | awk 'OFS="\t" {print $4,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | join -1 1 -2 1 - <(gzip -dc $6 | sort -k1,1) | awk 'OFS="\t" {print $1,$2,$6}' > $3/min_100_CpGs_read_pmd_overlap_correlation.tsv.gz

# PMD WINDOW LEVEL

# Extract the proportion that each window overlaps with PMDs

bedtools intersect -wao -a $7 -b $2 | sort -k1,1 -k2,2n -k3,3n -k4,4 -k5,5 | bedtools groupby -i - -g 1,2,3,4,5 -c 10 -o sum | awk 'OFS="\t" {print $4,$5,$6/($3-$2)}' > $4/pmd_proportions.tsv

# Combine the PMD proportion info with the correlation window info

gzip -dc $8 | sort -k1,1 | join -1 1 -2 1 <(sort -k1,1 $4/pmd_proportions.tsv) - | awk 'OFS="\t" {print $1,$3,$7}' > $4/windows_pmd_correlation.tsv


