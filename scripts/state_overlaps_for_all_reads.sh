#!/bin/bash

# Arguments
# $1: path to bed file containing genomic locations of all autosomal reads with at least 100 CpGs (gzip compressed)
# $2: path to output folder
# $3: path to chromHMM bed file
# $4: path to PMD bed file
# $5: path to non-PMD bed file

# Load required modules

module load igmm/apps/BEDTools/2.30.0



# Extract the overlap of reads in the group with chromHMM annotations, PMDs and non-PMDs

mkdir $2/specific_state_overlaps

# CHROMHMM

for i in $(awk '{print $4}' $3 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $3 | bedtools intersect -wao -a $1 -b - | awk 'OFS="\t" {print $4,$9}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/base_overlap_$i.txt.gz;
 grep -w $i $3 | bedtools intersect -wao -a $1 -b - | awk 'OFS="\t" {print $4,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/proportion_overlap_$i.txt.gz;
 done
 
 # Calculate the number of bases in the reads

num_bases=$(gzip -dc $1 | awk '{Tot+=$3-$2} END {print Tot}' )

# For all the reads in the group, calculate the % of bp cooresponding to each annotation

rm $2/chromHMM_overlap_percents.txt;
touch $2/chromHMM_overlap_percents.txt
for i in $(awk '{print $4}' $3 | sort -k1,1 | uniq | cat);
 do 
 gzip -dc $2/specific_state_overlaps/base_overlap_$i.txt.gz | awk -v tot_bases="${num_bases}" 'OFS="\t" {Tot+=$2} END {print 100*Tot/tot_bases}' >> $2/chromHMM_overlap_percents.txt;
 done
 
paste <(awk '{print $4}' $3 | sort -k1,1 | uniq | cat ) $2/chromHMM_overlap_percents.txt | awk 'OFS="\t" {print $1,$2}' > $2/chromHMM_overlap_percents_final.txt



### PMDS/NON-PMDS ###

# Extract the overlap and proportion overlap between the group reads and PMDs/non-PMDs

bedtools intersect -wao -a $1 -b $4 | awk 'OFS="\t" {print $4,$9}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/pmd_base_overlap.txt.gz
bedtools intersect -wao -a $1 -b $4 | awk 'OFS="\t" {print $4,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/pmd_proportion_overlap.txt.gz

bedtools intersect -wao -a $1 -b $5 | awk 'OFS="\t" {print $4,$9}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/nonpmd_base_overlap.txt.gz
bedtools intersect -wao -a $1 -b $5 | awk 'OFS="\t" {print $4,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 2 -o sum | gzip > $2/specific_state_overlaps/nonpmd_proportion_overlap.txt.gz

# Calculate the percentage of base pairs in the high RTS reads that overlap with PMDs/non-PMDs

pmd_overlap=$(bedtools intersect -a $4 -b $1 | awk '{Tot+=$3-$2} END {print Tot}')
nonpmd_overlap=$(bedtools intersect -a $5 -b $1 | awk '{Tot+=$3-$2} END {print Tot}') 
pmd_prop_overlap=$(echo "${pmd_overlap}/${num_bases}" | bc -l)
nonpmd_prop_overlap=$(echo "${nonpmd_overlap}/${num_bases}" | bc -l)

# Save to file
echo "PMD proportion overlap: $(echo "${pmd_prop_overlap}" | bc -l)" > $2/pmd_nonpmd_overall_prop_overlap.txt
echo "non-PMD proportion overlap: $(echo "${nonpmd_prop_overlap}" | bc -l)" >> $2/pmd_nonpmd_overall_prop_overlap.txt

