#!/bin/bash

# Arguments:
# $1: bed file containing all windows
# $2: bed file containing lowest correlation windows
# $3: bed file containig highest correlation windows
# $4: gene bed file (gzip compressed)
# $5: CGI bed file
# $6: chromHMM bed file
# $7: PMD bed file
# $8: non-PMD bed file
# $9: path to output folder

# Load required modules

module load igmm/apps/BEDTools/2.30.0

# Make directory to store proportion outputs

mkdir $9/specific_state_overlaps



### GENE AND CGI ENRICHMENT ###

# Calculate the number of bases in all windows
all_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $1)

# Calculate the overlap and proportion overlap of all windows with genes and CGIs 
all_gene_overlap=$(gzip -dc $4 | bedtools merge -i - | bedtools intersect -a $1 -b - | awk '{Tot+=$3-$2} END {print Tot}')
all_CGI_overlap=$(bedtools intersect -a $1 -b $5 | awk '{Tot+=$3-$2} END {print Tot}')
all_gene_prop_overlap=$(echo "${all_gene_overlap}/${all_bases}" | bc -l)
all_CGI_prop_overlap=$(echo "${all_CGI_overlap}/${all_bases}" | bc -l)



# Save list of the proportion overlaps (required for testing significance of enrichments)
gzip -dc $4 | bedtools merge -i - | bedtools intersect -wao -a $1 -b - | bedtools groupby -i - -g 1,2,3,4 -c 8 -o sum | awk '{print $5/($3-$2)}' | gzip > $9/specific_state_overlaps/window_proportion_overlap_genes.txt.gz
bedtools intersect -wao -a $1 -b $5 | bedtools groupby -i - -g 1,2,3,4 -c 9 -o sum | awk '{print $5/($3-$2)}' | gzip > $9/specific_state_overlaps/window_proportion_overlap_CGIs.txt.gz



# Calculate the number of bases in the lowest correlation windows
lowest_correlation_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $2)

# Calculate the overlap and proportion overlap of the lowest correlation windows with genes and CGIs
lowest_correlation_gene_overlap=$(gzip -dc $4 | bedtools merge -i - | bedtools intersect -a $2 -b - | awk '{Tot+=$3-$2} END {print Tot}')
lowest_correlation_CGI_overlap=$(bedtools intersect -a $2 -b $5 | awk '{Tot+=$3-$2} END {print Tot}')
lowest_correlation_gene_prop_overlap=$(echo "${lowest_correlation_gene_overlap}/${lowest_correlation_bases}" | bc -l)
lowest_correlation_CGI_prop_overlap=$(echo "${lowest_correlation_CGI_overlap}/${lowest_correlation_bases}" | bc -l)

# Save list of the proportion overlaps (required for testing significance of enrichments)
gzip -dc $4 | bedtools merge -i - | bedtools intersect -wao -a $2 -b - | bedtools groupby -i - -g 1,2,3,4 -c 8 -o sum | awk '{print $5/($3-$2)}' | gzip > $9/specific_state_overlaps/lowest_correlation_proportion_overlap_genes.txt.gz
bedtools intersect -wao -a $2 -b $5 | bedtools groupby -i - -g 1,2,3,4 -c 9 -o sum | awk '{print $5/($3-$2)}' | gzip > $9/specific_state_overlaps/lowest_correlation_proportion_overlap_CGIs.txt.gz



# Repeat this analysis for the highest correlation windows

# Calculate the number of bases in the highest correlation windows
highest_correlation_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3)

# Calculate the overlap and proportion overlap of the highest correlation windows with genes and CGIs
highest_correlation_gene_overlap=$(gzip -dc $4 | bedtools merge -i - | bedtools intersect -a $3 -b - | awk '{Tot+=$3-$2} END {print Tot}')
highest_correlation_CGI_overlap=$(bedtools intersect -a $3 -b $5 | awk '{Tot+=$3-$2} END {print Tot}')
highest_correlation_gene_prop_overlap=$(echo "${highest_correlation_gene_overlap}/${highest_correlation_bases}" | bc -l)
highest_correlation_CGI_prop_overlap=$(echo "${highest_correlation_CGI_overlap}/${highest_correlation_bases}" | bc -l)

# Save list of the proportion overlaps (required for testing significance of enrichments)
gzip -dc $4 | bedtools merge -i - | bedtools intersect -wao -a $3 -b - | bedtools groupby -i - -g 1,2,3,4 -c 8 -o sum | awk '{print $5/($3-$2)}' | gzip > $9/specific_state_overlaps/highest_correlation_proportion_overlap_genes.txt.gz
bedtools intersect -wao -a $3 -b $5 | bedtools groupby -i - -g 1,2,3,4 -c 9 -o sum | awk '{print $4,$5/($3-$2)}' | gzip > $9/specific_state_overlaps/highest_correlation_proportion_overlap_CGIs.txt.gz



# Calculate enrichments and save to file

echo "gene enrichment: $(echo "${lowest_correlation_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $9/low_correlation_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${lowest_correlation_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $9/low_correlation_gene_CGI_enrichments.txt

echo "gene enrichment: $(echo "${highest_correlation_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $9/high_correlation_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${highest_correlation_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $9/high_correlation_gene_CGI_enrichments.txt



### CHROMHMM ENRICHMENTS ###

# For each chromatin state extract the overlap and proportion overlap with all windows, low correlation windows and high correlation windows

for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat)
do
 grep -w $i $6 | bedtools intersect -wao -a $1 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/window_proportion_overlap_$i.txt.gz;
 grep -w $i $6 | bedtools intersect -wao -a $2 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/bottom_10_percent_correlation_proportion_overlap_$i.txt.gz;
 grep -w $i $6 | bedtools intersect -wao -a $3 -b - | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/top_10_percent_correlation_proportion_overlap_$i.txt.gz;
done

# Extract the overlap between all windows and the chromHMM annotations in one file
bedtools intersect -a $6 -b $1 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $9/window_overlap_annotations.bed

# For all the windows analysed, calculate the % of bp cooresponding to each annotation

rm $9/all_window_overlap_chromHMM_percents.txt;
touch $9/all_window_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $9/window_overlap_annotations.bed | awk -v a="${all_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $9/all_window_overlap_chromHMM_percents.txt;
 done



# Repeat for the lowest correlation windows

# Extract the overlap between the lowest correlation windows and the chromHMM annotations
bedtools intersect -a $6 -b $2 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $9/bottom_10_percent_correlation_overlap_annotations.bed

# For the lowest correlation windows, calculate the % of bp cooresponding to each annotation
rm $9/bottom_10_percent_correlation_overlap_chromHMM_percents.txt
touch $9/bottom_10_percent_correlation_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $9/bottom_10_percent_correlation_overlap_annotations.bed | awk -v a="${lowest_correlation_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $9/bottom_10_percent_correlation_overlap_chromHMM_percents.txt;
 done



# Repeat for the highest correlation windows

# Extract the overlap between the highest correlation windows and the chromHMM annotations
bedtools intersect -a $6 -b $3 | awk 'OFS="\t" {print $1,$2,$3,$4,$3-$2}' > $9/top_10_percent_correlation_overlap_annotations.bed

# For the highest correlations windows, calculate the % of bp cooresponding to each annotation
rm $9/top_10_percent_correlation_overlap_chromHMM_percents.txt
touch $9/top_10_percent_correlation_overlap_chromHMM_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w $i $9/top_10_percent_correlation_overlap_annotations.bed | awk -v a="${highest_correlation_bases}" 'OFS="\t" {Tot+=$3-$2} END {print $4,100*Tot/a}' >> $9/top_10_percent_correlation_overlap_chromHMM_percents.txt;
 done



# NOTE THAT ENRICHMENTS CAN THEN BE OBTAINED BY DIVIDING VALUES IN bottom_10_percent_correlation_overlap_chromHMM_percents.txt
# and top_10_percent_correlation_overlap_chromHMM_percents.txt BY VALUES IN all_windows_percent_correlation_overlap_chromHMM_percents.txt.
# THIS IS NOT DONE HERE AS INSTEAD MATHEMATICA WAS USED TO OBTAIN THE ENRICHMENTS.



### PMD/non-PMD ENRICHMENTS ###

# For PMD/non-PMD extract the overlap with all windows, the lowest correlation windows and the highest correlation windows

 bedtools intersect -wao -a $1 -b $7 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/window_proportion_overlap_pmd.txt.gz;
 bedtools intersect -wao -a $1 -b $8 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/window_proportion_overlap_nonpmd.txt.gz;
 
 bedtools intersect -wao -a $2 -b $7 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/bottom_10_percent_correlation_proportion_overlap_pmd.txt.gz;
 bedtools intersect -wao -a $2 -b $8 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/bottom_10_percent_correlation_proportion_overlap_nonpmd.txt.gz;

 bedtools intersect -wao -a $3 -b $7 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/top_10_percent_correlation_proportion_overlap_pmd.txt.gz;
 bedtools intersect -wao -a $3 -b $8 | awk 'OFS="\t" {print $4,$9,$9/($3-$2)}' | sort -k1,1 | bedtools groupby -i - -g 1 -c 3 -o sum | gzip > $9/specific_state_overlaps/top_10_percent_correlation_proportion_overlap_nonpmd.txt.gz;
 
 
# Extract the overlap and proportion overlap of all windows with PMDs and non-PMDs

all_pmd_overlap=$(bedtools intersect -a $7 -b $1 | awk '{Tot+=$3-$2} END {print Tot}')
all_nonpmd_overlap=$(bedtools intersect -a $8 -b $1 | awk '{Tot+=$3-$2} END {print Tot}')
all_pmd_prop_overlap=$(echo "${all_pmd_overlap}/${all_bases}" | bc -l)
all_nonpmd_prop_overlap=$(echo "${all_nonpmd_overlap}/${all_bases}" | bc -l)

# Extract the overlap and proportion overlap of the lowest correlation windows with PMDs and non-PMDs

lowest_correlation_pmd_overlap=$(bedtools intersect -a $7 -b $2 | awk '{Tot+=$3-$2} END {print Tot}') 
lowest_correlation_nonpmd_overlap=$(bedtools intersect -a $8 -b $2 | awk '{Tot+=$3-$2} END {print Tot}')
lowest_correlation_pmd_prop_overlap=$(echo "${lowest_correlation_pmd_overlap}/${lowest_correlation_bases}" | bc -l)
lowest_correlation_nonpmd_prop_overlap=$(echo "${lowest_correlation_nonpmd_overlap}/${lowest_correlation_bases}" | bc -l)

# Extract the overlap and proportion overlap of the highest correlation windows with PMDs and non-PMDs

highest_correlation_pmd_overlap=$(bedtools intersect -a $7 -b $3 | awk '{Tot+=$3-$2} END {print Tot}') 
highest_correlation_nonpmd_overlap=$(bedtools intersect -a $8 -b $3 | awk '{Tot+=$3-$2} END {print Tot}')
highest_correlation_pmd_prop_overlap=$(echo "${highest_correlation_pmd_overlap}/${highest_correlation_bases}" | bc -l)
highest_correlation_nonpmd_prop_overlap=$(echo "${highest_correlation_nonpmd_overlap}/${highest_correlation_bases}" | bc -l)
 
# Calculate enrichments and save to file

echo "PMD enrichment: $(echo "${lowest_correlation_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $9/low_correlation_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${lowest_correlation_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $9/low_correlation_pmd_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${highest_correlation_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $9/high_correlation_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${highest_correlation_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $9/high_correlation_pmd_nonpmd_enrichments.txt


 
