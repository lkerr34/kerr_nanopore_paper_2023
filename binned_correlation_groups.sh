#!/bin/bash

# Arguments:
# $1: path to bed file containing all windows analysed
# $2: path to tsv file containing window name and the binned correlation group that the window belongs to
# $3: path to folder where output will be stored
# $4: gene bed file (gzip compressed)
# $5: CGI bed file
# $6: chromHMM bed file
# $7: PMD bed file
# $8: non-PMD bed file

# Load required modules

module load igmm/apps/BEDTools/2.30.0

# Create a file containing the window location, name and the binned correlation group it belongs to

join -1 4 -2 1 <(sort -k4,4 $1) <(sed '1d' $2 | sort -k1,1) | awk 'OFS="\t" {print $2,$3,$4,$1,$5}' > $3/window_bin_info.bed

# Separate these into 5 subfiles

grep -w 'very_low' $3/window_bin_info.bed > $3/very_low_window_bin_info.bed
grep -w 'low' $3/window_bin_info.bed > $3/low_window_bin_info.bed
grep -w 'moderate' $3/window_bin_info.bed > $3/moderate_window_bin_info.bed
grep -w 'high' $3/window_bin_info.bed > $3/high_window_bin_info.bed
grep -w 'very_high' $3/window_bin_info.bed > $3/very_high_window_bin_info.bed

# Calculate the number of bases in all windows as well as in each of the five binned correlation groups
all_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $1)
very_low_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3/very_low_window_bin_info.bed)
low_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3/low_window_bin_info.bed)
moderate_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3/moderate_window_bin_info.bed)
high_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3/high_window_bin_info.bed)
very_high_bases=$(awk '{Tot+=$3-$2} END {print Tot}' $3/very_high_window_bin_info.bed)




### GENES AND CGIS ###

# Calculate the overlap and proportion overlap of all windows with genes and CGIs 
all_gene_overlap=$(gzip -dc $4 | bedtools merge -i - | bedtools intersect -a $1 -b - | awk '{Tot+=$3-$2} END {print Tot}')
all_CGI_overlap=$(bedtools intersect -a $1 -b $5 | awk '{Tot+=$3-$2} END {print Tot}')
all_gene_prop_overlap=$(echo "${all_gene_overlap}/${all_bases}" | bc -l)
all_CGI_prop_overlap=$(echo "${all_CGI_overlap}/${all_bases}" | bc -l)

# Check the overlap with genes and CGIs for windows in each of the five clusters

gzip -dc $4 | bedtools merge -i - | bedtools intersect -wao -a - -b $3/window_bin_info.bed | awk 'OFS="\t" {print $7,$8,$9}' | bedtools groupby -i - -g 1,2 -c 3 -o sum > $3/gene_overlaps.tsv
bedtools intersect -wao -a $3/window_bin_info.bed -b $5 | awk 'OFS="\t" {print $4,$5,$10}' | bedtools groupby -i - -g 1,2 -c 3 -o sum > $3/CGI_overlaps.tsv

# For each binned group, calculate the number of bases overlapping genes

very_low_gene_overlap=$(grep -w 'very_low' $3/gene_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
low_gene_overlap=$(grep -w 'low' $3/gene_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
moderate_gene_overlap=$(grep -w 'moderate' $3/gene_overlaps.tsv | awk '{Tot+=$3} END {print Tot}') 
high_gene_overlap=$(grep -w 'high' $3/gene_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
very_high_gene_overlap=$(grep -w 'very_high' $3/gene_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')

# For each binned group, calculate the proportion of bases overlapping genes

very_low_gene_prop_overlap=$(echo "${very_low_gene_overlap}/${very_low_bases}" | bc -l)
low_gene_prop_overlap=$(echo "${low_gene_overlap}/${low_bases}" | bc -l)
moderate_gene_prop_overlap=$(echo "${moderate_gene_overlap}/${moderate_bases}" | bc -l)
high_gene_prop_overlap=$(echo "${high_gene_overlap}/${high_bases}" | bc -l)
very_high_gene_prop_overlap=$(echo "${very_high_gene_overlap}/${very_high_bases}" | bc -l)

# For each binned group, calculate the number of bases overlapping CGIs

very_low_CGI_overlap=$(grep -w 'very_low' $3/CGI_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
low_CGI_overlap=$(grep -w 'low' $3/CGI_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
moderate_CGI_overlap=$(grep -w 'moderate' $3/CGI_overlaps.tsv | awk '{Tot+=$3} END {print Tot}') 
high_CGI_overlap=$(grep -w 'high' $3/CGI_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')
very_high_CGI_overlap=$(grep -w 'very_high' $3/CGI_overlaps.tsv | awk '{Tot+=$3} END {print Tot}')

# For each binned group, calculate the proportion of bases overlapping CGIs

very_low_CGI_prop_overlap=$(echo "${very_low_CGI_overlap}/${very_low_bases}" | bc -l)
low_CGI_prop_overlap=$(echo "${low_CGI_overlap}/${low_bases}" | bc -l)
moderate_CGI_prop_overlap=$(echo "${moderate_CGI_overlap}/${moderate_bases}" | bc -l)
high_CGI_prop_overlap=$(echo "${high_CGI_overlap}/${high_bases}" | bc -l)
very_high_CGI_prop_overlap=$(echo "${very_high_CGI_overlap}/${very_high_bases}" | bc -l)


# Calculate enrichments and save to file

echo "gene enrichment: $(echo "${very_low_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $3/very_low_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${very_low_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $3/very_low_gene_CGI_enrichments.txt

echo "gene enrichment: $(echo "${low_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $3/low_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${low_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $3/low_gene_CGI_enrichments.txt

echo "gene enrichment: $(echo "${moderate_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $3/moderate_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${moderate_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $3/moderate_gene_CGI_enrichments.txt

echo "gene enrichment: $(echo "${high_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $3/high_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${high_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $3/high_gene_CGI_enrichments.txt

echo "gene enrichment: $(echo "${very_high_gene_prop_overlap}/${all_gene_prop_overlap}" | bc -l)" > $3/very_high_gene_CGI_enrichments.txt
echo "CGI enrichment: $(echo "${very_high_CGI_prop_overlap}/${all_CGI_prop_overlap}" | bc -l)" >> $3/very_high_gene_CGI_enrichments.txt



### CHROMHMM ###

# Extract the overlap between all windows and the chromHMM annotations

bedtools intersect -wao -a $6 -b $3/window_bin_info.bed | awk 'OFS="\t" {print $8,$9,$4,$10}' > $3/window_overlap_bins.tsv

# For each binned group, calculate the % of bp cooresponding to each annotation

rm $3/correlation_very_low_bin_overlap_percents.txt
touch $3/correlation_very_low_bin_overlap_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w 'very_low' $3/window_overlap_bins.tsv | grep -w $i | awk -v tot_bases="${very_low_bases}" 'OFS="\t" {Tot+=$4} END {print $3,100*Tot/tot_bases}' >> $3/correlation_very_low_bin_overlap_percents.txt;
 done

rm $3/correlation_low_bin_overlap_percents.txt
touch $3/correlation_low_bin_overlap_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w 'low' $3/window_overlap_bins.tsv | grep -w $i | awk -v tot_bases="${low_bases}" 'OFS="\t" {Tot+=$4} END {print $3,100*Tot/tot_bases}' >> $3/correlation_low_bin_overlap_percents.txt;
 done

rm $3/correlation_moderate_bin_overlap_percents.txt
touch $3/correlation_moderate_bin_overlap_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w 'moderate' $3/window_overlap_bins.tsv | grep -w $i | awk -v tot_bases="${moderate_bases}" 'OFS="\t" {Tot+=$4} END {print $3,100*Tot/tot_bases}' >> $3/correlation_moderate_bin_overlap_percents.txt;
 done
 
rm $3/correlation_high_bin_overlap_percents.txt
touch $3/correlation_high_bin_overlap_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w 'high' $3/window_overlap_bins.tsv | grep -w $i | awk -v tot_bases="${high_bases}" 'OFS="\t" {Tot+=$4} END {print $3,100*Tot/tot_bases}' >> $3/correlation_high_bin_overlap_percents.txt;
 done

rm $3/correlation_very_high_bin_overlap_percents.txt
touch $3/correlation_very_high_bin_overlap_percents.txt
for i in $(awk '{print $4}' $6 | sort -k1,1 | uniq | cat);
 do 
 grep -w 'very_high' $3/window_overlap_bins.tsv | grep -w $i | awk -v tot_bases="${very_high_bases}" 'OFS="\t" {Tot+=$4} END {print $3,100*Tot/tot_bases}' >> $3/correlation_very_high_bin_overlap_percents.txt;
 done

# NOTE THAT ENRICHMENTS CAN THEN BE OBTAINED BY DIVIDING VALUES IN the ..._bin_overlap percents.txt files BY CORRESPONDING
# VALUES IN THE all_windows_percent_correlation_overlap_chromHMM_percents.txt FILE GENERATED IN THE enrichment_analysis.sh SCRIPT.
# THIS IS NOT DONE HERE AS INSTEAD MATHEMATICA WAS USED TO OBTAIN THE ENRICHMENTS.



### PMDs/non-PMDs ###

# Extract the overlap between all windows and the PMD/nonPMD annotations

bedtools intersect -wao -a $7 -b $3/window_bin_info.bed | awk 'OFS="\t" {print $8,$9,$4,$10}' > $3/pmd_overlap_bins.tsv
bedtools intersect -wao -a $8 -b $3/window_bin_info.bed | awk 'OFS="\t" {print $8,$9,$4,$10}' > $3/nonpmd_overlap_bins.tsv

# Calculate the total overlap of all windows with PMDs and non-PMDs

all_pmd_overlap=$(awk '{Tot+=$4} END {print Tot}' $3/pmd_overlap_bins.tsv)
all_nonpmd_overlap=$(awk '{Tot+=$4} END {print Tot}' $3/nonpmd_overlap_bins.tsv)

# Calculate the proportion overlap of all windows with PMDs and non-PMDs

all_pmd_prop_overlap=$(echo "${all_pmd_overlap}/${all_bases}" | bc -l)
all_nonpmd_prop_overlap=$(echo "${all_nonpmd_overlap}/${all_bases}" | bc -l)


# Extract the length of overlaps with PMDs for each binned group	

very_low_pmd_overlap=$(grep -w 'very_low' $3/pmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
low_pmd_overlap=$(grep -w 'low' $3/pmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
moderate_pmd_overlap=$(grep -w 'moderate' $3/pmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
high_pmd_overlap=$(grep -w 'high' $3/pmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
very_high_pmd_overlap=$(grep -w 'very_high' $3/pmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 

# For each binned group, calculate the proportion of bases overlapping PMDs

very_low_pmd_prop_overlap=$(echo "${very_low_pmd_overlap}/${very_low_bases}" | bc -l)
low_pmd_prop_overlap=$(echo "${low_pmd_overlap}/${low_bases}" | bc -l)
moderate_pmd_prop_overlap=$(echo "${moderate_pmd_overlap}/${moderate_bases}" | bc -l)
high_pmd_prop_overlap=$(echo "${high_pmd_overlap}/${high_bases}" | bc -l)
very_high_pmd_prop_overlap=$(echo "${very_high_pmd_overlap}/${very_high_bases}" | bc -l)

# Extract the length of overlaps with non-PMDs for each binned group

very_low_nonpmd_overlap=$(grep -w 'very_low' $3/nonpmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
low_nonpmd_overlap=$(grep -w 'low' $3/nonpmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
moderate_nonpmd_overlap=$(grep -w 'moderate' $3/nonpmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
high_nonpmd_overlap=$(grep -w 'high' $3/nonpmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}') 
very_high_nonpmd_overlap=$(grep -w 'very_high' $3/nonpmd_overlap_bins.tsv | awk '{Tot+=$4} END {print Tot}')

# For each binned group, calculate the proportion of bases overlapping non-PMDs

very_low_nonpmd_prop_overlap=$(echo "${very_low_nonpmd_overlap}/${very_low_bases}" | bc -l)
low_nonpmd_prop_overlap=$(echo "${low_nonpmd_overlap}/${low_bases}" | bc -l)
moderate_nonpmd_prop_overlap=$(echo "${moderate_nonpmd_overlap}/${moderate_bases}" | bc -l)
high_nonpmd_prop_overlap=$(echo "${high_nonpmd_overlap}/${high_bases}" | bc -l)
very_high_nonpmd_prop_overlap=$(echo "${very_high_nonpmd_overlap}/${very_high_bases}" | bc -l)

# Calculate enrichments and save to file

echo "PMD enrichment: $(echo "${very_low_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $3/very_low_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${very_low_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $3/very_low_pmd_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${low_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $3/low_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${low_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $3/low_pmd_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${moderate_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $3/moderate_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${moderate_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $3/moderate_pmd_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${high_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $3/high_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${high_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $3/high_pmd_nonpmd_enrichments.txt

echo "PMD enrichment: $(echo "${very_high_pmd_prop_overlap}/${all_pmd_prop_overlap}" | bc -l)" > $3/very_high_pmd_nonpmd_enrichments.txt
echo "non-PMD enrichment: $(echo "${very_high_nonpmd_prop_overlap}/${all_nonpmd_prop_overlap}" | bc -l)" >> $3/very_high_pmd_nonpmd_enrichments.txt



