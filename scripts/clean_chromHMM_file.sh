#!/bin/bash

# Arguments:
# $1: path to hg38 chromHMM file (obtained by downloading hg19 data from UCSC, converting to bed format and lifting over to hg38 genome)
# $2: path to folder where outputs will be stored

# Start an interactive session

qlogin -l h_vmem=16G

# Load required modules

module load igmm/apps/BEDTools/2.30.0

# Define important directories

genome_dir=/exports/igmm/eddie/sproul-nanopore/GM24385_R9.4.1/genome_annotations

# Create a bed file for each type of chromHMM annotation

grep '1_Active_Promoter' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"active_promoter"}' > $2/chromHMM_active_promoter.bed
grep '2_Weak_Promoter' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"weak_promoter"}' > $2/chromHMM_weak_promoter.bed
grep '3_Poised_Promoter' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"poised_promoter"}' > $2/chromHMM_poised_promoter.bed
grep 'Strong_Enhancer' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"strong_enhancer"}' > $2/chromHMM_strong_enhancer.bed
grep 'Weak_Enhancer' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"weak_enhancer"}' > $2/chromHMM_weak_enhancer.bed
grep '8_Insulator' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"insulator"}' > $2/chromHMM_insulator.bed
grep '9_Txn_Transition' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"translational_transition"}' > $2/chromHMM_translational_transition.bed
grep '10_Txn_Elongation' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"translational_elongation"}' > $2/chromHMM_translational_elongation.bed
grep '11_Weak_Txn' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"weakly_transcribed"}' > $2/chromHMM_weakly_transcribed.bed
grep '12_Repressed' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"repressed"}' > $2/chromHMM_repressed.bed
grep '13_Heterochrom/lo' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"heterochromatin"}' > $2/chromHMM_heterochromatin.bed
grep 'Repetitive/CNV' $1 | sort -k1,1 -k2,2n | bedtools merge -i - | awk 'OFS="\t" {print $1,$2,$3,"repetitive"}' > $2/chromHMM_repetitive.bed

# Remove any regions which are ambiguously defined in more than one chromatin state

cat $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_active_promoter.bed -b - | bedtools subtract -a $2/chromHMM_active_promoter.bed -b - > $2/chromHMM_active_promoter_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_weak_promoter.bed -b - | bedtools subtract -a $2/chromHMM_weak_promoter.bed -b - > $2/chromHMM_weak_promoter_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_poised_promoter.bed -b - | bedtools subtract -a $2/chromHMM_poised_promoter.bed -b - > $2/chromHMM_poised_promoter_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_strong_enhancer.bed -b - | bedtools subtract -a $2/chromHMM_strong_enhancer.bed -b - > $2/chromHMM_strong_enhancer_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_weak_enhancer.bed -b - | bedtools subtract -a $2/chromHMM_weak_enhancer.bed -b -> $2/chromHMM_weak_enhancer_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_insulator.bed -b - | bedtools subtract -a $2/chromHMM_insulator.bed -b - > $2/chromHMM_insulator_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_translational_transition.bed -b - | bedtools subtract -a $2/chromHMM_translational_transition.bed -b - > $2/chromHMM_translational_transition_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_translational_elongation.bed -b - | bedtools subtract -a $2/chromHMM_translational_elongation.bed -b - > $2/chromHMM_translational_elongation_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_weakly_transcribed.bed -b - | bedtools subtract -a $2/chromHMM_weakly_transcribed.bed -b - > $2/chromHMM_weakly_transcribed_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_heterochromatin.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_repressed.bed -b - | bedtools subtract -a $2/chromHMM_repressed.bed -b - > $2/chromHMM_repressed_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_repetitive.bed | bedtools intersect -a $2/chromHMM_heterochromatin.bed -b - | bedtools subtract -a $2/chromHMM_heterochromatin.bed -b - > $2/chromHMM_heterochromatin_clean.bed
cat $2/chromHMM_active_promoter.bed $2/chromHMM_weak_promoter.bed $2/chromHMM_poised_promoter.bed $2/chromHMM_strong_enhancer.bed $2/chromHMM_weak_enhancer.bed $2/chromHMM_insulator.bed $2/chromHMM_translational_transition.bed $2/chromHMM_translational_elongation.bed $2/chromHMM_weakly_transcribed.bed $2/chromHMM_repressed.bed $2/chromHMM_heterochromatin.bed | bedtools intersect -a $2/chromHMM_repetitive.bed -b - | bedtools subtract -a $2/chromHMM_repetitive.bed -b - > $2/chromHMM_repetitive_clean.bed

# Combine these together into a single chromHMM file

rm $2/hg38_chromHMM_autosomes_clean.bed # Remove file if it already exists
cat $2/*clean.bed > $2/hg38_chromHMM_autosomes_clean.bed

# Remove intermediate files

rm $2/chromHMM_*

logout
