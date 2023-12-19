# Sign into Eddie

ssh eddie.ecdf.ed.ac.uk


####################################################################################################



# Create an anaconda environment that contains the packages required for basecalling,
# alignment and methylation calling

make_anaconda_script= # path to make_anaconda_env.sh script
anaconda_dir= # path to main anaconda directory
conda_env_name= # name of anaconda environment to make

${make_anaconda_script} ${anaconda_dir} ${conda_env_name}



####################################################################################################



# DOWNLOAD AND PREPROCESS DATA

# Start an interactive session

qlogin -l h_vmem=16G

# Define the main directory to work in and the required script

main_dir= # path to main working directory
split_data_script= # path to split_data_into_batches.sh script

# Create directories to store data

mkdir ${main_dir}/fast5_files

# Load module to download data.

module load igmm/apps/awscli/2.1.6

# Move to directory and download the data.

cd ${main_dir}/fast5_files
aws s3 sync --no-sign-request s3://ont-open-data/gm24385_2020.09/flowcells/20200914_1354_6B_PAF27096_e7c9eae6/fast5_pass set_1
aws s3 sync --no-sign-request s3://ont-open-data/gm24385_2020.09/flowcells/20200914_1357_1-E11-H11_PAF27462_d3c9678e/fast5_pass set_2

# Split files into batches of 10

${split_data_script} ${main_dir}/fast5_files/set_1 86
${split_data_script} ${main_dir}/fast5_files/set_2 85

logout

####################################################################################################



# BASECALLING, ALIGNMENT AND METHYLATION CALLING

# Define directories and variables

bac_script= # path to basecall_align_call_methylation.sh script
main_dir= # path to main working directory
raw_data= # path to directory containing subfolders of fast5 files
anaconda_dir= # path to main anaconda directory
conda_env_name= # name of anaconda envirnoment with required packages
guppy_config=dna_r9.4.1_450bps_hac_prom.cfg # guppy configuration required for basecalling
reference_genome= # path to hg38 reference genome
mtsv_2_bed_script= # path to mtsv2bedGraph.py
LLR_threshold=2.0 # log likelihood ratio threshold used to call methylation

# Make directories to store outputs

mkdir ${main_dir}/guppy_basecalls
mkdir ${main_dir}/guppy_basecalls/set_1
mkdir ${main_dir}/guppy_basecalls/set_2

mkdir ${main_dir}/alignments
mkdir ${main_dir}/alignments/set_1
mkdir ${main_dir}/alignments/set_2

mkdir ${main_dir}/meth_calls
mkdir ${main_dir}/meth_calls/set_1
mkdir ${main_dir}/meth_calls/set_2

# Create lists of file set names

ls ${raw_data}/set_1 > ${main_dir}/file_set_names.txt
ls ${raw_data}/set_2 > ${main_dir}/file_set_names_2.txt

# Submit array jobs to bascall, align and call methylation

qsub -t 1-87 ${bac_script} ${main_dir}/file_set_names.txt ${raw_data}/set_1 ${reference_genome} ${conda_env_name} ${anaconda_dir} ${guppy_config} ${main_dir}/guppy_basecalls/set_1 ${main_dir}/alignments/set_1 ${main_dir}/meth_calls/set_1 ${mtsv_2_bed_script}
qsub -t 1-86 ${bac_script} ${main_dir}/file_set_names_2.txt ${raw_data}/set_2 ${reference_genome} ${conda_env_name} ${anaconda_dir} ${guppy_config} ${main_dir}/guppy_basecalls/set_2 ${main_dir}/alignments/set_2 ${main_dir}/meth_calls/set_2 ${mtsv_2_bed_script}

# Merge Nanopolish tsv outputs to one file 

awk 'FNR>1 || NR==1' ${main_dir}/meth_calls/set_1/meth_calls_file_set_*.tsv ${main_dir}/meth_calls/set_2/meth_calls_file_set_*.tsv > ${main_dir}/meth_calls/all_meth_calls.tsv


####################################################################################################



# BULK METHYLATION

# MAIN OUTPUTS: 
# tsv and bed files giving bulk methylation of CpGs within canonical chromosomes (inc. X,Y,M)
# bed file giving bulk methylation of CpGs within autosomes that have minimum coverage 10 (poorly sequenced regions excluded)

qlogin -l h_vmem=16G

# Define directories and variables

bulk_meth_script= # path to bulk_methylation.sh script
meth_tsv= # path to combined Nanopolish tsv file
bulk_meth_dir= # path to output directory for bulk meth info
anaconda_dir= # path to main anaconda folder
conda_env_name=guppyenv # name of anaconda environment with required packages
meth_freq_script= # path to calculate_methylation_frequency_include_ambiguous_calls_and_ignore_groupings.py script
excluded_regions= # path to bed file containing any regions to exclude from analysis
LLR_threshold=2.0 # log likelihood ratio threshold used to call methylation

# Submit job to extract bulk methylation levels 

qsub ${bulk_meth_script} ${meth_tsv} ${bulk_meth_dir} ${anaconda_dir} ${conda_env_name} ${meth_freq_script} ${excluded_regions} ${LLR_threshold}

logout



##############################################################################################################



# BULK ANALYSIS: GENOMIC WINDOWS

# MAIN OUTPUTS:
# tsv file containing bulk-level statistics in genomic windows

# Define directories and variables

window_bulk_stats_script= # path to window_population_statistics.sh script
chrom_sizes= # path to file containing chromosome sizes for hg38 reference genome
window_size=100 # genomic window size (in units of kb)
window_dir= # path to directory where output from window analysis shoud be stored
CpG_meth_summary= # path to bed file containing bulk autosomal CpG methylation levels (as outputted in bulk_methylation.sh script)
R_window_bulk_stats_script= # path to window_population_statistics.R R script

# Submit job to calculate population-level statistics for each genomic window

qsub ${window_bulk_stats_script} ${chrom_sizes} ${window_size} ${window_dir} ${CpG_meth_summary} ${R_window_bulk_stats_script}

# Check the number of windows containing methylation info for at least 100 CpGs

gzip -dc ${window_dir}/population_level/window_population_stats.tsv.gz | wc -l



####################################################################################################



# SINGLE-READ STATISTICS

# MAIN OUTPUTS:
# bed files containing read-level statistics
# bed files containing read-level statistics along with window info for all reads that align entirely within one window
# a single bed file containing the merged read-level statistics along with window info for all reads that align entirely within one window
# a "shuffled" version of reads which assigns each read entirely contained within one window to a random window

qlogin -l h_vmem=128G

# Load bedtools
module load igmm/apps/BEDTools/2.30.0 

# Define directories and variables

single_read_stats_script= # path to single_read_statistics.sh script
main_dir= # path to main working directory
file_names_1= # path to first text file with list of folder names where Nanopore fast5 files are stored
file_names_2= # path to second text file with list of folder names where Nanopore fast5 files are stored
R_single_read_stats_script= # path to single_read_statistics.R R script
meth_bed_files_1= # path to first directory containing methylation bed files
meth_bed_files_2= # path to second directory containing methylation bed files
window_dir= # path to window analysis directory
window_bed_file= # path to bed file containing genomic windows
LLR_threshold=2.0 # log likelihood ratio threshold used to call methylation

# Create folders to store stats output

mkdir ${main_dir}/read_stats
mkdir ${main_dir}/read_stats/set_1
mkdir ${main_dir}/read_stats/set_2

# Create folders to store window stats output

mkdir ${window_dir}/single_molecule_level
mkdir ${window_dir}/single_molecule_level/set_1
mkdir ${window_dir}/single_molecule_level/set_2

# Submit array jobs to calculate per-read stats.

qsub -t 1-87 ${single_read_stats_script} ${file_names_1} ${R_single_read_stats_script} ${meth_bed_files_1} ${main_dir}/read_stats/set_1 ${window_bed_file} ${window_dir}/single_molecule_level/set_1 ${LLR_threshold}
qsub -t 1-86 ${single_read_stats_script} ${file_names_2} ${R_single_read_stats_script} ${meth_bed_files_2} ${main_dir}/read_stats/set_2 ${window_bed_file} ${window_dir}/single_molecule_level/set_2 ${LLR_threshold}

# Extract stats information for canonical autosomal reads containing more than 100 CpGs
# Extract the genomic positions of these reads
# Extract the number of CpGs and mean methylation in these reads

cat ${main_dir}/read_stats/set_1/stats_file_set_*.bed.gz ${main_dir}/read_stats/set_2/stats_file_set_*.bed.gz | gzip -dc | grep -v -P 'chr.+_' | grep -v -E 'chrX|chrY|chrM|chromosome' | awk 'OFS="\t" {if($8>=100) {print $4,$9,$10,$11,$12,$13,$14}}' | gzip > ${main_dir}/read_stats/read_stats_min_100_CpGs.tsv.gz
cat ${main_dir}/read_stats/set_1/stats_file_set_*.bed.gz ${main_dir}/read_stats/set_2/stats_file_set_*.bed.gz | gzip -dc | grep -v -P 'chr.+_' | grep -v -E 'chrX|chrY|chrM|chromosome' | awk 'OFS="\t" {if($8>=100) {print $1,$2,$3,$4}}' | gzip > ${main_dir}/read_stats/min_100_CpG_reads_genomic_coords.bed.gz 
cat ${main_dir}/read_stats/set_1/stats_file_set_*.bed.gz ${main_dir}/read_stats/set_2/stats_file_set_*.bed.gz | gzip -dc | grep -v -P 'chr.+_' | grep -v -E 'chrX|chrY|chrM|chromosome' | awk 'OFS="\t" {if($8>=100) {print $4,$8,$9}}' | gzip >  ${main_dir}/read_stats/reads_with_over_100_CpGs_num_CpGs_means.tsv.gz


# Join all the window output files containing the stats and window info together

cat ${window_dir}/single_molecule_level/set_1/stats_file_set* ${window_dir}/single_molecule_level/set_2/stats_file_set* > ${window_dir}/single_molecule_level/stats_file_set.bed.gz

# Extract the read name and correlation associated with reads containing over 100 CpGs

gzip -dc ${window_dir}/single_molecule_level/stats_file_set.bed.gz | awk 'OFS="\t" {if($5>=100){print $4,$9}}' | gzip > ${window_dir}/single_molecule_level/reads_correlation.tsv.gz

# Shuffle the reads contained entirely within a single window that have more than 100 CpGs and randomly assign them to a window
# Then calculate the mean correlation for each window using the randomly assigned reads
# (This provides control data as to what would be the expected window correlation distribution would be if there was no associateion between a reads genomic location and it's correlation)

gzip -dc ${window_dir}/single_molecule_level/stats_file_set.bed.gz | awk '{if($5>=100){print $9}}' | shuf > ${window_dir}/single_molecule_level/shuffled_correlation.txt
gzip -dc ${window_dir}/single_molecule_level/stats_file_set.bed.gz | awk '{if($5>=100){print $13}}' | paste -d "\t" - ${window_dir}/single_molecule_level/shuffled_correlation.txt | grep -w -v 'NA' | sort -k1,1 > ${window_dir}/single_molecule_level/shuffled_correlation_window.tsv
bedtools groupby -i ${window_dir}/single_molecule_level/shuffled_correlation_window.tsv -g 1 -c 2 -o mean  > ${window_dir}/single_molecule_level/mean_shuffled_correlation.tsv

# Shuffle the reads contained entirely within a single window that have more than 100 CpGs and randomly assign them to a window
# Then calculate the mean RTS for each window using the randomly assigned reads
# (This provides control data as to what would be the expected window RTS distribution would be if there was no associateion between a reads genomic location and it's RTS)

gzip -dc ${window_dir}/single_molecule_level/stats_file_set.bed.gz | awk '{if($5>=100){print $10}}' | shuf > ${window_dir}/single_molecule_level/shuffled_RTS.txt
gzip -dc ${window_dir}/single_molecule_level/stats_file_set.bed.gz | awk '{if($5>=100){print $13}}' | paste -d "\t" - ${window_dir}/single_molecule_level/shuffled_RTS.txt | grep -w -v 'NA' | sort -k1,1 > ${window_dir}/single_molecule_level/shuffled_RTS_window.tsv
bedtools groupby -i ${window_dir}/single_molecule_level/shuffled_RTS_window.tsv -g 1 -c 2 -o mean  > ${window_dir}/single_molecule_level/mean_shuffled_RTS.tsv

logout


####################################################################################################



# MOST/LEAST HETEROGENEOUS WINDOWS (DEFINED VIA CORRELATION)

# MAIN OUTPUTS:
# tsv file containing mean single-read statistics for each window
# bedGraphs with the mean and correlation associated with each window
# tsv files with the stats associated with the 10% of windows with the highest correlation and the 10% of windows with the lowest correlation
# bed files with the genomic locations of the highest and lowest correlation windows
# bed files containing bulk methylation levels for CpGs in the highest and lowest correlation windows
# bed files containing bulk methylation levels for all CpGs excluding those in the highest correlation windows (similarly for lowest correlation windows)

qlogin -l h_vmem=128G

# Define directories and variables

high_low_correlation_script= # path to highest_lowest_correlation_windows.sh script
mean_window_stats_script= # path to calculate_mean_window_stats.R R script
window_dir= # path to window analysis directory
sm_window_dir= # path to single-molecule window analysis directory
window_bed_file= # path to bed file containing genomic windows
stats_window_file= # path to gzipped single-read stats bed file
CpG_meth_summary= # path to bed file containing bulk CpG methylation levels (as outputted in bulk_methylation.sh script)

# Submit script to calculate mean single-read statistics for each window and identify the windows associated with the highest/lowest correlation, 

qsub ${high_low_correlation_script} ${mean_window_stats_script} ${stats_window_file} ${sm_window_dir} ${window_bed_file} ${CpG_meth_summary}

# Combine the population statistics and single-molecule statistics into one file

join -1 7 -2 1 <(gzip -dc ${window_dir}/population_level/window_population_stats.tsv.gz | sort -k7,7) <(gzip -dc ${sm_window_dir}/mean_window_stats.tsv.gz | sort -k1,1) | awk 'OFS="\t" {print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13}'  | gzip > ${window_dir}/population_vs_single_molecule_stats.tsv.gz

logout



####################################################################################################



# IDENTIFY PMDS AND NON-PMDS

# MAIN OUTPUTS
# bed files containing the genomic locations of PMDs and non-PMDs (each with min length 200kb)

# Define directories and variables

pmd_nonpmd_identification_script= # path to identify_PMDs_nonPMDs.sh script
main_dir= # path to main working directory
CpG_bulk_meth_info= # path to tsv bulk methylation summary file 
chrom_sizes= # path to sorted chromosome size file (sort -k1,1)
excluded_regions= # path to bed file containing any regions to exclude from analysis

qsub ${pmd_nonpmd_identification_script} ${main_dir} ${CpG_bulk_meth_info} ${chrom_sizes} ${excluded_regions}



####################################################################################################



# CLEAN CHROMHMM ANNOTATION FILES

# MAIN OUTPUTS:
# a "clean" chromHMM file where all ambiguous regions (i.e. regions that are associated with 
# multiple chromatin states due to lift over from hg19 to hg38) are removed

# hg19 chromHMM annotations can be downloaded from UCSC and converted into bed format via
# gzip -dc ${UCSC_FILE} | sed '1d' | awk 'OFS="\t" {print $2,$3,$4,$5}' > ${OUTPUT_BED}
# Liftover to the hg38 genome and sort file
# Remove non-canonical chromosomes and X, Y, M chromosomes

# Define directories and variables

clean_chromHMM_script= # path to clean_chromHMM_file.sh script
original_chromHMM_file= # path to hg38 chromHMM bed file (generated as described above)
output_dir= # path to directory where all output will be stored

${clean_chromHMM_script} ${original_chromHMM_file} ${output_dir}
${clean_chromHMM_script} ${original_chromHMM_file} ${output_dir}



####################################################################################################



# ENRICHMENT ANALYSIS

# MAIN OUTPUTS:
# text files containing the enrichment of the lowest and highest correlation windows in genes and CGIs
# text files containing the percentage overlap of all windows, the lowest correlation windows and the highest correlation windows in chromHMM states 
# (can subsequently be used to obtain enrichments)
# text files containing the enrichment of the lowest and highest correlation windows in PMDs and non-PMDs

# Can make bed file from gencode gtf file
# gzip -dc ${INPUT_GTF.gz}  | grep -P "\tgene\t" - | cut -f1,4,5,7,9 | sed 's/[[:space:]]/\t/g' | sed 's/[;|"]//g' | awk -F $'\t' 'BEGIN { OFS=FS } { print $1,$2-1,$3,$6,".",$4,$10,$12,$14 }' | sort -k1,1 -k2,2n | gzip > ${OUTPUT}.bed.gz 

qlogin -l h_vmem=32G

# Define directories and variables

enrichment_script= # path to enrichment_analysis.sh script
all_windows= # path to bed file containing all genomic windows analysed
low_correlation_windows= # path to bed file containing lowest correlation windows
high_correlation_windows= # path to bed file containing highest correlation windows
gene_file= # path to gzipped bed file of gene locations
CGI_file= # path to bed file of CGI locations
chromHMM_file= # path to bed file of "clean" chromHMM annotations
PMD_file= # path to PMD location bed file
nonPMD_file= # path to non-PMD location bed file
output_dir= # path to single-molecule window analysis directory

# Run script to perform enrichment analysis

${enrichment_script} ${all_windows} ${low_correlation_windows} ${high_correlation_windows} ${gene_file} ${CGI_file} ${chromHMM_file} ${PMD_file} ${nonPMD_file} ${output_dir}

logout



####################################################################################################



# ENRICHMENT OF BINNED CORRELATION GROUPS IN HETEROCHROMATIN AND PMDs

# MAIN OUTPUTS:
# text files containing enrichment of each of the five binned correlation groups in genes and CGIs
# text files containing the percentage overlap of all windows and each of the five binned correlation groups with chromHMM states 
# (can subsequently be used to obtain enrichments)
# text files containing the enrichment of the lowest and highest correlation windows in PMDs and non-PMDs

qlogin -l h_vmem=16G

# Define directories and variables

binned_correlation_script= # path to binned_correlation_groups.sh script
all_windows= # path to bed file containing all genomic windows analysed
window_bin_info= # path to tsv file containing window and binned correlation group info (here R was used to create groups). It is assumed that five groups with names "very_low", "low", "moderate", "high" and "very_high" were created
output_dir= # path to directory to store output from binned correlation analysis 
gene_file= # path to gzipped bed file of gene locations
CGI_file= # path to bed file of CGI locations
chromHMM_file= # path to bed file of "clean" chromHMM annotations
PMD_file= # path to bed file of PMD locations
nonPMD_file= # path to bed file of non-PMD locations

${binned_correlation_script} ${all_windows} ${window_bin_info} ${output_dir} ${gene_file} ${CGI_file} ${chromHMM_file} ${PMD_file} ${nonPMD_file}

logout



####################################################################################################



# CORRELATION BETWEEN CORRELATION AND OVERLAP WITH HETEROCHROMATIN/PMDS (READ AND WINDOW LEVEL)

# MAIN OUTPUTS:
# tsv file containing the correlation and overlap with heterochromatin/PMDs for each individual read
# tsv file containing the mean correlation and 

qlogin -l h_vmem=16G

correlation_overlap_with_heterochromatin_PMDs_script= # path to correlation_and_overlap_with_heterochromatin_PMDs.sh script
chromHMM_file= # path to bed file of "clean" chromHMM annotations
PMD_file= # path to PMD location bed file
read_output_dir= # path to output directory for read correlation info
window_output_dir= # path to output folder for window correlation info
min_100CpG_reads= # path to gzipped bed file containing genomic locations of reads with minimum 100 CpGs
min_100CpG_reads_stats= # path to gzipped file containing stats for autosomal reads with minimum 100 CpGs
window_bin_info= # path to bed file containing binned correlation info for all windows
mean_window_stats= # path to gzipped file containing mean stats for genomic windows

${correlation_overlap_with_heterochromatin_PMDs_script} ${chromHMM_file} ${PMD_file} ${read_output_dir} ${window_output_dir} ${min_100CpG_reads} ${min_100CpG_reads_stats} ${window_bin_info} ${mean_window_stats}

logout



####################################################################################################



# ANALYSIS OF HIGH/LOW CORRELATION READS OBTAINED VIA MIXTURE MODELLING

# MAIN OUTPUTS:
# tsv files containing stats for low and high correlation read groups (gzipped)
# bed files containing genomic location of low and high correlation reads (gzipped)
# text files containing percent overlap of all reads, low correlation reads and high correlation reads with chromHMM annotations
# (can be used to calculate enrichments of low correlation reads and high correlation reads in each chromHMM annotation)
# text files containing enrichment of low correlation reads and high correlation reads in PMDs and non-PMDs

# It is assumed that text files containing the read names associated with the low and high correlation groups have been created
# In our case a mixed model approach was used locally in Mathematica to identify reads associated with each group and output results

qlogin -l h_vmem=16G

# EXRTACT THE OVERLAP OF ALL READS WITH AT LEAST 100 CPGS WITH DIFFERENT STATES

# Define directories and variables

all_reads_overlap_script= # path to state_overlaps_for_all_reads.sh script
min_100CpG_reads= # path to gzipped bed file containing genomic locations of reads with minimum 100 CpGs
output= # path to output directory for read correlation info
chromHMM_file= # path to bed file of "clean" chromHMM annotations
PMD_file= # path to PMD locations bed file
nonPMD_file= # path to non-PMD locations bed file

${all_reads_overlap_script} ${min_100CpG_reads} ${output} ${chromHMM_file} ${PMD_file} ${nonPMD_file}


# EXRACT STATS AND OVERLAPS OF READS WITH DIFFERENT STATES FOR THE LOW AND HIGH CORRELATION READS

# Define directories and variables

group_stats_and_overlaps_script= # path to basic_stats_and_overlaps_for_read_group.sh script
low_correlation_reads= # path to text file containing names of low correlation reads (see above)
high_correlation_reads= # path to text file containing names of high correlation reads (see above)
min_100CpG_reads_stats= # path to gzipped file containing stats for autosomal reads with minimum 100 CpGs
low_correlation_output_folder= # path to folder where output for low correlation reads will be stored
high_correlation_output_folder= # path to folder where output for high correlation reads will be stored
min_100CpG_reads= # path to gzipped bed file containing genomic locations of reads with minimum 100 CpGs
chromHMM_file= # path to bed file of "clean" chromHMM annotations
PMD_file= # path to PMD location bed file
nonPMD_file= # path to non-PMD location bed files

${group_stats_and_overlaps_script} ${low_correlation_reads} ${min_100CpG_reads_stats} ${low_correlation_output_folder} ${min_100CpG_reads} ${chromHMM_file} ${PMD_file} ${nonPMD_file}
${group_stats_and_overlaps_script} ${high_correlation_reads} ${min_100CpG_reads_stats} ${high_correlation_output_folder} ${min_100CpG_reads} ${chromHMM_file} ${PMD_file} ${nonPMD_file}

# Enrichments were then calculated via mathematica

logout



####################################################################################################



# EXTRACT INFO REQUIRED TO CALCULATE DISTANCE-DEPENDENT CORRELATIONS FROM ALL READS

# MAIN OUTPUT:
# bed files containing the information required to calculate distance-dependent correlations for each read (gzip compressed)
# bed files containing the information required to calculate distance-dependent correlations for reads with minimum 100 CpGs (gzip compressed)
# single file containing all information required to calucalte distance-dependent correlations for reads with minimum 100 CpGs (gzip compressed)

qlogin -l h_vmem=16G

# Define directories and variables

split_RTS_dist_dep_info_script= # path to read_distance_info_min_100_CpGs.sh script 
file_names_1= # path to first text file with list of folder names where Nanopore fast5 files are stored
file_names_2= # path to second text file with list of folder names where Nanopore fast5 files are stored
meth_bed_files_1= # path to directory containing first set of methylation bed files
meth_bed_files_2= # path to directory containing second set of methylation bed files
output_dir_1= # path to directory to store outputs for all reads in first set
output_dir_2= # path to directory to store outputs for all reads in second set
min_100CpG_reads= # path to gzipped bed file containing genomic locations of reads with minimum 100 CpGs
output_dir_min_100CpGs_1= # path to directory to store outputs for reads with minimum 100 CpGs from first set
output_dir_min_100CpGs_2= # path to directory to store outputs for reads with minimum 100 CpGs from second set

# Submit array jobs to extract info required to compute distance-dependent correlations from all reads and reads with minimum 100 CpGs

qsub -t 1-87 ${split_RTS_dist_dep_info_script} ${file_names_1} ${meth_bed_files_1} ${output_dir_1} ${min_100CpG_reads} ${output_dir_min_100CpGs_1}
qsub -t 1-86 ${split_RTS_dist_dep_info_script} ${file_names_2} ${meth_bed_files_2} ${output_dir_2} ${min_100CpG_reads} ${output_dir_min_100CpGs_2}

logout



###################################################################################################



# CALCULATE DISTANCE-DEPENDENT CORRELATIONS (all CpGs vs. all CpGs)

# MAIN OUTPUTS:
# tsv files containing the mean distance-dependent correlations for reads in the low correlation group and for reads in the high correlation group

qlogin -l h_vmem=128G

# Define directories and variables

dist_dep_correlation_script= # path to all_vs_all_batch_individual_molecule_distance_dependent_correlations.sh script
file_names_1= # path to first text file with list of folder names where Nanopore fast5 files are stored
file_names_2= # path to second text file with list of folder names where Nanopore fast5 files are stored
dist_info_min_100CpG_reads_1= # path to directory containing gzipped bed file with info required to compute distance-dependent correlations for reads with minimum 100 CpGs from first set
dist_info_min_100CpG_reads_2= # path to directory containing gzipped bed file with info required to compute distance-dependent correlations for reads with minimum 100 CpGs from second set
low_correlation_reads= # path to gzipped bed file of low correlation reads
high_correlation_reads= # path to gzipeed bed file of high correlation reads
low_corr_dist_output_dir_1= # path to output directory for distance info associated with low correlation reads in the first set
low_corr_dist_output_dir_2= # path to output directory for distance info associated with low correlation reads in the second set
high_corr_dist_output_dir_1= # path to output directory for distance info associated with high correlation reads in the first set
high_corr_dist_output_dir_2= # path to output directory for distance info associated with high correlation reads in the second set
R_dist_dep_correlation_script= # path to all_vs_all_individual_molecule_distance_dependent_correlations.R R script
LLR_threshold=2.0 # log likelihood ratio threshold used to call methylation
low_corr_max_dist=2000 # max distance considered 
high_corr_max_dist=600 # max distance considered
low_corr_corr_output_dir_1= # path to output directory for correlations associated with low correlation reads in the first set
low_corr_corr_output_dir_2= # path to output directory for correlations associated with low correlation reads in the second set
high_corr_corr_output_dir_1= # path to output directory for correlations associated with high correlation reads in the first set
high_corr_corr_output_dir_2= # path to output directory for correlations associated with high correlation reads in the second set


# Submit array jobs to calculate individual distance-dependent correlations

qsub -t 1-87 ${dist_dep_correlation_script} ${file_names_1} ${dist_info_min_100CpG_reads_1} ${low_correlation_reads} ${low_corr_dist_output_dir_1} ${R_dist_dep_correlation_script} ${LLR_threshold} ${low_corr_max_dist} ${low_corr_corr_output_dir_1}
qsub -t 1-86 ${dist_dep_correlation_script} ${file_names_2} ${dist_info_min_100CpG_reads_2} ${low_correlation_reads} ${low_corr_dist_output_dir_2} ${R_dist_dep_correlation_script} ${LLR_threshold} ${low_corr_max_dist} ${low_corr_corr_output_dir_2}

qsub -t 1-87 ${dist_dep_correlation_script} ${file_names_1} ${dist_info_min_100CpG_reads_1} ${high_correlation_reads} ${high_corr_dist_output_dir_1} ${R_dist_dep_correlation_script} ${LLR_threshold} ${high_corr_max_dist} ${high_corr_corr_output_dir_1}
qsub -t 1-86 ${dist_dep_correlation_script} ${file_names_2} ${dist_info_min_100CpG_reads_2} ${high_correlation_reads} ${high_corr_dist_output_dir_2} ${R_dist_dep_correlation_script} ${LLR_threshold} ${high_corr_max_dist} ${high_corr_corr_output_dir_2}


# Calculate the mean distance-dependent correlations over all reads
# Since there are so many files we first calculate the mean over each batch and then take the mean of the batches

# Load bedtools
module load igmm/apps/BEDTools/2.30.0 

low_corr_dir= # path to directory to store overall low correlation read outputs
high_corr_dir= # path to directory to store overall high correlation read outputs

# Calculate the mean distance-dependent correlations in each low correlation batch

for i in {0..86};\
 do echo $i; awk 'FNR>1' ${low_corr_corr_output_dir_1}/file_set_$i/*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${low_corr_dir}/low_corr_mean_correlations_set_1_$i.tsv;
done

for i in {0..85};\
 do echo $i; awk 'FNR>1' ${low_corr_corr_output_dir_2}/file_set_$i/*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${low_corr_dir}/low_corr_mean_correlations_set_2_$i.tsv;
done

# Calculate the mean distance-dependent correlations in each high correlation batch

for i in {0..86};\
 do echo $i; awk 'FNR>1' ${high_corr_corr_output_dir_1}/file_set_$i/*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${high_corr_dir}/high_corr_mean_correlations_set_1_$i.tsv;
done

for i in {0..85};\
 do echo $i; awk 'FNR>1' ${high_corr_corr_output_dir_2}/file_set_$i/*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${high_corr_dir}/high_corr_mean_correlations_set_2_$i.tsv;
done

# Calculate overall mean distance-dependent correlations

cat ${low_corr_dir}/low_corr_mean_correlations_set_*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${low_corr_dir}/all_vs_all_low_corr_reads_mean_distance_dependent_correlations.tsv
rm ${low_corr_dir}/low_corr_mean_correlations_set_*.tsv

cat ${high_corr_dir}/high_corr_mean_correlations_set_*.tsv | sort -k1,1n | bedtools groupby -g 1 -c 2 -o mean > ${high_corr_dir}/all_vs_all_high_corr_reads_mean_distance_dependent_correlations.tsv
rm ${high_corr_dir}/high_corr_mean_correlations_set_*.tsv

logout