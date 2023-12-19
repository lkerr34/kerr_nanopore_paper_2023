#!/bin/bash

# Arguments:
# $1: path to anaconda directory
# $2: name of anaconda environment to create

mkdir $1
mkdir $1/envs
mkdir $1/pkgs

# Load required modules

module load anaconda/5.0.1

# Configure enviroments and package directories

conda config --add envs_dirs $1/envs
conda config --add pkgs_dirs $1/pkgs

conda --add channel default
conda --add channel bioconda
conda --add channel conda-forge

# Create and activate environment

conda create -n $2 python=3.8 pip=20.3.4
source activate $2

# Install software from conda 

conda install -c bioconda samtools
conda install -c bioconda bedtools
conda install -c bioconda minimap2
conda install -c bioconda nanopolish
conda install -c bioconda ont_vbz_hdf_plugin
conda install -c bioconda pysam
conda install -c anaconda pandas