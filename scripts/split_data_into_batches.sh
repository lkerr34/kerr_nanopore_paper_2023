#!/bin/bash

# Arguments
# $1: path to the fast5 files
# $2: one less than the number of batches desired

for ((i=0; i<=$2; i++))
do

echo $i

mkdir $1/file_set_$i
output_dir=$1/file_set_$i

l=$((i*10))
m=$((i*10+9))

for ((j=$l; j<=$m; j++))
do

mv $1/*_$j.fast5 $output_dir

done

done
