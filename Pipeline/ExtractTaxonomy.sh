#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to concatenate all contigs from a run and blast them against a chosen database (default = nr)
# Input:    -d --data_dir  path to directory containing raw data
#           -b --database  database to be used (default = nr)
#           -t --threads  number of cpu threads to use

# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -b|--database) database="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -s|--file_string) file_string="$2"; shift;;
    esac
    shift
done

install_path=$(dirname -- "$0")/

[[ -z $database ]] && database=~/4tHardDrive/db/ncbi_taxonomy # Default database is nr
[[ -z $threads ]] && threads=10 # Default number of threads is 10

for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    hit_contigs=$(cat ${data_dir}/${accession}/${accession}.diamond.txt | cut -f1  | uniq)
    for contig in $hit_contigs; do
        top_hit=$(grep $contig ${data_dir}/${accession}/${accession}.diamond.txt | head -n1 | cut -f2)
        echo $top_hit
    done
done