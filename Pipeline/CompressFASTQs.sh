#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to compress useful FASTQ files and remove unnecessary files

# Input:    -r --reads path to host filtered, trimmed reads
#           -t --threads number of threads to use

# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -r|--reads) reads="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
    esac
    shift
done

[[ -z $threads ]] && threads=10 # Default number of threads is 10

a=($(ls $reads))

for file in ${a[@]}; do 
    [[ ${file: -2} == "fq" ]] || [[ ${file: -5} == "fastq" ]] && pigz -vf -p $threads ${reads}/${file}
done