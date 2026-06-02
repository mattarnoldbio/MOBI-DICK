#!/usr/bin/bash

# Matt Arnold 2023

# Check progress

# Input: -r --raw_reads path to raw reads

# Parse arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -r|--raw_reads) raw_reads="$2"; shift;;
    esac
    shift
done

counter=0

for directory in $raw_reads/*; do
    unmapped=$(ls $directory | grep contigs_out)
    [[ -z $unmapped ]]  || ((counter++)) 
done

n_files=$(ls -ld  ${raw_reads}/*/ | wc -l)

echo Progress: $counter / $n_files files processed