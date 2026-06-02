#!/usr/bin/bash

# Matt Arnold 2024

# This script is used to convert a csv of contigs to bowtie indices, in order to map raw reads to the contigs

# Input:    -d --data_dir  path to directory containing data
#           -c --contigs  path to contigs file  
#           -m --metadata  path to metadata file containing information about the contigs
#           -o --output   output directory for bowtie indices
#           -l --length   minimum length of contigs to include in bowtie indices



# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -c|--contigs) contigs="$2"; shift;;
      -m|--metadata) metadata="$2"; shift;;
      -o|--output) output="$2"; shift;;
      -l|--length) length="$2"; shift;;
    esac
    shift
done

[[ -z $length ]] && length=250 # Default minimum length of contigs is 250


awk -v len=$length -F ',' ' { if ($8 > len) print } ' $contigs > all_virus_hits_filtered_length${length}.csv

while IFS= read -r contig; do
    contig_id=$(echo $contig | cut -d',' -f2)
    grep $contig_id $contigs -A1 | tee ./all_viral_contigs.fasta >> ${data_dir}/all_viral_contigs.fasta
done < all_virus_hits_filtered_length${length}.csv


# Create bowtie indices for the contigs
bowtie2-build ${data_dir}/all_viral_contigs.fasta ${output}/all_viral_contigs

# Run bowtie2 to map the raw reads to the contigs
