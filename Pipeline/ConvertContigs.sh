#!/usr/bin/bash

# Matt Arnold 2024

# This script is used to convert contigs in a fasta file from SPAdes to megahit format for downstream processing

# Input:    -c --contigs  path to contigs file

# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -c|--contigs) contigs="$2"; shift;;
    esac
    shift
done

[[ -z $contigs ]] && echo "No contigs file provided" && exit 1


contig_ids=($(grep ">" $contigs))

for id in ${contig_ids[@]}; do
    intermediate_id=$(echo $id | sed 's/_length_/ len=/g' | sed 's/_cov_/ multi=/g')
    name=$(echo $intermediate_id | cut -d" " -f1)
    len=$(echo $intermediate_id | cut -d" " -f2)
    multi=$(echo $intermediate_id | cut -d" " -f3)
    new_id=$(echo $name flag=1 $multi $len)
    #echo $new_id
    sed -i "s/$id/$new_id/g" $contigs
done