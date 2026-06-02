#!/usr/bin/bash

# Matt Arnold 2023

# This script indexes a host reference genome for host read mapping using Bowtie2

# Input:    -g --genome    host reference genome fasta file (or multiple concatenated genomes)
#           -o --out_file  output file name

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -g|--genome) genome="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
    esac
    shift
done


bowtie2-build $genome $out_file # index host reference genome for host read mapping using Bowtie2