#!/usr/bin/bash
# Matt Arnold 2023
# This script is used to run the whole pipeline
# Input:    -d --data_dir  path to directory containing raw data
#           -g --genome    path to host reference genome (or multiple concatenated genomes)
#                          (GENOME MUST BE INDEXED WITH BOWTIE2, see Pipeline/IndexHostGenome.sh)
#           -r --raw_reads path to trimmed raw reads
#           -o --out_file  output file name
#           -m --metadata  path to metadata file containing species names for each accession number
#           -s --species_column column number of species in metadata file
#           -h --host_ref_genome_dir path to directory containing host reference genomes
# Parse command line arguments
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -g|--genome) genome="$2"; shift;;
      -r|--raw_reads) raw_reads="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
      -m|--metadata) metadata="$2"; shift;;
      -s|--species_column) species_column="$2"; shift;;
      -h|--host_ref_genome_dir) host_ref_genome_dir="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -n|--no_host_filtering) no_host_filtering=true; shift;;
    esac
    shift
done
install_path=$(dirname -- "$0")/
[[ -z $threads ]] && threads=10 # Default number of threads is 10
for directory in $(ls -d $data_dir/*/); do
  ${install_path}ReadTrimming.sh -r $directory -o $directory -t $threads
  [[ -z $genome ]] && [[ ! -z $metadata ]] && ${install_path}HostFiltering.sh -r $directory  -h $host_ref_genome_dir -o $directory -t $threads -m $metadata -s $species_column
  [[ -z $metadata ]]  && [[ ! -z $genome ]] && ${install_path}HostFiltering.sh -r $directory  -h $host_ref_genome_dir -o $directory -t $threads -g $genome
  [[ -z $metadata ]] && [[ -z $genome ]] && echo No host genome or metadata file specified, skipping host filtering
  ${install_path}Assembly.sh -r $directory -o $directory -t $threads --no_host_filtering $no_host_filtering
done

