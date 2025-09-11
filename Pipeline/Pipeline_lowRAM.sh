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

# PARSE COMMAND LINE ARGUMENTS
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -g|--genome) genome="$2"; shift;;
      -m|--metadata) metadata="$2"; shift;;
      -s|--species_column) species_column="$2"; shift;;
      -h|--host_ref_genome_dir) host_ref_genome_dir="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -n|--no_host_filtering) no_host_filtering=true; shift;;
      -b|--database) database="$2"; shift;;
      -k|--krona_tools_db) krona_tools_db="$2"; shift;;
      -p|--no_prinseq) no_prinseq=true; shift;;
    esac
    shift
done

# SET DEFAULTS
install_path=$(dirname -- "$0")/ # Get path to install directory

[[ -z $threads ]] && threads=10 # Default number of threads is 10
[[ -z $database ]] && database=/db/diamond/nr.dmnd # Default database is nr
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy

# RUN PIPELINE

for directory in $(ls -d $data_dir/*/); do # Loop through all directories in data_dir

  if [[ -d ${directory}/contigs_out ]]; then
    contigs_out=$(ls ${directory}/contigs_out | grep contigs.fa)
    contig_file_length=$(wc -l $directory/contigs_out/$contigs_out | awk '{print $1}') # Get number of lines in contig file
    [[ contig_file_length -ne 0 ]] && echo contigs found for ${directory}, skipping processing  && continue # Skip if contigs found
  fi
  
  touch ${directory}/log.txt # Create log file
  echo "Processing $directory" >> ${directory}/log.txt # Print directory name to log file

  ## TRIM READS
  ${install_path}ReadTrimming.sh -r $directory -o $directory -t $threads # Run read trimming

  ## DEDUPLICATE READS
  [[ -z $no_prinseq ]] && ${install_path}ReadDeduplication.sh -r $directory -o $directory -t $threads || echo "Skipping prinseq for sample ${directory}" # Run read deduplication

  ## FILTER HOST READS  
  [[ -z $genome ]] && [[ ! -z $metadata ]] && ${install_path}HostFiltering.sh -r $directory  -h $host_ref_genome_dir -o $directory -t $threads -m $metadata -s $species_column # Run host filtering if a metadata file was provided
  [[ -z $metadata ]]  && [[ ! -z $genome ]] && ${install_path}HostFiltering.sh -r $directory  -o $directory -t $threads -g $genome # Run host filtering if a host genome was provided
  [[ -z $metadata ]] && [[ -z $genome ]] && echo No host genome or metadata file specified, skipping host filtering # If neither a metadata file or host genome was provided, skip host filtering

  ## ASSEMBLE READS
  ${install_path}Assembly.sh -r $directory -o $directory -t $threads --no_host_filtering $no_host_filtering # Run assembly


  ## COMPRESS FASTQ FILES
  ${install_path}CompressFASTQs.sh -r $directory -t $threads # Compress FASTQ files

  ## CHECK PROGRESS
  ${install_path}CheckProgress.sh -r $data_dir # Check progress
done

## RUN DIAMOND
${install_path}RunDiamond_lowRAM.sh -d $data_dir -t $threads -b $database -k $krona_tools_db # Run diamond on all contigs (for perfomrance reasons, this is done on a file of all contigs concatenated together)