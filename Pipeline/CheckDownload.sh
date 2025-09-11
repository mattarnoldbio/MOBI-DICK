#!/usr/bin/bash
# Matt Arnold 2023

# This script is used to check the files have correctly downloaded from the SRA

# Input:    -a --accession accession number
#           -d --download_dir path to directory containing downloaded files

# Parse args
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -a|--accession) accession="$2"; shift;;
        -d|--download_dir) download_dir="$2"; shift;;
    esac
    shift
done

function download_missing {
    echo Processing transcriptome with ID: $record
    echo Fetching SRA file
    prefetch -p $record # download SRA file
    [[ -d $download_dir/$record ]] || mkdir $download_dir/$record
    cd $download_dir/$record # change directory to the SRA file
    echo Downloading fastq file\(s\) 
    fasterq-dump -p --split-3 -e 8 $record # download fastq file(s) from SRA file
    
    gzip *.fastq # compress fastq file(s)
    rm ./*.sra
    echo Record $record done
    cd $download_dir # change directory back to the main directory
}

echo dir is $download_dir
for directory in $(ls -d $download_dir/*); do # Loop through all directories in download directory
    
    correct_files=0
    for file in $(ls $directory); do # Loop through all files in directory
        [[ $file == *.fq* ]] || [[ $file == *.fastq* ]] && correct_files=1 && break # If a fastq file is found, set correct_files to 1
    done
    record=$(basename $directory) # Extract accession number from path
    [[ $correct_files == 1 ]] && echo $record already downloaded # If a fastq file was found, print message
    [[ $correct_files == 0 ]] && download_missing  # If no fastq files were found, print error
done