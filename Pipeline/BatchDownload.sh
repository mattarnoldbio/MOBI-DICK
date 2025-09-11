#!/usr/bin/bash

# Matt Arnold 2023

# This script downloads the SRA files from NCBI using the SRA toolkit and converts them to fastq files using the fasterq-dump tool
# The IDs are shuffled at the start of the script and checked before downloading, so download can be parellelised by running multiple 
# instances of this script (making sure to send the output to different files)

# Input: transcriptome csv file with SRA IDs in the column n columns from the last (change n variable if needed)
# i.e. n=1 means ID is in the last column, n=2 means ID is in the second last column etc.

# Tip: This will take a long time to run, so it is recommended to run it in the background using nohup

time {
    n=1 # specify column number of SRA ID in the csv file (numbering is inverse, from the last column)
    transcriptome_file=$1 # specify transcriptome csv file as input to this script

    transcriptomes=$( tail -n +2 $transcriptome_file | rev | cut -f$n -d, | rev | shuf ) # extract SRA IDs from the csv file


    for record in $transcriptomes # loop through the SRA IDs
        do if [ -d $record ]
            then echo $record already downloaded
            continue
        fi
        echo Processing transcriptome with ID: $record
        echo Fetching SRA file
        prefetch -p $record # download SRA file
        [[ -d $record ]] || mkdir $record
        cd ./$record # change directory to the SRA file
        echo Downloading fastq file\(s\) 
        fasterq-dump -p --split-3 -e 8 $record # download fastq file(s) from SRA file
        
        gzip *.fastq # compress fastq file(s)
        rm ./*.sra
        echo Record $record done
        cd ../ # change directory back to the main directory
    done
}