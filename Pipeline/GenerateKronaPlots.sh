#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to generate krona plot files for the outputs of a diamond/blast run

# Input:    -d --data_dir  path to directory containing raw data
#           -k --krona_tools_db  path to krona tools taxonomy database

# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
        -d|--data_dir) data_dir="$2"; shift;;
        -k|--krona_tools_db) krona_tools_db="$2"; shift;;
        -w|--which_db) which_db="$2"; shift;;
    esac
    shift
done

[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy



for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    cat ${data_dir}/all_contigs_hits.${which_db}.txt | grep ${accession} > ${data_dir}/${accession}/${accession}.${which_db}.txt
    ktImportBLAST ${data_dir}/${accession}/${accession}.${which_db}.txt -o ${data_dir}/${accession}/${accession}.${which_db}.krona.html -tax $krona_tools_db
done