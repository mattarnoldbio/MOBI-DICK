#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to concatenate all contigs from a run and blast them against a chosen database (default = nr)
# Input:    -d --data_dir  path to directory containing raw data
#           -b --database  database to be used (default = nr)
#           -t --threads  number of cpu threads to use

# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -b|--database) database="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -k|--krona_tools_db) krona_tools_db="$2"; shift;;
      -s|--file_string) file_string="$2"; shift;;
    esac
    shift
done

install_path=$(dirname -- "$0")/

[[ -z $database ]] && database=/db/diamond/nr.dmnd # Default database is nr
[[ -z $threads ]] && threads=10 # Default number of threads is 10
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy

# if [[ ! -f ${data_dir}/all_contigs.fa ]]; then
#   for directory in $(ls -d $data_dir/*/); do
#       accession=$(basename $directory)
#       contigs=${directory}/contigs_out/${accession}.contigs.fa
#       [[ ! -f $contigs ]] && echo WARNING: contigs not found for file ${accession}, skipping this sample && continue
#       cat $contigs | sed "s/^>/>${accession}|/" >> ${data_dir}/all_contigs.fa
#   done
# fi


# diamond blastx -d $database -q ${data_dir}/all_contigs.fa -o ${data_dir}/all_contigs.${file_string}diamond.txt --outfmt 6 -p $threads -b 3

for directory in $(ls -d $data_dir/*/); do
    echo CHANGE THIS BACK!!!!
    accession=$(basename $directory)
    cat ${data_dir}/all_contigs.${file_string}diamond | grep ${accession} > ${data_dir}/${accession}/${accession}.${file_string}diamond.txt
    ktImportBLAST ${data_dir}/${accession}/${accession}.${file_string}diamond.txt -o ${data_dir}/${accession}/${accession}.${file_string}diamond.krona.html -tax $krona_tools_db
done