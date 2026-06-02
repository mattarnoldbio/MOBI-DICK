#!/usr/bin/bash

# Matt Arnold 2023

# This script is to produce user-reeadable output tables from the Krona plots

# Input:    -d --data_dir path to data directory

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -i|--install_path) install_path="$2"; shift;;
      -w|--which_db) which_db="$2"; shift;;
      -s|--score_filter) score_filter="$2"; shift;;
      -r|--read_mode) read_mode="$2"; shift;;
    esac
    shift
done

[[ -z $which_db ]] && which_db=diamond # Default database is nr
[[ -z $score_filter ]] && score_filter=10 && echo Filtering blast results using default cutoff of 1x10-10 # Default score filter is 0
[[ -z $read_mode ]] && read_mode=False # Default read mode is false

for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    if [[ $read_mode != True ]]; then 
      python ${install_path}/ParseKrona.py -k ${data_dir}/${accession}/${accession}.${which_db}.krona.html -c ${data_dir}/${accession}/contigs_out/${accession}.contigs.fa -o ${data_dir}/${accession}/ -a $accession -w $which_db -s $score_filter -a $accession --contig_mode
    elif [[ $read_mode == True ]]; then
      python ${install_path}/ParseKrona.py -k ${data_dir}/${accession}/${accession}.${which_db}.krona.html -c ${data_dir}/${accession}/contigs_out/${accession}.contigs.fa -o ${data_dir}/${accession}/ -a $accession -w $which_db -s $score_filter  --read_mode
    fi
    tail -n+2 ${data_dir}/${accession}/${which_db}_${accession}_all_virus_hits.csv >> ${data_dir}/${which_db}_all_virus_hits.csv
done

# python ${install_path}/JoinResults.py -d $data_dir