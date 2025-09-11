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
      -f|--score_filter) score_filter="$2"; shift;;
      -w|--which_db) which_db="$2"; shift;;

    esac
    shift
done

install_path=$(dirname -- "$0")/

[[ -z $database ]] && database=/db/diamond/nr.dmnd # Default database is nr
[[ -z $threads ]] && threads=10 # Default number of threads is 10
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy
[[ -z $score_filter ]] && score_filter=10 && echo Filtering blast results using default cutoff of 1x10-10 # Default score filter is 0
[[ -z $which_db ]] && which_db=diamond # Default database is nr


for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    pigz -d $directory/*_unmap_*.fastq.gz
    cat $directory/*_unmap_1.fastq $directory/*_unmap_2.fastq > ${data_dir}/${accession}/${accession}_unmap.fastq
    seqtk seq -a ${data_dir}/${accession}/${accession}_unmap.fastq > ${data_dir}/${accession}/${accession}_unmap.fasta
    rm ${data_dir}/${accession}/${accession}_unmap.fastq
    reads=${data_dir}/${accession}/${accession}_unmap.fasta
    [[ ! -f $reads ]] && echo WARNING: reads not found for file ${accession}, skipping this sample && continue
    sed -i "s/^>/>${accession}|/" $reads
    cat $reads >> ${data_dir}/all_unmapped_reads.fa
    pigz $directory/*_unmap_*.fastq
done

diamond blastx -d $database -q ${data_dir}/all_unmapped_reads.fa -o ${data_dir}/all_reads.${file_string}diamond.txt --outfmt 6 -p $threads 


for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    cat ${data_dir}/all_reads.${file_string}diamond.txt | grep ${accession} > ${data_dir}/${accession}/${accession}.reads.${file_string}diamond.txt
    ktImportBLAST ${data_dir}/${accession}/${accession}.reads.${file_string}diamond.txt -o ${data_dir}/${accession}/${accession}.reads.${file_string}diamond.krona.html -tax $krona_tools_db
    accession=$(basename $directory)
    python ${install_path}/ParseKrona.py -k ${data_dir}/${accession}/${accession}.reads.${file_string}diamond.krona.html -c ${data_dir}/${accession}/${accession}_unmap.fasta -o ${data_dir}/${accession} -w $which_db -s $score_filter -r true
    tail -n+2 ${data_dir}/${accession}/${which_db}_${accession}_read_level_virus_hits.csv  >> ${data_dir}/${which_db}_read_level_virus_hits.csv
done