#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to exclude low-quality reads from the trimmed reads

# Input:    -r --reads path to raw reads
#           -o --out_file  output file name
#           -t --threads number of threads to use

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -r|--reads) reads="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
    esac
    shift
done

[[ -z $threads ]] && threads=10 # Default number of threads is 10

a=($(ls $reads))

paired=false

for file in ${a[@]}; do # check if paired end or single end
    #[[ ${file: -16} == "_1_trimmed.fq.gz" ]] || [[ ${file: -13} == "_1_trimmed.fq" ]] || 
    [[ ${file: -11} == "val_1.fq.gz" ]] || [[ ${file: -8} == "val_1.fq" ]] && R1=$reads/$file
    #[[ ${file: -16} == "_2_trimmed.fq.gz" ]] || [[ ${file: -13} == "_2_trimmed.fq" ]] || 
    [[ ${file: -11} == "val_2.fq.gz" ]] || [[ ${file: -8} == "val_2.fq" ]] && paired=true && R2=$reads/$file && echo "Paired-end reads found"
done 

if [[ -z $R1 ]]; then # If no reads have been found yet, check for reads with different name format
    for file in ${a[@]}; do 
        [[ ${file: -14} == "_trimmed.fq.gz" ]] || [[ ${file: -11} == "_trimmed.fq" ]] && R1=$reads/$file && echo "R1 file found: $R1"
    done 
fi

[[ -z $R1 ]] && { printf '%s\n' "No trimmed reads found for file $accession" >&2; exit 1; } # If no reads have been found, exit with error



[[ $paired == true ]] &&  prinseq++ -derep -fastq $R1 -fastq2 $R2 -out_good ${reads}/prinseq_good -out_bad ${reads}/prinseq_bad > ${reads}/prinseq_out.txt 2>&1 ||  prinseq++ -derep -fastq $R1 -out_good ${reads}/prinseq_good -out_bad ${reads}/prinseq_bad > ${reads}/prinseq_out.txt 2>&1 # Run prinseq


if [[ $paired == true ]]; then
    mv ${reads}/prinseq_good_1.fastq $R1
    mv ${reads}/prinseq_good_2.fastq $R2
    rm -f ${reads}/prinseq_good_1_singletons.fastq
    rm -f ${reads}/prinseq_good_2_singletons.fastq
    rm -f ${reads}/prinseq_bad_1.fastq
    rm -f ${reads}/prinseq_bad_2.fastq
    else
    mv ${reads}/prinseq_good.fastq $R1
    rm -f ${reads}/prinseq_good_singletons.fastq
    rm -f ${reads}/prinseq_bad.fastq
fi

accession=$(basename `readlink -f $reads`) # Extract accession number from path
echo accession is $accession

echo "prinseq reads" >> ${reads}/log.txt
[[ $paired == true ]]  && expr `(wc -l ${reads}/*1_val_1.fq |cut -f1 -d " ")` / 4 >> ${reads}/log.txt || expr `(wc -l ${reads}/*_trimmed.fq |cut -f1 -d " ")` / 4 >> ${reads}/log.txt # Print number of reads to log file



