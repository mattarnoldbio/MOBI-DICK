#/usr/bin/bash

# Matt Arnold 2023

# This script is to trim the raw reads using trim_galore

# Input:    -r --raw_reads path to raw reads
#           -o --out_file  output file name

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -r|--raw_reads) raw_reads="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
    esac
    shift
done

a=($(ls $raw_reads))

for file in ${a[@]}; do # check if paired end or single end
    pigz -d -f -p $threads $raw_reads/$file
done 


for file in ${a[@]}; do # check if trimmed read files already exist
    #echo ${file: -13}
    [[ ${file: -13} == "trimmed.fq.gz" ]] || [[ ${file: -10} == "trimmed.fq" ]] || [[ ${file: -8} == "val_1.fq" ]] || [[ ${file: -11} == "val_1.fq.gz" ]] && { printf '%s\n' "Reads already trimmed" >&2; exit 1; }
done 

paired=false # Set paired to false by default


a=($(ls $raw_reads))
accession=$(basename `readlink -f $raw_reads`) # Extract accession number from path

for file in ${a[@]}; do # check if paired end or single end
    #echo file is $file
    [[ ${file: -11} == "_1.fastq.gz" ]] || [[ ${file: -8} == "_1.fastq" ]] || [[ ${file: -13} == "_R1_001.fastq" ]] && R1=$raw_reads/$file
    [[ ${file: -11} == "_2.fastq.gz" ]] || [[ ${file: -8} == "_2.fastq" ]] || [[ ${file: -13} == "_R2_001.fastq" ]]  && R2=$raw_reads/$file && paired=true && break 
done 

if [[ -z $R1 ]]; then # If no reads have been found yet, check for reads with different name format
    for file in ${a[@]}; do 
      [[ ${file: -9} == ".fastq.gz" ]] || [[ ${file: -6} == ".fastq" ]] && R1=$raw_reads/$file
    done 
fi

[[ $paired == true ]] && echo Processing reads as paired end || echo Processing reads as single end. If this is incorrect, please check read names.

echo "Raw reads" >> ${raw_reads}/log.txt
expr `(wc -l $R1 |cut -f1 -d " ")` / 4 >> ${raw_reads}/log.txt  # Print number of reads to log file

[[ $paired == true ]] && trim_galore $R1 $R2 -o $out_file -j $threads --paired -q 25 ||  trim_galore $R1 -o $out_file -j $threads -q 25 # Trim reads 

echo "trim_galore reads" >> ${raw_reads}/log.txt
[[ $paired == true ]]  && expr `(wc -l ${raw_reads}/*1_val_1.fq |cut -f1 -d " ")` / 4 >> ${raw_reads}/log.txt || expr `(wc -l ${raw_reads}/*_trimmed.fq |cut -f1 -d " ")` / 4 >> ${raw_reads}/log.txt # Print number of trimmed reads to log file
