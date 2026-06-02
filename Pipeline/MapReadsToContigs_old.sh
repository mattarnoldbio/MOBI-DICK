#!/usr/bin/bash

# Matt Arnold 2024

# This script maps the raw reads back to the contigs

# Input:    -d --data_dir  path to directory containing the raw reads and contigs
#           -t --threads   number of threads to use
#           -c --contigs_csv path to csv file containing contigs to map reads to
#                               (must be in MOBI-DICK filtered hit format, see viral_hits.csv
#                               for example)



# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -c|--contigs_csv) contigs_csv="$2"; shift;;
    esac
    shift
done

[[ -z $threads ]] && threads=8 # Default number of threads is 8
[[ -z $contigs_csv ]] && contigs_csv=${data_dir}/vertebrate_viral_hits.csv

mkdir -p ${data_dir}/contig_mappings

declare -a read_files=()


echo reads_mapped > ${data_dir}/contig_mappings/tmp

while IFS="" read -r hit || [ -n "$hit" ]
do
    sequence=$(echo ${hit} | rev | cut -d, -f1 | rev)
    [[ $sequence == "sequence" ]] && continue
    accession=$(echo ${hit} | cut -d, -f1)
    contig=$(echo ${hit} | cut -d, -f2)
    id=${accession}_${contig}

    mkdir -p ${data_dir}/contig_mappings/${accession}/${id}

    echo ">${id}" > ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.fasta
    echo $sequence >> ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.fasta

    bowtie2-build -q ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.fasta ${data_dir}/contig_mappings/${accession}/${id}/${id}

    R1=${data_dir}/${accession}/$(ls ${data_dir}/${accession} | grep -E "1.fastq|1.fq" | grep -vE "val|.txt|unmap")
    R2=${data_dir}/${accession}/$(ls ${data_dir}/${accession} | grep -E "2.fastq|2.fq" | grep -vE "val|.txt|unmap")

    echo R1 is $R1, R2 is $R2
    [[ ${read_files[-1]} != $R2 ]] && read_files+=($R1 $R2)

    for file in {$R1,$R2}; do
       [[ ${file: -3} == ".gz" ]] && pigz -d $file -p $threads
    done

    R1=${R1%.gz}
    R2=${R2%.gz}
    
    bowtie2 -x ${data_dir}/contig_mappings/${accession}/${id}/${id} -1 $R1 -2 $R2 -S ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.sam -p $threads
    samtools view -@ $threads -bS ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.sam > ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.bam
    samtools sort -@ $threads ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.bam -o ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.bam
    samtools index ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.bam
    samtools stats ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.bam | grep ^SN | cut -f 2-  > ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.stats

    grep "^reads mapped:" ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.stats | cut -f2 >> ${data_dir}/contig_mappings/tmp

    echo For $id, reads mapped: $(grep "^reads mapped:" ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.stats | cut -f2)
    rm ${data_dir}/contig_mappings/${accession}/${id}/sequence.${id}.sam
    
done < $contigs_csv


paste -d, ${contigs_csv} ${data_dir}/contig_mappings/tmp > ${data_dir}/contig_mappings/tmp_

mv ${data_dir}/contig_mappings/tmp_ ${contigs_csv}

rm ${data_dir}/contig_mappings/tmp

for file in ${read_files[@]}; do
    pigz $file -p $threads -q
done
