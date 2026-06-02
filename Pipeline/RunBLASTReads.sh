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
      -x|--xtra_fast) xtra_fast=$2; shift;;

    esac
    shift
done

install_path=$(dirname -- "$0")/

[[ -z $database ]] && database=/db/diamond/nr.dmnd # Default database is nr
[[ -z $threads ]] && threads=8 # Default number of threads is 10
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy
[[ -z $score_filter ]] && score_filter=10 && echo Filtering blast results using default cutoff of 1x10-10 # Default score filter is 0
[[ -z $which_db ]] && which_db=diamond # Default database is nr
[[ -z $file_string ]] && file_string='unmapped_' # Default file string is 'unmapped_
[[ -z $xtra_fast ]] && xtra_fast=false


if [[ ! -f $data_dir/all_unmapped_reads.fa ]] || [[ -f $data_dir/all_unmapped_reads.fa ]] && [[ ! -s $data_dir/all_unmapped_reads.fa ]] ; then
    for directory in $(ls -d $data_dir/*/); do
            accession=$(basename $directory)
            pigz -d $directory/*_unmap*.fastq.gz
            paired_end=false
            for f in $(ls ${directory}/*unmap_2*); do 
                if [ -e "$f" ]; then
                    cat $directory/*_unmap_1.fastq $directory/*_unmap_2.fastq > ${data_dir}/${accession}/${accession}_unmap.fastq
                    paired_end=true
                else
                    echo "Single end/Nanopore reads detected, skipping concatenation"
                fi
            done 
            #cat $directory/*_unmap_1.fastq $directory/*_unmap_2.fastq > ${data_dir}/${accession}/${accession}_unmap.fastq

            seqtk seq -a ${data_dir}/${accession}/${accession}_unmap.fastq > ${data_dir}/${accession}/${accession}_unmap.fasta
            [[ $paired_end == true ]] && rm ${data_dir}/${accession}/${accession}_unmap.fastq
            reads=${data_dir}/${accession}/${accession}_unmap.fasta
            [[ ! -f $reads ]] && echo WARNING: reads not found for file ${accession}, skipping this sample && continue
            sed -i "s/^>/>${accession}|/" $reads
            cat $reads >> ${data_dir}/all_unmapped_reads.fa
            echo "Finished processing ${accession} with" $(grep -c "^>" $reads) "reads"
            pigz $directory/*_unmap*.fastq
    done
fi

[[ -f ${data_dir}/all_reads_hits.fa ]] && mkdir -p ${data_dir}/.archive && mv ${data_dir}/all_reads_hits.fa ${data_dir}/.archive/all_reads_hits.fa.$(date +%Y%m%d%H%M%S)

if [[ $which_db == blastn ]] ; then # || [[ $(($(wc -l all_reads_hits.fa | cut -f1 -d" "))) -eq 0 ]] ; then

    touch ${data_dir}/all_reads_hits.fa
    unmapped_reads=${data_dir}/all_unmapped_reads.fa
    for directory in $(ls -d $data_dir/*/); do # For each sample
        accession=$(basename $directory) # Get the accession number
        if [[ ! -f ${directory}/${accession}_no_virus_hits.txt ]]; then # If there are virus hits
             
            virus_hits=$(cat ${directory}/diamond_${accession}_read_level_all_virus_hits.csv | cut -f 2 -d "," | tail -n+2 ) # Extract the information about which sequences contained hits to viruses
            echo "Processing ${accession} with $(echo $virus_hits | wc -w) hits"
            for hit in $virus_hits; do # For each hit
                read=$(grep "${accession}|${hit}" ${unmapped_reads} -A 1 -m1) # Extract the contig header and sequence from the contigs file
                echo "$read" >> ${data_dir}/all_reads_hits.fa # Add the accession number to the contig header and append to the all_contigs_hits.fa file
            done
        fi
    done
fi


if [[ $xtra_fast == true ]]; then
    blocks=6
    chunks=1
elif [[ $xtra_fast == false ]]; then
    blocks=1
    chunks=4
fi



if [[ $which_db == diamond ]]; then
    diamond blastx -d $database -q ${data_dir}/all_unmapped_reads.fa -o ${data_dir}/all_reads.${file_string}diamond.txt --outfmt 6 -p $threads --sensitive -b $blocks -c $chunks
elif [[ $which_db == blastn ]]; then
    database='/db/blast_v5/nt'
    blastn -db $database -query ${data_dir}/all_reads_hits.fa -out ${data_dir}/all_reads.${file_string}blastn.txt -outfmt 6 -num_threads $threads
fi
#diamond blastx -d $database -q ${data_dir}/all_unmapped_reads.fa -o ${data_dir}/all_reads.${file_string}diamond.txt --outfmt 6 -p $threads 

[[ -f ${data_dir}/${which_db}_read_level_all_virus_hits.csv ]] && [[ which_db==diamond ]] && mkdir -p ${data_dir}/.archive && mv ${data_dir}/${which_db}_read_level_all_virus_hits.csv ${data_dir}/.archive/${which_db}_read_level_all_virus_hits.csv.$(date +%Y%m%d%H%M%S)

for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    cat ${data_dir}/all_reads.${file_string}${which_db}.txt | grep ${accession} > ${data_dir}/${accession}/${accession}.reads.${file_string}${which_db}.txt
    ktImportBLAST ${data_dir}/${accession}/${accession}.reads.${file_string}${which_db}.txt -o ${data_dir}/${accession}/${accession}.reads.${file_string}${which_db}.krona.html -tax $krona_tools_db
    accession=$(basename $directory)
    python ${install_path}/ParseKrona.py -k ${data_dir}/${accession}/${accession}.reads.${file_string}${which_db}.krona.html -c ${data_dir}/${accession}/${accession}_unmap.fasta -o ${data_dir}/${accession} -w $which_db -s $score_filter --read_mode
    tail -n+2 ${data_dir}/${accession}/${which_db}_${accession}_read_level_all_virus_hits.csv  >> ${data_dir}/${which_db}_read_level_all_virus_hits.csv
done