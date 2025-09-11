#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to blast a fasta of contigs against a chosen database (default = nt)

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
      -s|--score_filter) score_filter="$2"; shift;;
    esac
    shift
done

install_path=$(dirname -- "$0")/

[[ -z $database ]] && database=/db/blast_v5/nt # Default database is nr
[[ -z $threads ]] && threads=8 # Default number of threads is 8
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy
[[ -z $score_filter ]] && score_filter=-10 # Default score filter is -10

if [[ ! -f ${data_dir}/all_contigs.fa ]]; then
    for directory in $(ls -d $data_dir/*/); do
        accession=$(basename $directory)
        contigs=${directory}/contigs_out/${accession}.contigs.fa
        [[ ! -f $contigs ]] && echo WARNING: contigs not found for file ${accession}, skipping this sample && continue
        cat $contigs | sed "s/^>/>${accession}|/" >> ${data_dir}/all_contigs.fa
    done
fi

if [[ ! -f ${data_dir}/all_contigs_hits.fa ]] || [[ $(($(wc -l all_contigs_hits.fa | cut -f1 -d" "))) -eq 0 ]] ; then
    touch ${data_dir}/all_contigs_hits.fa

    for directory in $(ls -d $data_dir/*/); do # For each sample
        accession=$(basename $directory) # Get the accession number
        if [[ ! -f ${directory}/${accession}_no_virus_hits.txt ]]; then # If there are virus hits
            contigs=${directory}/contigs_out/${accession}.contigs.fa # Get the contigs file
            virus_hits=$(cat ${directory}/diamond_${accession}_virus_hits.csv | cut -f 2 -d "," | tail -n+2 ) # Extract the information about which sequences contained hits to viruses
            for hit in $virus_hits; do # For each hit
                contig=$(grep "${hit} " ${directory}/contigs_out/${accession}.contigs.fa -A 1) # Extract the contig header and sequence from the contigs file
                echo "$contig" | sed "s/^>/>${accession}|/" >> ${data_dir}/all_contigs_hits.fa # Add the accession number to the contig header and append to the all_contigs_hits.fa file
            done
        fi
    done
fi

blastn -db $database -query ${data_dir}/all_contigs_hits.fa -out ${data_dir}/all_contigs_hits.blastn.txt -outfmt 6 -num_threads $threads

for directory in $(ls -d $data_dir/*/); do
    accession=$(basename $directory)
    cat ${data_dir}/all_contigs_hits.blastn.txt | grep ${accession} > ${data_dir}/${accession}/${accession}.blastn.txt
    ktImportBLAST ${data_dir}/${accession}/${accession}.blastn.txt -o ${data_dir}/${accession}/${accession}.blastn.krona.html -tax $krona_tools_db
done

## POST PROCESS BLAST RESULTS
${install_path}PostProcessKrona.sh -d $data_dir -i ${install_path} -w blastn -s $score_filter # Post process diamond results
