#!/usr/bin/bash

# Matt Arnold 2023

# This script filters blast hits for viruses

# Input:    -b --blast_file: the blast file to be filtered
#           -v --virus_ids: a file containing the list of all virus taxids
#           -o --out_file: the file to write the output to

# Parse command line arguments

while [[ $# -gt 0 ]]
    do
        case $1 in
            -b|--blast_file) blast_file="$2"; shift;;
            -v|--virus_ids) virus_ids="$2"; shift;;
            -o|--out_file) out_file="$2"; shift ;;
        esac
        shift
    done


# Filter blast file for virus hits

if [[ ! -f ${blast_file}.krona.tab ]]; then
    echo Finding taxids for BLAST hits 
    #ktClassifyBLAST $blast_file -o ${blast_file}.krona.tab
else
    echo ${blast_file}.krona.tab already exists, skipping this step
fi

head $virus_ids

while IFS= read -r row; do
    
    read -ra hit <<<"$row"
    id=${hit[0]}
    taxid=${hit[1]}
    #echo searching for taxid $taxid
    score=${hit[2]}
    if grep -qE '\s'${taxid}$ $virus_ids; then
        echo -e $id'\t'$taxid'\t'$score >> $out_file
    fi
done < ${blast_file}.krona.tab