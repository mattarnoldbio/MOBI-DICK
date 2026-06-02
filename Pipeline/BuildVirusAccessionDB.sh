#!/usr/bin/bash

# Matt Arnold 2023

# This script takes a list of virus taxids and a directory of accession files
# and outputs a file containing all accessions for the given taxids

# Input:    -v --virus_ids: a file containing the list of all virus taxids
#           -a --accession_dir: a directory containing all accession files
#           -o --out_file: the file to write the output to


# Parse command line arguments
while [[ $# -gt 0 ]]
    do
        case $1 in
            -v|--virus_ids) virus_ids="$2"; shift;;
            -a|--accession_file) accession_file="$2"; shift;;
            -o|--out_file) out_file="$2"; shift ;;
        esac
        shift
    done


# query_id () {
#     for accessionfile in $accession_dir/*; do
#         #echo Searching $accessionfile for $1
#         taxid__=$1
#         taxid_=$(echo "\b"${taxid__}"\b")
#         grep -E $taxid_ $accessionfile >> $out_file
#     done
# }



while IFS= read -r accession_row; do
    read -ra taxid <<<"$accession_row"
    accession=${taxid[0]}
    taxid_=${taxid[1]}
    if [[ $taxid_ == $taxid__ ]]; then
        [[ $check -eq 1 ]] && echo -e $accession'\t'$taxid_ >> $out_file
        continue
    fi
    taxid__=$taxid_ #$(echo $'\s'$taxid_$'\s')
    if grep -qE '\s'${taxid__}$ $virus_ids; then
        echo -e $accession'\t'$taxid_  >> $out_file
        check=1
    else
        check=0
    fi
done < $accession_file #/home/arno01m/4tHardDrive/db/ncbi_taxonomy/accession2taxid/dummy.txt #prot.accession2taxid.FULL.1_short