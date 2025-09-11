#!/usr/bin/bash

# Matt Arnold 2023

# This script takes a list of virus taxids and a directory of accession files

# Input:    -v --virus_ids: a file containing the list of all virus taxids
#           -a --accession_dir: a directory containing all accession files
#           -o --out_file_string: the file to write the output to

# Parse command line arguments
while [[ $# -gt 0 ]]
    do
        case $1 in
            -v|--virus_ids) virus_ids="$2"; shift;;
            -a|--accession_dir) accession_dir="$2"; shift;;
            -o|--out_file_string) out_file_string="$2"; shift ;;
        esac
        shift
    done

for file in /${accession_dir}/prot.accession2taxid.FULL.*;  do
    #echo $file
    #touch ${file}_viruses
    /home/arno01m/Repos/SRAmining/Pipeline/BuildVirusAccessionDB.sh -v $virus_ids -a $file -o ${file}_viruses &
done