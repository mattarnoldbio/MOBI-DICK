#!/usr/bin/bash

# Input: -d --data_dir path to directory containing raw data
#        -o --outdir path to output directory

# Parse arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -o|--outdir) outdir="$2"; shift;;
      -r|--read_mode) read_mode="$2"; shift;;
    esac
    shift
done

[[ -z $outdir ]] && outdir=$data_dir
[[ $read_mode == "true" ]] && filestring="_read_level" || filestring=""

which_db=blastn

count=0

for directory in $(ls -d $data_dir/*/); do
    count=$((count+1))
    echo Processing $directory $count of $(ls -d $data_dir/*/ | wc -l)
    accession=$(basename $directory)
    tail -n+2 ${data_dir}/${accession}/${which_db}_${accession}${filestring}_non_virus_hits.csv >> $outdir/${which_db}${filestring}_non_virus_hits.csv
done