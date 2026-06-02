#!/usr/bin/bash


for directory in ./*/; do 
    [[ -f ${directory}/log.txt ]] || continue
    accession=$(basename $directory)
    echo -ne ${accession}"\t">> readcounts.tsv
    arr=($(egrep ^[0-9]  ${directory}/log.txt))
    for item in "${arr[@]}"; do   
        echo -ne "$item""\t" >> readcounts.tsv
    done
    echo -e $(( $(wc -l ${directory}/contigs_out/${accession}.contigs.fa | cut -f1 -d" ") / 2 ))"\t" >> readcounts.tsv
    #echo -ne "\n" >> readcounts.tsv 
done
