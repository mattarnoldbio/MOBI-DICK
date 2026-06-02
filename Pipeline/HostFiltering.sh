#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to filter the host reads from the raw reads
# Two ways of running this script are supported:
#       1. If a metadata file is provided, the script will use the species name in the metadata file to find the corresponding host reference genome
#       2. If no metadata file is provided, the script will use the host reference genome specified in the command line arguments

# Input:    -g --genome    path to host reference genome (or multiple concatenated genomes)
#                          (GENOME MUST BE INDEXED WITH BOWTIE2, see Pipeline/IndexHostGenome.sh)
#
#           -r --raw_reads path to trimmed raw reads
#           -o --out_file  output file name
#           -m --metadata  path to metadata file containing species names for each accession number
#           -s --species_column column number of species in metadata file
#           -h --host_ref_genome_dir path to directory containing host reference genomes

# Output:   File containing reads NOT mapping to the host genome

# Parse command line arguments
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -g|--genome) genome="$2"; shift;;
      -r|--raw_reads) raw_reads="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
      -m|--metadata) metadata="$2"; shift;;
      -s|--species_column) species_column="$2"; shift;;
      -h|--host_ref_genome_dir) host_ref_genome_dir="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -i|--install_path) install_path="$2"; shift;;
    esac
    shift
done

# [[ $(conda info --envs | tail -n+3 | grep SRA_mining | wc -l) -e 0 ]] && ~/SRA_mining/Pipeline/Install.sh # Check if conda environment exists, if not, create it

# conda activate SRA_mining

[[ -z $threads ]] && threads=8 # Default number of threads is 10
[[ -z $install_path ]] && install_path=$(dirname -- "$0")/ # Get path to install directory
[[ -z species_column ]] && species_column=3 # Default species column is 3 (this happened to be which column it was in the metadata file I was using)

accession=$(basename `readlink -f $raw_reads`) # Extract accession number from path

reference=""

if [[ -z $genome ]]; then
    echo Running in SRA mode, using possible reference genomes in $host_ref_genome_dir
    ref_genomes=($(ls $host_ref_genome_dir)) # Get host reference genomes


    # If a metadata file was provided, use the species name to find the corresponding host reference genome
    if [[ $metadata ]] && [[ -z $species_column ]] 
        then species_column=3 # default species column is 3 (this happened to be which column it was in the metadata file I was using)
        echo "No species column specified, using column 3 as default" # Warning for others to check this is correct
        else
        if [[ $metadata ]] 
            then echo "Species column specified as $species_column" # If a species column was specified, print it
        fi
    fi


    if [[ ! -z $metadata ]] # If a metadata file was provided, use the species name to find the corresponding host reference genome
        then [[ -z $host_ref_genome_dir ]] && echo "No host reference genome directory specified. Please specify a directory containing the indexed host reference genomes you wish to use using the -h flag" && exit 1 # If no host reference genome directory was specified, exit with error
        species=$(cat $metadata | grep $accession | cut -f$species_column -d, ) # Find the species associated with the accession number
        echo Species is $species
        taxonomy=$(${install_path}famdb.py lineage -a "$species") # Get the taxonomy of the species
        for host_taxon in ${ref_genomes[@]}; do
            if [[ $(echo $taxonomy | tr ' ' '_' | grep -c $host_taxon) -ne 0 ]] # If the taxonomy contains the host taxon, use that reference genome
                then reference=$(ls ${host_ref_genome_dir}/${host_taxon}/*.4.bt2 | rev | cut -d"/" -f1 | rev | cut -d"." -f1)  # Set the reference genome
                echo "Using $reference as host reference genome"
                break
            fi
        done
        [[ -z $reference ]] && echo "No host reference genome found for $accession . Make sure you have a drectory containing the indexed reference genome you wish to use, and that this is saved in the directory you issued using the -h flag. This directory must be named as the species you are mappping to, or as a higher taxnomic level in the lineage of the species of interest (e.g. a directory called /Phocidae for the species Halichoerus grypus)" && exit 1 # If no host reference genome was found, exit with error
    fi
    else
    echo Running in metagenomics mode, using $genome as host reference genome
fi


a=($(ls $raw_reads)) # Get files in specified directory


echo Aligning $accession to host genome

paired=false

for file in ${a[@]}; do # check if paired end or single end
    #[[ ${file: -16} == "_1_trimmed.fq.gz" ]] || [[ ${file: -13} == "_1_trimmed.fq" ]] || 
    [[ ${file: -11} == "val_1.fq.gz" ]] || [[ ${file: -8} == "val_1.fq" ]] && R1=$raw_reads/$file
    #[[ ${file: -16} == "_2_trimmed.fq.gz" ]] || [[ ${file: -13} == "_2_trimmed.fq" ]] || 
    [[ ${file: -11} == "val_2.fq.gz" ]] || [[ ${file: -8} == "val_2.fq" ]] && paired=true && R2=$raw_reads/$file && echo "Paired-end reads found"
done 

if [[ -z $R1 ]]; then # If no reads have been found yet, check for reads with different name format
    for file in ${a[@]}; do 
        [[ ${file: -14} == "_trimmed.fq.gz" ]] || [[ ${file: -11} == "_trimmed.fq" ]] && R1=$raw_reads/$file && echo "R1 file found: $R1"
    done 
fi

[[ -z $R1 ]] && { printf '%s\n' "No trimmed reads found for file $accession" >&2; exit 1; } # If no reads have been found, exit with error



[[ $paired == true ]] && echo Processing R1: $R1 R2: $R2 || echo Processing unpaired reads: $R1

[[ ! -z $reference ]] && genome=${host_ref_genome_dir}/${host_taxon}/$reference

echo genome is $genome

# UNCOMMENT HERE 


# Run the alignment to the host reference genome using bowtie2 (conditional statement for paired end or single end reads)
[[ $paired == true ]] && (bowtie2 -p $threads --local -x $genome -1 $R1 -2 $R2 -S ${raw_reads}/${accession}.sam)  || (bowtie2 -p $threads --local -x $genome -U $R1  -S ${raw_reads}/${accession}.sam)
#[[ $paired ]] && bowtie2 -p 10 --local -x $genome -1 $R1 -2 $R2 --quiet || bowtie2 -p 10 --local -x $genome -U $R1  

# Post-process the alignment file and format it for downstream analysis
samtools view -@ $threads -bS ${raw_reads}/${accession}.sam > ${raw_reads}/${accession}.bam
samtools sort -@ $threads ${raw_reads}/${accession}.bam -o ${raw_reads}/${accession}.bam
samtools index ${raw_reads}/${accession}.bam

# Extract the reads that did not map to the host genome
[[ $paired == true ]] && bam2fastq -f --unaligned --no-aligned -o ${raw_reads}/${accession}_unmap#.fastq ${raw_reads}/${accession}.bam || bam2fastq -f --unaligned --no-aligned -o ${raw_reads}/${accession}_unmap.fastq ${raw_reads}/${accession}.bam

rm ${raw_reads}/${accession}.sam

# Write the number of host reads to a log file
echo "host-filtered reads" >> ${raw_reads}/log.txt
[[ $paired == true ]] && expr `(wc -l ${raw_reads}/*_unmap_1.fastq |cut -f1 -d " ")` / 4 >> ${raw_reads}/log.txt || expr `(wc -l ${raw_reads}/*_unmap.fastq |cut -f1 -d " ")` / 4 >> ${raw_reads}/log.txt


#conda deactivate_
