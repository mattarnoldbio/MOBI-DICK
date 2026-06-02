#!/usr/bin/bash

# Matt Arnold 2023

# This script is used to assemble the transcriptomes using megaHIT

# Input:    -r --reads path to host filtered, trimmed reads
#           -o --out_file  output file name
#           -t --threads  number of cpu threads to use
#           -i --install_path path to install directory
#           -a --assembler assembler to use (megahit or spades)
#           -n --no_host_filtering  skip host filtering


# Parse command line arguments

while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -r|--reads) reads="$2"; shift;;
      -o|--out_file) out_file="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -i|--install_path) install_path="$2"; shift;;
      -a|--assembler) assembler="$2"; shift;;
      -n|--no_host_filtering) no_host_filtering=false; shift;;
      -s|--sensitive) sensitive=false; shift;;
    esac
    shift
done

#conda activate SRA_mining

accession=$(basename $reads) # Extract accession number from path
a=($(ls $reads))

paired=false # Set paired to false by default

echo Assembling reads for sample $accession. 





for file in ${a[@]}; do # check if paired end or single end
  #echo file is $file
  if [[ $no_host_filtering != true ]]
    then
      echo Processing unmapped reads
      [[ ${file: -14} == "_unmap_1.fastq" ]] || [[ ${file: -17} == "_unmap_1.fastq.gz" ]]  && R1=$reads/$file && echo R1 is $R1
      [[ ${file: -14} == "_unmap_2.fastq" ]] || [[ ${file: -17} == "_unmap_2.fastq.gz" ]] && paired=true && R2=$reads/$file && echo "Paired-end reads found. R2 is $R2"
      [[ ${file: -14} == "_unmap_1.fastq" ]] || [[ ${file: -17} == "_unmap_1.fastq.gz" ]]  && R1=$reads/$file && echo R1 is $R1
      [[ ${file: -14} == "_unmap_2.fastq" ]] || [[ ${file: -17} == "_unmap_2.fastq.gz" ]] && paired=true && R2=$reads/$file && echo "Paired-end reads found. R2 is $R2"
  elif [[ $no_host_filtering == true ]]
    then
      [[ ${file: -16} == "_1_trimmed.fq.gz" ]] || [[ ${file: -13} == "_1_trimmed.fq" ]] || [[ ${file: -11} == "val_1.fq.gz" ]] || [[ ${file: -8} == "val_1.fq" ]] && R1=$reads/$file && echo R1 is $R1
      [[ ${file: -16} == "_2_trimmed.fq.gz" ]] || [[ ${file: -13} == "_2_trimmed.fq" ]] || [[ ${file: -11} == "val_2.fq.gz" ]] || [[ ${file: -8} == "val_2.fq" ]] && paired=true && R2=$reads/$file && echo "Paired-end reads found. R2 is $R2"
  fi
done 

if [[ $paired == false ]]; then
    for file in ${a[@]}; do
        [[ ${file: -12} == "_unmap.fastq" ]] || [[ ${file: -14} == "_unmap_M.fastq" ]] && R1=$reads/$file && echo R1 is $R1
    done
fi

# Make sure --out-prefix is configured correctly to show accession; this is needed for concatenation before BLASTing
# Check this working for both single- and paired-end reads
if [[ $assembler == "megahit" ]]; then
    [[ $paired == true ]] && megahit -1 $R1 -2 $R2 -m 0.1 -t $threads -o ${reads}/contigs_out --out-prefix $accession  --min-contig-len 10 || megahit -r $R1 -m 0.1 -t $threads --out-prefix $accession -o ${reads}/contigs_out  --min-contig-len 10 # Run megahit
elif [[ $assembler == "megahit_sensitive_1" ]]; then
    [[ $paired == true ]] && megahit -1 $R1 -2 $R2 -m 0.1 -t $threads --preset meta-sensitive -o ${reads}/contigs_out --out-prefix $accession --min-contig-len 10 || megahit -r $R1 -m 0.1 -t $threads --preset meta-sensitive --out-prefix $accession -o ${reads}/contigs_out  --min-contig-len 10 # Run megahit
elif [[ $assembler == "megahit_sensitive_2" ]]; then
    [[ $paired == true ]] && megahit -1 $R1 -2 $R2 -m 0.1 -t $threads --preset meta-sensitive --prune-level 0 -o ${reads}/contigs_out --out-prefix $accession  --min-contig-len 10 || megahit -r $R1 -m 0.1 -t $threads --preset meta-sensitive --prune-level 0 --out-prefix $accession -o ${reads}/contigs_out --min-contig-len 10 # Run megahit
elif [[ $assembler == "megahit_sensitive_3" ]]; then
    echo "Running megahit with --min-count 1"
    [[ $paired == true ]] && megahit -1 $R1 -2 $R2 -m 0.1 -t $threads --min-count 1 --prune-level 0 --prune-depth 1 --max-tip-len 400 -o ${reads}/contigs_out --out-prefix $accession  --min-contig-len 10 || megahit -r $R1 -m 0.1 -t $threads  --min-count 1 --prune-level 0 --prune-depth 1 --out-prefix $accession -o ${reads}/contigs_out  --min-contig-len 10 # Run megahit
elif [[ $assembler == "megahit_sensitive_4" ]]; then
    echo "Running megahit with --min-count 1"
    [[ $paired == true ]] && megahit -1 $R1 -2 $R2 -m 0.1 -t $threads --min-count 1 --prune-level 0 --prune-depth 1 --max-tip-len 50 -o ${reads}/contigs_out --out-prefix $accession  --min-contig-len 10 || megahit -r $R1 -m 0.1 -t $threads  --min-count 1 --prune-level 0 --prune-depth 1 --out-prefix $accession -o ${reads}/contigs_out  --min-contig-len 10 # Run megahit

else
    if [[ $assembler == "spades" ]]; then
      [[ $paired == true ]] && spades.py -1 $R1 -2 $R2 -t $threads -o ${reads}/contigs_out || spades.py -s $R1 -t $threads -o ${reads}/contigs_out  # Run spades
    elif [[ $assembler == "metaspades" ]]; then
      [[ $paired == true ]] && spades.py -1 $R1 -2 $R2 -t $threads --meta -o ${reads}/contigs_out || spades.py -s $R1 -t $threads --meta -o ${reads}/contigs_out  # Run metaspades
    fi
    ${install_path}/ConvertContigs.sh -c ${reads}/contigs_out/contigs.fasta # Convert contigs to megahit format
    mv ${reads}/contigs_out/contigs.fasta ${reads}/contigs_out/${accession}.contigs.fa # Rename contigs file
fi