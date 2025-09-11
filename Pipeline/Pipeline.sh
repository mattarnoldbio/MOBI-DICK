#!/usr/bin/bash


################################################################################
# Help                                                                         #
################################################################################
Help()
{

echo MOBI-DICK PIPELINE
echo
echo Matt Arnold 2023
echo
echo This script is used to run the whole pipeline
echo
echo ARGUMENTS
echo -e '\t '          -h --help '\t\t\t'  display this help and exit
echo
echo required
echo -e Input:'\t '    -d --data_dir '\t\t'  path to directory containing raw data
echo
echo recommended
echo -e '\t '          -b --database '\t\t'  path to protein database to be used "(default = /db/diamond/nr.dmnd)"
echo -e '\t '          -u --nuc_database '\t\t'  path to database to be used for blastn "(default = /db/blast_v5/nt)"
echo -e '\t '          -k --krona_tools_db '\t\t'  path to krona tools taxonomy database "(default = /db/kronatools/taxonomy)"
echo -e '\t '          -t --threads '\t\t\t'  number of cpu threads to use "(Default = 8)"
echo -e '\t '          -g --genome '\t\t\t'    path to host reference genome "(or multiple concatenated genomes)"
echo -e '\t \t\t\t\t '                         "(GENOME MUST BE INDEXED WITH BOWTIE2, see Pipeline/IndexHostGenome.sh)"
echo -e '\t \t\t\t\t '                         "path must end with directory/basename (e.g. home/Stuff/YourDirectory/YourGenome)"
echo -e '\t \t\t\t\t '                         "where the .bt2 files are in YourDirectory and named YourGenome.1.bt2 etc."
#echo -e '\t '          -o --out_file '\t\t'  output file name
echo
echo SRA mining mode
echo -e '\t '          -r --host_ref_genome_dir '\t' path to directory containing host reference genomes
echo -e '\t '          -m --metadata '\t\t'  path to metadata file containing species names for each accession number
echo -e '\t '          -s --species_column '\t\t'  column number of species in metadata file
echo
echo -e optional i.e. not recommended
echo -e '\t '          -n --no_host_filtering '\t'  skip host filtering "(Default = false)"
echo -e '\t '          -f --score_filter '\t\t'  score filter for krona plots "(default = -10)"
echo -e '\t '          -a --assembler '\t\t'  assembler to use "(megahit, spades or metaspades default = megahit)"
echo -e '\t'

}

################################################################################
################################################################################

# PARSE COMMAND LINE ARGUMENTS
while [[ "$#" -gt 0 ]]
  do
    case $1 in
      -d|--data_dir) data_dir="$2"; shift;;
      -g|--genome) genome="$2"; shift;;
      -m|--metadata) metadata="$2"; shift;;
      -s|--species_column) species_column="$2"; shift;;
      -r|--host_ref_genome_dir) host_ref_genome_dir="$2"; shift;;
      -t|--threads) threads="$2"; shift;;
      -n|--no_host_filtering) no_host_filtering=false; shift;;
      -b|--database) database="$2"; shift;;
      -u|--nuc_databse) nuc_database="$2"; shift;;
      -k|--krona_tools_db) krona_tools_db="$2"; shift;;
      -f|--score_filter) score_filter="$2"; shift;;
      -a|--assembler) assembler="$2"; shift;;
      -h|--help) Help; exit 1;;
    esac
    shift
done



if [[ ! $(echo $CONDA_DEFAULT_ENV) == "MOBI-DICK" ]]; then
  echo "###########################################################################################################"
  echo "WARNING!"
  echo"Your current conda env is" $(echo $CONDA_DEFAULT_ENV) "If the pipeline crashes, it might be because of this!"
  echo "###########################################################################################################"

fi

# SET DEFAULTS
install_path=$(dirname -- "$0")/ # Get path to install directory

[[ -z $threads ]] && threads=8 # Default number of threads is 10
[[ -z $database ]] && database=/db/diamond/nr.dmnd # Default database is nr
[[ -z $nuc_database ]] && nuc_database=/db/blast_v5/nt # Default database is nt
[[ -z $krona_tools_db ]] && krona_tools_db=/db/kronatools/taxonomy # Default krona tools database is taxonomy
[[ -z $score_filter ]] && score_filter=-10 # Default score filter is -10
[[ -z $assembler ]] && assembler=megahit # Default assembler is megahit



log_file=${data_dir}/mobi_dick_$(date +"%d_%m_%y_%H_%M_%S").log
echo Output logged to $log_file
touch $log_file
echo Processing files in $data_dir 

exec 3>&1 1> $log_file 2>&1

trap "date -Is" DEBUG

# RUN PIPELINE


for directory in $(ls -d $data_dir/*/); do # Loop through all directories in data_dir

  if [[ -d ${directory}/contigs_out ]]; then
    contigs_out=$(ls ${directory}/contigs_out | grep contigs.fa)
    contig_file_length=$(wc -l $directory/contigs_out/$contigs_out | awk '{print $1}') # Get number of lines in contig file
    [[ contig_file_length -ne 0 ]] && echo contigs found for ${directory}, skipping processing    && continue # Skip if contigs found
  fi

  touch ${directory}/log.txt # Create log file
  echo "Processing $directory"  ${directory}/log.txt # Print directory name to log file
  echo Trimming reads  
  ## TRIM READS
  ${install_path}ReadTrimming.sh -r $directory -o $directory -t $threads   # Run read trimming
  echo Deduplicating reads  
  ## DEDUPLICATE READS
  ${install_path}ReadDeduplication.sh -r $directory -o $directory -t $threads  # Run read deduplication
  echo Attempting host filtering  
  ## FILTER HOST READS  
  [[ -z $genome ]] && [[ ! -z $metadata ]] && ${install_path}HostFiltering.sh -r $directory  -h $host_ref_genome_dir -o $directory -t $threads -m $metadata -s $species_column -i $install_path   # Run host filtering if a metadata file was provided
  [[ -z $metadata ]]  && [[ ! -z $genome ]] && ${install_path}HostFiltering.sh -r $directory  -o $directory -t $threads -g $genome -i $install_path  # Run host filtering if a host genome was provided
  [[ -z $metadata ]] && [[ -z $genome ]] && echo No host genome or metadata file specified, skipping host filtering  # If neither a metadata file or host genome was provided, skip host filtering
  echo Assembling reads using assembly preset $assembler 
  ## ASSEMBLE READS
  ${install_path}Assembly.sh -r $directory -o $directory -t $threads -i $install_path -a $assembler --no_host_filtering $no_host_filtering  # Run assembly
  
  ## COMPRESS FASTQ FILES
  ${install_path}CompressFASTQs.sh -r $directory -t $threads # Compress FASTQ files

  ## CHECK PROGRESS
  ${install_path}CheckProgress.sh -r $data_dir  # Check progress
done

echo All reads processed. Searching database for hits to viral genomes 

## RUN DIAMOND
${install_path}RunDiamond.sh -d $data_dir -t $threads -b $database -k $krona_tools_db  # Run diamond on all contigs (for perfomrance reasons, this is done on a file of all contigs concatenated together)

## POST PROCESS DIAMOND RESULTS
${install_path}PostProcessKrona.sh -d $data_dir -i ${install_path} -w diamond -s $score_filter   # Post process diamond results

## RUN BLAST
${install_path}RunBLAST.sh -d $data_dir -b $nuc_database -t $threads -k $krona_tools_db  # Run blast on all contigs (for perfomrance reasons, this is done on a file of all contigs concatenated together)

## POST PROCESS BLAST RESULTS
#${install_path}PostProcessKrona.sh -d $data_dir -i ${install_path} -w blastn -s $score_filter # Post process diamond results

echo Pipeline complete 
echo Pipeline complete >&3