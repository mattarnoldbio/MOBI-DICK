# MetagenOmic BIoinformatics: Diamond Informed Classification and Krona plot (MOBI DICK) 

![mobi-dick logo](https://github.com/mattarnoldbio/MOBI-DICK/blob/master/mobi-dick.png)

The MOBI DICK pipeline processes sequencing reads, trimming, deduplicating, and optionally filtering reads from a provided host genome, before classifying using Diamond BLASTX and presenting as a Krona plot.

## Usage

- The scripts required to run the pipeline are in the `/Pipeline` directory. Only this directory is a required download for install. On Alpha2, the `Pipeline` directory is available at `/home3/2574106a/SRAmining/Pipeline` (no download required). *WARNING*: this version will be the bleeding edge beta; usually this will be fine but downloading your own copy might be more stable, especially if reproducibility is important for your results
- If you do not have `mamba` installed, install it (`bioconda` dpependencies will not resolve if conda is used to install the necessary packages). See the `miniconda` GitHub page for install istructions ([here](https://github.com/conda-forge/miniforge#mambaforge))
- run `mamba env create -f {PathToInstallOnYourSystem}/Pipeline/environment.yml`
- run `conda activate MOBI-DICK`. *You will need to do this at the start of every session where you want to run any pipeline scripts, as the pipeline has some python dependencies. If you forget, it WILL crash.*
- Runnning the pipeline with host filtering also requires a `bowtie2` indexed host genome. (Instructions below on how to do this.) Do this and take a note of the path to the directory containing the `.bt2` files.
  - run `IndexHostGenome.sh -g {path/to/host/genome} -o {location/to/save/indexed/genome}` to index your host genome. 
- Organise your raw data in a single directory, with a subdirectory for each sample containing the raw `fastq` files.
- To run the pipeline in metagenomics mode: 
  - `Pipeline.sh -d {directory containing subdirectories, each with one set of paired end reads, or one single-end read in} -g {path to basename of host genome in host genome directory, containing set of bt2 index files} -t {number of CPU cores to use. Recommended = 8}` (n.b. `-g` path must end with directory/basename (e.g. /home/Stuff/YourDirectory/YourGenome) where the .bt2 files are in YourDirectory and named YourGenome.1.bt2 etc.)
