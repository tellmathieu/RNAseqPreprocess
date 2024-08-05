# Preprocess for RNAseq DE Analysis
## Pipeline

Adapted from Patricia Baldrich's pipeline

Script for prepDE.py3 was retrieved from https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 with only minor changes made

This preprocess pipeline is specifically for paired end reads, so you must have 2 reads per sample.

## 1. Clone this repository to your working space

```
git clone https://github.com/tellmathieu/RNAseqAnalysis.git
```

## 2. Set up an environment to run this pipeline
- a) install mamba if you don't already have it - you'll have to go through the installation process
```
wget "https://github.com/conda-forge/miniforge/releases/latest/download/Miniforge3-$(uname)-$(uname -m).sh"
bash Miniforge3-$(uname)-$(uname -m).sh
```
- b) Initialize mamba and set up an environment - this can take 5-10 minutes - hisat2 and stringtie need to be these specific versions for recreating our results from original pipeline
```
source ~/.bashrc #or you can log out and then log back in
mamba create -c conda-forge -c bioconda -y -n preProcessRNAseq snakemake=8.15.2 hisat2=2.2.1 samtools=1.20 stringtie=2.1.4 
mamba activate preProcessRNAseq
```

## 3. Give your fastq files a consistent naming convention.

- `{sample_name}_R1.fastq` or `{sample_name}_R1.fq`
- `{sample_name}_R2.fastq` or `{sample_name}_R2.fq`

The `{sample_name}` is what you choose (don't put any spaces in this name). The files must have the suffixes, including the read number (1 or 2).

An example of how my files were named is:

- `Col0_AWF_NT_R4_R1.fastq`
- `Col0_AWF_NT_R4_R2.fastq`


## 4. Change filepaths to your filepaths in the analysis_pipeline_RNAseq.smk

```
cd RNAseqAnalysis
nano preprocess_RNAseq_DE.smk
```
In the nano editor (or whatever your preferred editor is), change the following 8 variables:
```
############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###

fastqDirectory="/home/RNAseqAnalysis/fastq"
genome="/home/RNAseqAnalysis/input/Ath_genome.fa"
annotation="/home/RNAseqAnalysis/input/Ath_genes.gff"
```
Make sure you have read/write privileges in these files.

The `fastqDirectory` variable is the folder where you have all your fastq files (with the naming convention I specified).

The `genome` variable is the filepath to the `.fasta` or `.fa` file you are using to map the RNAseq data against.

The `annotation` variable is the filepath to the `.gff` file you are using to map the RNAseq data against. You usually get this from where you get the genome file from. 


## 7. Run the pipeline.

Now, you're ready to run the pipeline. Most of the time is taken with mapping the RNAseq (fastq) files to the genome and annotation. With 12 samples and 20 threads, it took a little more than a day for me. I usually run this using the `screen` software (so you can detach but the command continues to run).

The `-j 20` specification is how many threads I told the computer to allocate to this job. I made it the same as what I specified in the file itself, so each "job" the computer runs gets the full 20 threads.
```
snakemake --snakefile preprocess_RNAseq_DE.smk -j 20 -p
```
