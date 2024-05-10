# RNAseqAnalysis
RNAseq Analysis Pipeline
(Adapted from Patricia Baldrich's pipeline)

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
- b) Initialize mamba and set up an environment - this can take 5-10 minutes
```
source ~/.bashrc #or you can log out and then log back in
mamba create -c conda-forge -c bioconda -y -n RNAseq snakemake hisat2=2.2.1 samtools stringtie=2.2.3 bioconductor-deseq2 bioconductor-edger
mamba activate RNAseq
```

## 3. Change filepaths to your filepaths in the analysis_pipeline_RNAseq.smk
```
cd RNAseqAnalysis
vim analysis_pipeline_RNAseq.smk
```
In the vim editor (or whatever your preferred editor is), change the following 4 variables:
```
############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###

fastqDirectory="/home/tmathieu/RNAseq_analysis_pipeline/fastq"
genome="/home/tmathieu/RNAseq_analysis_pipeline/input/Ath_genome.fa"
annotation="/home/tmathieu/RNAseq_analysis_pipeline/input/Ath_genes.gff"
threads=20
```
Make sure you have read/write privileges in these files.

## 4. Run the pipeline.
```
snakemake --snakefile analysis_pipeline_RNAseq.smk -j 20 -p
```
