# RNAseqAnalysis
## RNAseq Analysis Pipeline

Adapted from Patricia Baldrich's pipeline

Script for prepDE.py3 was retrieved from https://ccb.jhu.edu/software/stringtie/dl/prepDE.py3 with only minor changes made

This analysis is specifically for paired end reads, so you must have 2 reads per sample.

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
mamba create -c conda-forge -c bioconda -y -n RNAseq snakemake=8.15.2 hisat2=2.2.1 samtools=1.20 stringtie=2.1.4 r-upsetr=1.4.0 bioconductor-deseq2=1.42.0 bioconductor-apeglm=1.24.0 r-dplyr=1.1.4 r-ggplot2=3.5.1 r-ggrepel=0.9.5 r-rcolorbrewer=1.1_3 r-rlang=1.1.4 r-purrr=1.0.1 r-gplots=3.1.3.1 r-devtools=2.4.5
mamba activate RNAseq
```

## 3. Give your fastq files a consistent naming convention.

### `{sample_name}_R1.fastq` or `{sample_name}_R1.fq`
### `{sample_name}_R2.fastq` or `{sample_name}_R2.fq`

The `{sample_name}` is what you choose (don't put any spaces in this name). The files must have the suffixes, including the read number (1 or 2).

An example of how my files were named is:

## `Col0_AWF_NT_R4_R1.fastq`
## `Col0_AWF_NT_R4_R2.fastq`

## 4. Create metadata file (csv)

This is your description of all the different RNAseq samples you have, and the axes you want to compare when you run the differential expression analysis.

The two required columns are `sample` and `DE`. You can add more, but this program doesn't read them. The `sample` column must be the first column.

The `sample` column must match `{sample_name}`. The `DE` column should have your treatment names.
```
sample,DE
sample1,treatment1
sample2,treatment1
sample3,treatment1
sample4,treatment2
sample5,treatment2
sample6,treatment2
sample7,control
sample8,control
sample9,control
```

Here is an example of what my file looked like:
```
sample,treatment,replicate,DE
Col0_AWF_NT_R4,Col0_AWF_NT,R4,Col0_AWF_NT
Col0_AWF_NT_R5,Col0_AWF_NT,R5,Col0_AWF_NT
Col0_AWF_NT_R6,Col0_AWF_NT,R6,Col0_AWF_NT
Col0_CL_R1,Col0_CL,R1,Col0_CL
Col0_CL_R2,Col0_CL,R2,Col0_CL
Col0_CL_R3,Col0_CL,R3,Col0_CL
Col0_CL_ribo_R1,Col0_CL_ribo,R1,Col0_CL_ribo
Col0_CL_ribo_R2,Col0_CL_ribo,R2,Col0_CL_ribo
Col0_CL_ribo_R3,Col0_CL_ribo,R3,Col0_CL_ribo
```

## 5. Create treatment vs comparisons (treatmentComparisons) file (csv)

The `num` column is a sequentially numbered column.

The `treatment` column are the treatments that you treated your samples with. These must be exactly the same as the `DE` column in step 4.

The `control` column are what you consider controls in your experiment. These must be exactly the same as the `DE` column in step 4.

```
num,treatment,control
1,treatment1,control
2,treatment2,control
```

Here is an example of my file that I used:
```
num,treatment,control
1,Col0_AWF_NT,Col0_CL
2,Col0_LSW,Col0_CL
3,Col0_LSW,Col0_AWF_NT
4,Col0_CL_ribo,Col0_CL
```

## 6. Change filepaths to your filepaths in the analysis_pipeline_RNAseq.smk

```
cd RNAseqAnalysis
nano analysis_pipeline_RNAseq.smk
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
metadata="/home/RNAseqAnalysis/input/colData.csv"
treatmentComparisons="/home/RNAseqAnalysis/input/treatVsControl.csv"
alpha=0.05
threads=20
ctrl="Col0_CL"
```
Make sure you have read/write privileges in these files.

The `fastqDirectory` variable is the folder where you have all your fastq files (with the naming convention I specified).

The `genome` variable is the filepath to the `.fasta` or `.fa` file you are using to map the RNAseq data against.

The `annotation` variable is the filepath to the `.gff` file you are using to map the RNAseq data against. You usually get this from where you get the genome file from. 

The `metadata` variable is the filepath to the csv file you created in step 4.

The `treatmentComparisons` variable is the filepath to the csv file you created in step 5.

The `alpha` variable is your cutoff p-value for significance.

The `threads` variable is for the computer to know how many threads to allocate each step. More threads usually means its quicker, but computer systems have limits, so that's why I have you specify it.

The `ctrl` variable is the main control of your experiment. It must match exactly one of the options in the `DE` column in step 4 and one of the options in the `control` column in step 5. The reason we specify it in both step 5 and this step is because when you run differential expression analysis you need the software (DEseq2) to know which is the control, but I also wanted to allow for more comparisons to be specified if you wanted, hence step 5.


## 7. Run the pipeline.

Now, you're ready to run the pipeline. Most of the time is taken with mapping the RNAseq (fastq) files to the genome and annotation. With 12 samples and 20 threads, it took a little more than a day for me. I usually run this using the `screen` software (so you can detach but the command continues to run).

The `-j 20` specification is how many threads I told the computer to allocate to this job. I made it the same as what I specified in the file itself, so each "job" the computer runs gets the full 20 threads.
```
snakemake --snakefile analysis_pipeline_RNAseq.smk -j 20 -p
```
