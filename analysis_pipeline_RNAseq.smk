import os, sys, glob

############# --- INPUT VARIABLES --- #############
### --- Change file paths here to your data --- ###
##### --- Full filepath starting with "/"" --- ####
### --- No slashes at the end of a filepath --- ###

fastqDirectory="/home/tmathieu/RNAseqAnalysis/fastq"
genome="/home/tmathieu/RNAseqAnalysis/input/Ath_genome.fa"
annotation="/home/tmathieu/RNAseqAnalysis/input/Ath_genes.gff"
threads=20

###################### --- PIPELINE--- #######################
### --- Don't change unless you know what you're doing --- ###
mainDirectory = os.getcwd()

SAMPLES, NUM, = glob_wildcards(os.path.join(mainDirectory, 'fastq', '{sample}_R{num}.fq'))

index = ''.join(genome.split("/")[-1].split(".")[:-1])
genomePreDir = genome.split("/")
genomeDir = '/'.join(genomePreDir[:len(genomePreDir)-1])
print(genomeDir)

mainDirectory = os.getcwd()

rule all:
	input: 
		geneCount = os.path.join(mainDirectory, 'summaryResults', 'merged_transcript_count.csv'),
		transcriptCount = os.path.join(mainDirectory, 'summaryResults', 'merged_gene_count.csv')

#rule uncompressGZ:
#	input:
#		os.path.join(mainDirectory, 'fastq', '{compressedfiles}q.gz')
#	output:
#		os.path.join(mainDirectory, 'fastq', '{compressedfiles}q')
#	shell: '''
#			gunzip -k {input}
#		'''

rule convertFASTQtoFQ:
	input:
		os.path.join(fastqDirectory, '{fastq}.fastq')
	threads: threads
	output:
		os.path.join(fastqDirectory, '{fastq}.fq')
	shell: '''
			mv {input} {output}
		'''

rule indexGenome:
	input:
		genome = genome
	threads: threads
	params:
		index = index,
		genomeDir = genomeDir
	output:
		expand(os.path.join(genomeDir, "{index}.{indexnum}.ht2"), index=index, indexnum = [1,2,3,4,5,6,7,8])
	shell: '''
			hisat2-build \
				-p {threads} \
				{input.genome} \
				{params.genomeDir}/{params.index}
		'''

rule alignFASTQs:
	threads: threads
	input:
		expand(os.path.join(genomeDir, "{index}.{indexnum}.ht2"), index=index, indexnum = [1,2,3,4,5,6,7,8]),
		r1 = os.path.join(fastqDirectory, '{sample}_R1.fq'),
		r2 = os.path.join(fastqDirectory, '{sample}_R2.fq')
	params:
		index = index,
		samDir = os.path.join(mainDirectory, 'sam'),
		genomeDir = genomeDir
	output:
		os.path.join(mainDirectory, 'sam', '{sample}.sam')
	shell: '''
			mkdir -p {params.samDir}

			hisat2 \
				-q \
				--dta \
				--no-unal \
				-p {threads} \
				--summary-file {params.samDir}/{wildcards.sample}.summary.txt \
				-x {genomeDir}/{params.index} \
				-1 {input.r1} \
				-2 {input.r2} \
				-S {output}
		'''

rule SAMtoBAM:
	input:
		os.path.join(mainDirectory, 'sam', '{sample}.sam')
	threads: threads
	output:
		os.path.join(mainDirectory, '{sample}.unsorted.bam')
	shell: '''
			samtools view \
				-@ {threads} \
				-Su \
				{input} \
				> {output}
		'''

rule sortBAMs:
	input:
		os.path.join(mainDirectory, '{sample}.unsorted.bam')
	threads: threads
	params:
		bamDir = os.path.join(mainDirectory, 'bam')
	output:
		os.path.join(mainDirectory, 'bam', '{sample}.bam')
	shell: '''
			mkdir -p {params.bamDir}

			samtools sort \
    			-@ {threads} \
    			{input} \
    			-o {output}

    		rm {input}
		'''

rule stringtieAssemble:
	input:
		bam = os.path.join(mainDirectory, 'bam', '{sample}.bam')
	threads: threads
	params:
		annotation = annotation,
		gtfDir = os.path.join(mainDirectory, 'gtf')
	output:
		gtf = os.path.join(mainDirectory, 'gtf', '{sample}', '{sample}_ST.gtf')
	shell: '''
			mkdir -p {params.gtfDir}/{wildcards.sample}

			stringtie \
				-B \
				-e \
				-p {threads} \
				-G {params.annotation} \
				-o {output.gtf} \
				{input.bam}
		'''

rule stringtieMerge:
	input:
		expand(os.path.join(mainDirectory, 'gtf', '{sample}', '{sample}_ST.gtf'), zip, sample = SAMPLES)
	threads: threads
	params:
		listGTF = os.path.join(mainDirectory, 'gtf', 'gtflist.txt'),
		annotation = annotation,
		gtfDir = os.path.join(mainDirectory, 'gtf'),
		summaryDir = os.path.join(mainDirectory, 'summaryResults')
	output:
		os.path.join(mainDirectory, 'summaryResults', 'merged_ST_1.gtf')
	shell: '''
			mkdir -p {params.summaryDir}

			ls {params.gtfDir}/*/*.gtf > {params.listGTF}

			stringtie \
				--merge \
				-p {threads}  \
				-G {params.annotation} \
				-o {output} \
				{params.listGTF}

			rm {params.listGTF}
		'''

rule stringtieAssembleWithMergedAnnotation:
	input:
		mergedAnnotation = os.path.join(mainDirectory, 'summaryResults', 'merged_ST_1.gtf'),
		bam = os.path.join(mainDirectory, 'bam', '{sample}.bam')
	threads: threads
	params:
		gtfMergeDir = os.path.join(mainDirectory, 'gtfMerge')
	output:
		gtf = os.path.join(mainDirectory, 'gtfMerge', '{sample}', '{sample}_ST_MERGE.gtf')
	shell: '''
		mkdir -p {params.gtfMergeDir}/{wildcards.sample}

		stringtie \
			-p {threads} \
			-e \
			-B \
			-G {input.mergedAnnotation} \
			-o {output.gtf} \
			{input.bam}
		'''

rule prepDE:
	input:
		expand(os.path.join(mainDirectory, 'gtfMerge', '{sample}', '{sample}_ST_MERGE.gtf'), zip, sample = SAMPLES)
	threads: threads
	params:
		gtfMergeDir = os.path.join(mainDirectory, 'gtfMerge'),
		prepDE = os.path.join(mainDirectory, 'script', 'prepDE.py3')
	output:
		geneCount = os.path.join(mainDirectory, 'summaryResults', 'merged_gene_count.csv'),
		transcriptCount = os.path.join(mainDirectory, 'summaryResults', 'merged_transcript_count.csv')
	shell: '''
			python3 {params.prepDE} \
				-i {params.gtfMergeDir} \
				-g {output.geneCount} \
				-t {output.transcriptCount} \
				-s AT
		'''
