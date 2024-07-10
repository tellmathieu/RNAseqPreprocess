import pandas
from gtfparse import read_gtf
from numpy import nan

#******************************
#**** Not used in snakemake ***
#**** Used for debugging ******

samples = ['Col0_AWF_NT_R4','Col0_AWF_NT_R5','Col0_AWF_NT_R6','Col0_CL_R1','Col0_CL_R2','Col0_CL_R3','Col0_CL_ribo_R1','Col0_CL_ribo_R2','Col0_CL_ribo_R3','Col0_LSW_R1','Col0_LSW_R2','Col0_LSW_R3']

allSamples = pandas.DataFrame(columns=['transcript_id'])
mismatchAllSamples = pandas.DataFrame(columns=['transcript_id'])
mismatchAllSamplesSig = pandas.DataFrame(columns=['transcript_id'])

for sample in samples:
	origGtf = read_gtf('/home/tmathieu/RNAseqAnalysis/gtf/%s/%s_ST.gtf' % (sample,sample),column_converters={"FPKM": float})
	newGtf = read_gtf('/home/tmathieu/RNAseqAnalysis/gtfMerge/%s/%s_ST_MERGE.gtf' % (sample,sample),column_converters={"FPKM": float})

	#joinedGtf = pandas.merge(origGtf, newGtf, how="outer", on="transcript_id")
	#origGenes = origGtf[(origGtf['feature']=='transcript')]
	#newGenes = newGtf[newGtf["feature"] == "gene"]

	origGenes = pandas.DataFrame(origGtf, columns=origGtf.columns)
	origGenes = origGenes[origGenes['feature'] == 'transcript']
	origGenes = origGenes[['gene_id','transcript_id','FPKM']]

	newGenes = pandas.DataFrame(newGtf, columns=newGtf.columns)
	newGenes = newGenes[newGenes['feature'] == 'transcript']
	newGenes = newGenes[['gene_id','transcript_id','FPKM']]

	joinedGenes = pandas.merge(origGenes, newGenes, how="outer", on="transcript_id", suffixes=('_orig%s' % sample, '_merged%s' % sample))
	query = 'FPKM_orig%s != FPKM_merged%s' % (sample,sample)
	mismatchGenes = joinedGenes.query(query)
	origFPKM = 'FPKM_orig' + sample
	mergedFPKM = 'FPKM_merged' + sample
	diffFPKM = 'FPKM_diff' + sample
	mismatchGenes[diffFPKM] = mismatchGenes[origFPKM] - mismatchGenes[mergedFPKM]
	mismatchGenesSig = mismatchGenes[mismatchGenes[diffFPKM] > 0.01]
	mismatchGenesSigMinus = mismatchGenes[mismatchGenes[diffFPKM] < -0.01]
	mismatchGenesSig = pandas.merge(mismatchGenesSig, mismatchGenesSigMinus, how="outer", on="transcript_id")
	allSamples = pandas.merge(allSamples, joinedGenes, how="outer", on="transcript_id")
	mismatchAllSamples = pandas.merge(mismatchAllSamples, mismatchGenes, how="outer", on="transcript_id")
	mismatchAllSamplesSig = pandas.merge(mismatchAllSamplesSig, mismatchGenesSig, how="outer", on="transcript_id")


allSamples.to_csv('FPKM.csv')
mismatchAllSamples.to_csv('mismatchFPKM.csv')
mismatchAllSamplesSig.to_csv('sigMismatchFPKM0.01.csv')

mismatchAllSamples = mismatchAllSamples.replace(' ',nan)

summarynull = mismatchAllSamples.isnull().sum()
summaryna = mismatchAllSamples.isna().sum()

summaryna.to_csv('na.txt')
summarynull.to_csv('null.txt')
