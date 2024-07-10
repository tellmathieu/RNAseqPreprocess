import pandas

#******************************
#**** Not used in snakemake ***
#**** Used for debugging ******


#origTable = pandas.read_csv('/Users/tell/Desktop/test_collapse.csv')
#newTable = pandas.read_csv('/Users/tell/Desktop/test_collapse_out.csv')
origTable = pandas.read_csv('/home/tmathieu/RNAseqAnalysis/2genecount.csv')
newTable = pandas.read_csv('/home/tmathieu/RNAseqAnalysis/summaryResults/merged_gene_count.csv')

#origTable['gene'] = origTable['gene_id'].str.split('.').str[0]

#origTable = origTable.drop['gene_id']

#origGenes = origTable['gene_id']

#allGenes = origGenes.append(newTable['gene_id'])

allGenes = pandas.concat([origTable, newTable])

allGenes['Counts'] = allGenes['gene_id'].map(allGenes['gene_id'].value_counts())

#allGenesCounts = allGenes.value_counts()

uniqueGenes = allGenes[allGenes['Counts'] <= 1]

#print(allGenes.value_counts())

print(uniqueGenes)
#print(clippedGeneCounts)
