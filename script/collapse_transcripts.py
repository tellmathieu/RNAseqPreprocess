import pandas

#******************************
#**** Not used in snakemake ***
#**** Used for debugging ******

#geneCountTable = pandas.read_csv('/Users/tell/Desktop/test_collapse.csv')
geneCountTable = pandas.read_csv('/home/tmathieu/RNAseqAnalysis/1genecount.csv')

geneCountTable['gene'] = geneCountTable['gene_id'].str.split('.').str[0]

duplicates = geneCountTable.duplicated(subset=['gene'])

grouped = geneCountTable.groupby('gene').agg({"sum"})

grouped2 = grouped.drop(columns=['gene_id'])
#print(geneCountTable['gene'][4][0])

grouped2.to_csv('/home/tmathieu/RNAseqAnalysis/2genecount.csv')
