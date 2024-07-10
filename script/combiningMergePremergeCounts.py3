import os, sys, pandas

def combineCountFiles():
	premergefile, mergedfile, combinedfile = getArgs()
	premerge = pandas.read_csv(premergefile, index_col=[0])
	merged = pandas.read_csv(mergedfile, index_col=[0])
	merged = merged.add_suffix('_merged')
	combined = pandas.concat([premerge, merged], axis=1, join="outer")
	combined.to_csv(combinedfile)
	return


def getArgs():
	premergefile = sys.argv[1]
	mergedfile = sys.argv[2]
	combinedfile = sys.argv[3]
	return premergefile, mergedfile, combinedfile

combineCountFiles()
