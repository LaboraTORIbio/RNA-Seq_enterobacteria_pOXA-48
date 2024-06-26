#!/usr/bin/python3

import os, argparse, sys
import pandas as pd
#pd.set_option('display.max_rows', None)


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''Generates the TERM2GENE table for clusterProfiler''')

parser.add_argument('-g', '--golist', help="List of unique GO terms", type=str, dest="GOLIST", required=True)
parser.add_argument('-a', '--annot', help="Table of GO annotations retrived from UniProt", type=str, dest="ANNOTATION", required=True)
parser.add_argument('-m', '--mappings', help="Mapping of GeneIDs (1st column) and RefSeq IDs (2nd column)", type=str, dest="MAPPINGS", required=True)
parser.add_argument('-o', '--output', help="Path and name for the output file", type=str, dest="OUTPUT", required=True)

args = parser.parse_args()

GOList = args.GOLIST
Annotation = args.ANNOTATION	
Mappings = args.MAPPINGS
Output = args.OUTPUT


### Store the GO terms in a list

GOs = []
with open(GOList, 'r') as GOlist:
	for line in GOlist:
		line = line.rstrip()
		GOs.append(line)


### Create table with GO ID - RefSeq mappings

tableGO_RS = []
with open(Annotation, 'r') as annotation:
	for line in annotation:
		line = line.rstrip()
		linesplit = line.split("\t")
		for go in GOs:
			if go in line:
				tableGO_RS.append(go + "\t" + linesplit[0])

#for i in tableGO_RS:
#	print(i)


### Read table with Gene ID - RefSeq mappings

tableGen_RS = []
with open(Mappings, 'r') as mappings:
	for line in mappings:
		if not line.startswith("GeneID"):
			line = line.rstrip()
			tableGen_RS.append(line)

#for i in tableGen_RS:
#	print(i)


### Create TERM2GENE table (GO ID - Gene ID)

dfGO_RS = pd.DataFrame([x.split('\t') for x in tableGO_RS], columns = ["GOs", "RefSeq"])
dfGen_RS = pd.DataFrame([x.split('\t') for x in tableGen_RS], columns = ["GeneID", "RefSeq"])

#print(dfGO_RS)
#print(dfGen_RS)

df_merged = pd.merge(dfGO_RS, dfGen_RS, on="RefSeq")
df_merged = df_merged.drop('RefSeq', axis=1)
df_merged.sort_values('GOs')
#print(df_merged)


### Export table
df_merged.to_csv(Output, sep="\t", index=False, header=False)

