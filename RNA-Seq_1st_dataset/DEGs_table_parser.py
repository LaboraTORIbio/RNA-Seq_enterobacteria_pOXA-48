#! /usr/bin/python3

import openpyxl
import pandas as pd
import os, argparse, sys
from functools import reduce


parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description='''DEGs TABLE PARSER''')

parser.add_argument('-d', nargs="+", help="Directories to parse", type=str, dest="DIRECTORIES", required=True)
parser.add_argument('-n', help="Base name for the output file", type=str, dest="NAME", required=True)

args = parser.parse_args()

DirNames = args.DIRECTORIES
DirNames = [i.replace("/", "") for i in DirNames]
DirNamesSheet = [i.replace("RNAseq_", "") for i in DirNames]
DirNamesSheet1 = [i+"_chr" for i in DirNamesSheet]
DirNamesSheet2 = [i+"_plas" for i in DirNamesSheet]
DirNamesSheet2.remove("MG1655_plas")
OutName = args.NAME
OutName1 = OutName+".xlsx"


fulldfListChr = []
dfListChr = []
fulldfListPlas = []
dfListPlas = []

for directory in DirNames:
	if os.path.isdir(directory):
		files = os.listdir(directory)
		for file in files:
			if "DE_results_filtered_chromosome_padj_GOannot" in file:
				filepath = directory+"/"+file
				fulldf = []
				df = []
				with open(filepath, 'r') as file:
					for line in file:
						line = line.rstrip("\n")
						line = line.split("\t")
						line = [i.replace('"', '') for i in line]
						fulldf.append(line)
						if "GeneID" in line or (("CDS" in line[-4] or "CDS" in line[-5]) and "hypothetical protein" not in line):
							df.append(line)
				
				fulldf = pd.DataFrame(fulldf)
				fulldf.columns = fulldf.iloc[0]
				fulldf = fulldf[1:]
				fulldfListChr.append(fulldf)
				
				df = pd.DataFrame(df)
				df.columns = df.iloc[0]
				df = df[1:]
				
				dfsub = df[["Type", "RefSeq", "Gene", "Product", "log2FoldChange"]]
				dfsub = dfsub.astype({"log2FoldChange":float})
				dfsub["log2FoldChange"] = dfsub["log2FoldChange"].round(decimals=2)
				strain = directory.replace("RNAseq_", "")
				strain = strain.replace("/", "")
				dfsub = dfsub.rename({"log2FoldChange": strain}, axis=1)
				dfListChr.append(dfsub)
			
			
			elif "DE_results_filtered_plasmids_padj_GOannot" in file:
				filepath = directory+"/"+file
				fulldf = []
				df = []
				with open(filepath, 'r') as file:
					for line in file:
						line = line.rstrip("\n")
						line = line.split("\t")
						line = [i.replace('"', '') for i in line]
						fulldf.append(line)
						if "GeneID" in line or (("CDS" in line[-4] or "CDS" in line[-5]) and "hypothetical protein" not in line):
							df.append(line)
				
				fulldf = pd.DataFrame(fulldf)
				fulldf.columns = fulldf.iloc[0]
				fulldf = fulldf[1:]
				fulldfListPlas.append(fulldf)
				
				df = pd.DataFrame(df)
				df.columns = df.iloc[0]
				df = df[1:]
				#print(df)
				
				dfsub = df[["Type", "RefSeq", "Gene", "Product", "log2FoldChange"]]
				dfsub = dfsub.astype({"log2FoldChange":float})
				dfsub["log2FoldChange"] = dfsub["log2FoldChange"].round(decimals=2)
				strain = directory.replace("RNAseq_", "")
				strain = strain.replace("/", "")
				dfsub = dfsub.rename({"log2FoldChange": strain}, axis=1)
				dfListPlas.append(dfsub)
				#print(dfsub)
	
	else:
		print(directory + " is an invalid directory name")


dfChr_merged = reduce(lambda left,right: pd.merge(left,right,on=["Type", "RefSeq", "Gene", "Product"], how='outer'), dfListChr).fillna("0")
dfPlas_merged = reduce(lambda left,right: pd.merge(left,right,on=["Type", "RefSeq", "Gene", "Product"], how='outer'), dfListPlas).fillna("0")



writer = pd.ExcelWriter(OutName1, engine='xlsxwriter')

dfChr_merged.to_excel(writer, sheet_name='All_chromosome_DEGs', index=False)
dfPlas_merged.to_excel(writer, sheet_name='All_plasmids_DEGs', index=False)

for df, n in zip(fulldfListChr, DirNamesSheet1):
	df.to_excel(writer, index=False, sheet_name=n)
for df, n in zip(fulldfListPlas, DirNamesSheet2):
	df.to_excel(writer, index=False, sheet_name=n)
	
writer.save()

