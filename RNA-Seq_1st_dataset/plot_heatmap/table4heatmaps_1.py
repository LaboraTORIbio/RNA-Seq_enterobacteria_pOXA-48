#!/usr/bin/python3

import os, re, sys
import pandas as pd
import numpy as np



# GENERATE A DATA FRAME OF STRAIN NAME, REFSEQ ID, GENE NAME, PRODUCT, LOG2FOLDCHANGE AND PADJ

# Directories with differential expression analysis data
DirNames = ["RNAseq_C325", "RNAseq_CF13", "RNAseq_H53", "RNAseq_J57", "RNAseq_K147", "RNAseq_EC10", "RNAseq_KPN04", "RNAseq_KPN07", "RNAseq_KPN10", "RNAseq_KPN15", "RNAseq_KPN18", "RNAseq_MG1655"]

df = []
for directory in DirNames:
	strain = directory.replace("RNAseq_", "")
	if os.path.isdir(directory):
		files = os.listdir(directory)
		for file in files:
			if "DE_results_raw_chromosome" in file:
				filepath = directory+"/"+file
				with open(filepath, 'r') as file:
					for line in file:
						if "CDS" in line:
							line = line.rstrip("\n")
							line = line.replace("\"", "")
							line = line.split("\t")
							refseq = re.sub(r'\.\d', '', line[-3])
							df.append(strain + "\t" + refseq + "\t" + line[-2] + "\t" + line[-1] + "\t" + line[2] + "\t" + line[5])
							


# GENERATE A DATA FRAME OF REFSEQ ID AND ENRICHED BIOLOGICAL PROCESSES

annotfiles = ["GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv", "GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv"]

enrichedbp = ["GO:0007155", "GO:0007049", "GO:0071973", "GO:0071978", "GO:0051301", "GO:0044780", "GO:0044781", "GO:0071555", "GO:0009252", "GO:0017004", "GO:0016226", "GO:0008360", "GO:0009437", "GO:0009061", "GO:0016310", "GO:0006099", "GO:0009103", "GO:0000105", "GO:0019439", "GO:0042128", "GO:0006412", "GO:0005975", "GO:0006526", "GO:0009097", "GO:0009228", "GO:0006189", "GO:0022904", "GO:0006096", "GO:0009245", "GO:0006457", "GO:0006355", "GO:0006935", "GO:0007165", "GO:0046677", "GO:0006281", "GO:0009306", "GO:0006865", "GO:0015628", "GO:0055085", "GO:0009401"]

annot = []
for files in annotfiles:
	with open(files, 'r') as file:
		for line in file:
			if any(bp in line for bp in enrichedbp):
				line = line.rstrip("\n")
				line = line.split("\t")
				for bp in enrichedbp:
					if bp in line[2]:
						go_term = line[2].split("; ")
						for term in go_term:
							if bp in term:
								annot.append(line[0] + "\t" + term.split("[")[1].split("]")[0] + "\t" + term.split("[")[0])
								break
							break



# MERGE DATA FRAMES

df = pd.DataFrame([x.split('\t') for x in df], columns = ["Strain", "RefSeq", "Gene", "Product", "log2FC", "padj"])
annot = pd.DataFrame([x.split('\t') for x in annot], columns = ["RefSeq", "GO_BP", "GO_Description"])
annot = annot.drop_duplicates()

merged_df = df.merge(annot, on='RefSeq', how='inner')
merged_df = merged_df.sort_values("RefSeq")


# Save table
merged_df.to_csv("heatmaps/DEGs4heatmap_full.tsv", sep="\t", index=False, header=True)


