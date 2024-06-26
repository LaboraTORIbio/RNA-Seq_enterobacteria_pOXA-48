#!/usr/bin/python3

import os, re, sys
import pandas as pd
import numpy as np


# Read table

merged_df = pd.read_csv("heatmaps/DEGs4heatmap_full.tsv", sep='\t')


# Table of genes common to all strains with at least one significant log2FC:

# Get the count of unique factors (strains) for each Gene identifier
Gene_counts = merged_df.groupby('Gene')['Strain'].nunique()
# Find the Gene identifiers that have counts equal to 12 (present in all strains)
Genes_with_all_factors = Gene_counts[Gene_counts == 12].index
# Filter the DataFrame to keep only the rows with the Gene identifiers having all factors
filtered_df = merged_df[merged_df['Gene'].isin(Genes_with_all_factors)]

# Create a copy of the DataFrame
df_copy = filtered_df.copy()
# Convert 'padj' column to numeric type
df_copy['padj'] = pd.to_numeric(df_copy['padj'], errors='coerce')
# Group by Gene and check if any p-value is significant
significant_Genes = df_copy.groupby('Gene')['padj'].apply(lambda x: any(x < 0.05))
# Get the Gene identifiers with at least one significant p-value
Genes_with_significance = significant_Genes[significant_Genes].index
# Filter the data frame to keep only the rows with significant Gene identifiers
filtered_df = df_copy[df_copy['Gene'].isin(Genes_with_significance)]



# Add column of parental grouping

parental = []
with open('GSEA/enriched_GOs_parental_BP.tsv', 'r') as file:
	for line in file:
		line = line.rstrip("\n")
		if line.startswith("GO"):
			line = line.split("\t")
			parental.append(line[0] + "\t" + line[3])

parental = pd.DataFrame([x.split('\t') for x in parental], columns = ["GO_BP", "Parental"])
final_df = filtered_df.merge(parental, on='GO_BP', how='inner')



# Save to tsv
final_df.to_csv("heatmaps/DEGs4heatmap.tsv", sep="\t", index=False, header=True)

