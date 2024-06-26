#! /usr/bin/python3

import glob
import os
import re
import pandas as pd
pd.set_option('display.max_rows', None)


rnaseq_dirs = ["../RNAseq_C325", "../RNAseq_CF13", "../RNAseq_H53", "../RNAseq_J57", "../RNAseq_K147", "../RNAseq_MG1655", "../RNAseq_EC10", "../RNAseq_KPN04", "../RNAseq_KPN07", "../RNAseq_KPN10", "../RNAseq_KPN15", "../RNAseq_KPN18"]
strain_names = ["C325", "CF13", "H53", "J57", "K147", "MG1655", "EC10", "KPN04", "KPN07", "KPN10", "KPN15", "KPN18"]

DEGs_files = glob.glob("../RNAseq_*/DE_results_raw_chromosome.tsv")

df_list = []

for file in DEGs_files:
    if os.path.dirname(file) in rnaseq_dirs:
        df = pd.read_csv(file, sep='\t', header=0)
        
        for strain in strain_names:
            if any(df.columns[df.columns.str.contains(strain)]):
                df["Strain"] = strain
        
        df_subset = df[["Strain", "GeneID", "log2FoldChange", "padj", "Gene", "Product"]]
        df_copy = df_subset.copy()
        df_copy["GeneID"] = df_copy["GeneID"].str.replace("cds-", "").str.replace("rna-", "")
        df_list.append(df_copy)

df_concat = pd.concat(df_list, axis=0)
df_concat.reset_index(drop=True, inplace=True)

df_SCO_mappings = pd.read_csv("SCO_mappings.tsv", sep='\t', header=0)

table = pd.merge(df_SCO_mappings, df_concat, on=["GeneID", "Strain"])
table['padj'] = table['padj'].fillna(False)
table['padj'] = table['padj'] < 0.05
table.to_csv('SCO_log2FC_full.tsv', sep="\t", index=False)

table_subset = table[["Strain", "SCO", "log2FoldChange", "padj"]]
table_subset.to_csv('SCO_log2FC.tsv', sep="\t", index=False)
