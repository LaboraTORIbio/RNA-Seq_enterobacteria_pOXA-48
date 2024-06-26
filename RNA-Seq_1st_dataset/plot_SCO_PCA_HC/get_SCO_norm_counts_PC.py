#! /usr/bin/python3

import glob
import re
import pandas as pd


vst_files = glob.glob("./vst_*.tsv")

strain_names = ["C325", "CF13", "H53", "J57", "K147", "MG1655", "EC10", "KPN04", "KPN07", "KPN10", "KPN15", "KPN18"]
cols_to_keep = ["GeneID", "Strain", "C325\..*", "CF13\..*", "H53\..*", "J57\..*", "K147\..*", "MG1655p\..*", "TC_EC10*", "TC_KPN04*", "TC_KPN07*", "TC_KPN10*", "TC_KPN15*", "TC_KPN18*"]
regex_patterns = [re.compile(pattern) for pattern in cols_to_keep]

df_list = []


for file in vst_files:
    df = pd.read_csv(file, sep='\t', header=0)
    
    for strain in strain_names:
        if any(df.columns[df.columns.str.contains(strain)]):
            df["Strain"] = strain
    
    matching_columns = [col for col in df.columns if any(pattern.match(col) for pattern in regex_patterns)]
    df_subset = df[matching_columns]
    df_copy = df_subset.copy()
    df_copy["GeneID"] = df_copy["GeneID"].str.replace("cds-", "")
    df_list.append(df_copy)

df_concat = pd.concat(df_list, axis=0)
df_concat.reset_index(drop=True, inplace=True)


df_SCO_mappings = pd.read_csv("SCO_mappings.tsv", sep='\t', header=0)

table = pd.merge(df_SCO_mappings, df_concat, on=["GeneID", "Strain"])
table = table.drop("GeneID", axis=1)
table = table.drop("Strain", axis=1)
collapsed_table = table.groupby("SCO").first().reset_index()

collapsed_table.to_csv('SCO_norm_counts_PC.tsv', sep="\t", index=False) 
