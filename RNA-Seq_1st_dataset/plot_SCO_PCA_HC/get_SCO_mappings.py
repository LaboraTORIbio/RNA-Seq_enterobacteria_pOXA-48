#! /usr/bin/python3

import os
import re
import pandas as pd


def list_files_with_path(directory):
    files_with_path = []
    for filename in os.listdir(directory):
        file_path = os.path.join(directory, filename)
        if os.path.isfile(file_path):
            files_with_path.append(file_path)
    return files_with_path


def parse_SCO(directory):
    fa_files = list_files_with_path(directory)    
    strain_names = ["C325", "CF13", "H53", "J57", "K147", "MG1655", "EC10", "KPN04", "KPN07", "KPN10", "KPN15", "KPN18"]    
    df_list = []
    
    for file in fa_files:
        orthogroup = file.replace(directory, "").replace(".fa", "")
        with open(file, 'r') as f:
            file_content = f.read()
            geneids = re.findall(r'pgaptmp_\d+', file_content)
            if len(geneids) == 12:
                df = pd.DataFrame(list(zip(geneids, strain_names)), columns =['GeneID', 'Strain'])
                df["SCO"] = orthogroup
                df_list.append(df)
            else:
                print("File " + file + " does not contain 12 single copy orthologs")
    
    df_concat = pd.concat(df_list, axis=0)
    df_concat.reset_index(drop=True, inplace=True)
    df_concat.to_csv('SCO_mappings.tsv', sep="\t", index=False) 


parse_SCO("./OrthoFinder/Single_Copy_Orthologue_Sequences/")
