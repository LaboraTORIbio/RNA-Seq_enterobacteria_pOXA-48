#!/usr/bin/python3

import re
import os


# Change here the path of MacSyFinder results (gamma, alpha or beta):
classdir = "macsyfinder_operon/gamma/"


# Iterating over each MacSyFinder directories
for strain in os.listdir(classdir):
    if os.path.isdir(classdir):
        #print(f"Reading directory: {strain}")
        
        resfile = os.path.join(classdir, strain, "all_systems.txt")
        
        # Opening result file
        if os.path.isfile(resfile):
            with open(resfile, "r") as f:
                for line in f:
                    
                    # Selecting line that contains the order of genes
                    if line.startswith("clusters = "):
                        line = line.strip()
                        line = re.sub('clusters = ', '', line)
                        line = re.sub('\[|\]|\(|\)|\'', '', line)
                        lsplit = line.split(", ")
                        genes = [i for i in lsplit if ('PF' in i)]
                        genes_str = "_".join(genes)
                        genes_str = re.sub('_PF\d{5}', '', genes_str)
                        
                        # Printing strains with genes in abnormal order for further manual inspection
                        if not genes_str.startswith("lysR") and not genes_str.endswith("lysR") or len(genes_str) > 26:
                            print(strain + "\t" + genes_str)
