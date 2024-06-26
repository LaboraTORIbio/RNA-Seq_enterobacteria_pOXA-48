#!/usr/bin/python3

import re
import os


# Change here the path of MacSyFinder results and NCBI datasets (gamma, alpha or beta):
classdir = "macsyfinder_operon/gamma/"
dbdir = "ncbi_dataset_gamma/"



### Iterating over each MacSyFinder directories

for strain in os.listdir(classdir):
    if os.path.isdir(classdir):
        
        # Saving paths
        resfile = os.path.join(classdir, strain, "all_systems.txt")
        annotfile = os.path.join(dbdir, strain, "genomic.gff")
        
        # Here we will store the gene ID of the first and last genes of the system
        firstgene = ""
        lastgene = ""
        
        
        ### Opening the MacSyFinder result file
        
        if os.path.isfile(resfile):
            #print(f"Reading MacSyFinder file of: {strain}")
            with open(resfile, "r") as f:
                for line in f:
                    
                    # Selecting the line that contains the order of genes of the system
                    if line.startswith("clusters = "):
                        line = line.strip()
                        line = re.sub('clusters = ', '', line)
                        line = re.sub('\[|\]|\(|\)|\'', '', line)
                        lsplit = line.split(", ")
                        regex = re.compile(r'^\d+')
                        genes = [i for i in lsplit if not regex.match(i)]
                        firstgene = genes[0]
                        lastgene = genes[-2]
                        

        ### Opening the GFF annotation file
        
        if os.path.isfile(annotfile):
            #print(f"\nParsing annotation file for: {strain}")
            with open(annotfile, "r") as f:

                # Flag to track if we are inside the system region
                inside_system = False
                strand = []
                
                for line in f:
                    # Set the flag to True when the first gene is found
                    if firstgene in line:
                        inside_system = True
                    
                    # Printing the lines of the system
                    if inside_system:
                        if 'CDS' in line:
                            line = line.strip()
                            lsplit = line.split("\t")
                            info = lsplit[8].split(";")
                            strand.append(lsplit[6])
                            if "product=" in lsplit[8]:
                                proditem = list(filter(lambda x: "product=" in x, info))
                                proditem = proditem[0]
                                prodsplit = proditem.split("=")
                                product = prodsplit[-1]
                                if "Name=" in lsplit[8]:
                                    nameitem = list(filter(lambda x: "Name=" in x, info))
                                    nameitem = nameitem[0]
                                    namesplit = nameitem.split("=")
                                    name = namesplit[-1]
                                
                    # Set the flag to False when the last gene is found
                    if lastgene in line:
                        inside_system = False
                        break # Exit the loop since we found the last gene
                
                
                ### Check strains that parsed incorrectly or has duplicated genes
                if len(strand) > 3:
                    print(strain + "\t" + str(len(strand)))

