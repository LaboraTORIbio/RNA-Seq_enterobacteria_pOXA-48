#!/usr/bin/python3

import argparse

# Redirect output to a file

parser = argparse.ArgumentParser()
parser.add_argument('filename')
args = parser.parse_args()

with open(args.filename, 'r') as inputfile:
    header = "GeneID\tChr\tStart\tEnd\tStrand\tType\tRefSeq\tGene\tProduct"
    print(header)
    for line in inputfile:
        if line[0].isdigit():
            line = line.rstrip()
            columns = line.split("\t")
            # Saving the gene IDs from the features' info column
            featuresall = columns[8]
            featuresplit = featuresall.split(";")
            genefield1 = featuresplit[0]
            genefield2 = genefield1.split("=")
            geneid = genefield2[1]
            items = {'tRNA', 'ncRNA', 'tmRNA', 'antisense_RNA', 'RNase_P_RNA', 'SRP_RNA'} # Do not include rRNA to remove them from DE analysis
            
            # Reading only lines that contain CDSs
            if columns[2] == "CDS":
                # Reading lines that contain RefSeq IDs
                if "RefSeq" in featuresall:
                    refseqitem = list(filter(lambda x: "RefSeq" in x, featuresplit))
                    refseqitem = refseqitem[0]
                    refseqsplit = refseqitem.split(":")
                    refseq = refseqsplit[-1]
                    # Reading lines that contain gene products
                    if "product=" in featuresall:
                        proditem = list(filter(lambda x: "product=" in x, featuresplit))
                        proditem = proditem[0]
                        prodsplit = proditem.split("=")
                        product = prodsplit[-1]
                        # Reading lines that contain gene names
                        if "gene=" in featuresall:
                            geneitem = list(filter(lambda x: "gene=" in x, featuresplit))
                            geneitem = geneitem[0]
                            genesplit = geneitem.split("=")
                            gene = genesplit[-1]
                            print(geneid + "\t" + columns[0] + "\t" + columns[3] + "\t" + columns[4] + "\t" + columns[6] + "\t" + "CDS" + "\t" + refseq + "\t" + gene + "\t" + product)
                        # Printing CDS that don't contain a gene name
                        else:
                            print(geneid + "\t" + columns[0] + "\t" + columns[3] + "\t" + columns[4] + "\t" + columns[6] + "\t" + "CDS" + "\t" + refseq + "\t" + "-" + "\t" + product)
                    # Printing CDS that don't contain gene product
                    else:
                        print(geneid + "\t" + columns[0] + "\t" + columns[3] + "\t" + columns[4] + "\t" + columns[6] + "\t" + "CDS" + "\t" + refseq + "\t" + "-" + "\t" + "-")
                # Printing CDS that don't contain a RefSeq ID
                else:
                    # Reading lines that contain gene products
                    if "product=" in featuresall:
                        proditem = list(filter(lambda x: "product=" in x, featuresplit))
                        proditem = proditem[0]
                        prodsplit = proditem.split("=")
                        product = prodsplit[-1]
                        print(geneid + "\t" + columns[0] + "\t" + columns[3] + "\t" + columns[4] + "\t" + columns[6] + "\t" + "CDS" + "\t" + "-" + "\t" + "-" + "\t" + product)	
               
            # Lines that contain other types of features (e.g. tRNA)
            elif columns[2] in items:
                print(geneid + "\t" + columns[0] + "\t" + columns[3] + "\t" + columns[4] + "\t" + columns[6] + "\t" + columns[2] + "\t" + "-" + "\t" + "-" + "\t" + "-")


