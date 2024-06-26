#!/usr/bin/python3

import re
import os
from Bio import SeqIO


### Change here the path of MacSyFinder results NCBI datasets and output directory (gamma, alpha or beta):

classdir = "macsyfinder_operon/gamma/"
dbdir = "ncbi_dataset_gamma/"

output_dir = "operon_seqs/gamma/"
os.makedirs(output_dir, exist_ok=True)



### Iterating over each MacSyFinder directories

for strain in os.listdir(classdir):
    if os.path.isdir(classdir):
        
        resfile = os.path.join(classdir, strain, "all_systems.txt")
        faafile = dbdir + strain + "/" + strain + ".faa"
        
        dict_id_gene = {}
        
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
                        dict_id_gene = {lsplit[0] : re.sub('_PF\d{5}', '', lsplit[1]), lsplit[3] : re.sub('_PF\d{5}', '', lsplit[4]), lsplit[6] : re.sub('_PF\d{5}', '', lsplit[7])}
        
        
        #### Opening and reading the protein fasta file
        
        concat_seqs = {}
        
        if os.path.isfile(faafile):
            for record in SeqIO.parse(faafile, "fasta"):
                sequence_id = record.id
                if sequence_id in dict_id_gene:  # If the sequence_id is in the dictionary, store sequence
                    gene_name = dict_id_gene[sequence_id]
                    sequence = str(record.seq)
                    concat_seqs[gene_name] = sequence
                    
                    # Create an individual FASTA file for each sequence
                    output_file1 = os.path.join(output_dir, f"{strain}_{sequence_id}_{gene_name}.fasta")
                    with open(output_file1, "w") as output_fasta1:
                        
                        # Change here also the Proteobacteria class!!
                        output_fasta1.write(">gamma_" + strain + "\n" + sequence)
        
        # Concatenate seqs of lysR-pirin-isochorismatase
        output_file2 = os.path.join(output_dir, f"{strain}_concat.fasta")
        if concat_seqs:
            with open(output_file2, "w") as output_fasta2:
                output_fasta2.write(">" + strain + "\n" + concat_seqs["lysR"] + concat_seqs["pirin"] + concat_seqs["isochorismatase"])
        
        
        print(f"Stored sequences of {strain} in {output_dir}")

