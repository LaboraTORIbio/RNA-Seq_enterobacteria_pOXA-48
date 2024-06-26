#!/usr/bin/python3

import re
import os


### Iterating over all directories

for directory in os.listdir("."):
    if os.path.isdir(os.path.join(".", directory)):
        print(f"Reading directory: {directory}")
        
        annot_file = os.path.join(".", directory, "genomic.gff")
        prot_file = os.path.join(".", directory, "protein.faa")
        
        
        ### Parsing the annotation gff file
        
        sequence_order = {}  # Dictionary to store the sequence identifiers and their order
        chromosome_order = {}  # Dictionary to store the chromosome order
        
        if os.path.isfile(annot_file):
            with open(annot_file, "r") as f:
                for line in f:
                    if not line.startswith("#") and "CDS" in line:
                        fields = line.strip().split("\t")
                        info = fields[8].split(";")
                        seq_id = info[0].replace("ID=cds-", "")
                        seq_id = re.sub('-\d+', '', seq_id)
                        order = int(fields[3])  # the order will be based on the start position of CDSs
                        chromosome = fields[0]
                        
                        # Using a tuple (order, chromosome) as the sorting key
                        sequence_order[seq_id] = (order, chromosome)
                        chromosome_order[chromosome] = None  # store chromosome names

        
        ### Reading the protein fasta file
        
        sequences = {}  # Dictionary to store the sequences and their sequence identifiers
        
        if os.path.isfile(prot_file):
            with open(prot_file, "r") as f:
                sequence_id = ""
                sequence = ""
            
                for line in f:
                    if line.startswith(">"):
                        if sequence_id != "":
                            sequences[sequence_id] = sequence
                        sequence_id = line.strip().split(" ")
                        sequence_id = sequence_id[0].replace(">", "")
                        sequence = ""
                    else:
                        sequence += line.strip()

                # Add the last sequence to the dictionary
                sequences[sequence_id] = sequence
               
        
        ### Sorting the sequence identifiers based on genomic position and chromosome name
        sorted_ids = sorted(sequence_order, key=lambda x: sequence_order[x])
                
        ### Ordering the chromosome names alphabetically
        sorted_chromosomes = sorted(chromosome_order.keys())
        
        
        ### Writing the sorted sequences to a new file

        updated_fasta_file = directory + "/" + directory + ".faa"

        with open(updated_fasta_file, "w") as f:
            for chromosome in sorted_chromosomes:
                for seq_id in sorted_ids:
                    order, chrom = sequence_order[seq_id]
                    if chrom == chromosome:
                        seq = sequences.get(seq_id, None)
                        if seq is not None:
                            f.write(f">{seq_id}\n")
                            f.write(f"{sequences[seq_id]}\n")

        print(f"***Sequences in the protein.fasta file have been sorted. The updated sequences are stored in {updated_fasta_file}.")

