#!/usr/bin/python3

import os
import re


def read_all_systems(macsyfinder_path):
    """Reads the all_systems.txt file and retrieves the ID before 'pirin_PF02678'"""
    with open(macsyfinder_path, 'r') as file:
        for line in file:
            if line.startswith('clusters = '):
                line = line.strip()
                line = re.sub('clusters = ', '', line)
                line = re.sub('\[|\]|\(|\)|\'', '', line)
                lsplit = line.split(", ")
                ind = lsplit.index('pirin_PF02678')
                id_pirin = lsplit[ind-1]
                return id_pirin
    return None


def extract_context_from_gff(gff_file, id_pirin, lines_before=5, lines_after=5):
    """Searches the GFF file for the target ID, extracts surrounding lines, and checks if LysR is adjacent."""
    with open(gff_file, 'r') as file:
        lines = file.readlines()

    for i, line in enumerate(lines):
        if id_pirin in line and 'CDS' in line:
            # Extract surrounding lines
            start = max(0, i - lines_before)
            end = min(len(lines), i + lines_after + 1)
            surrounding_lines = lines[start+1:end:2]
            
            # Find the line with LysR in surrounding lines
            lysr_line = next((l for l in surrounding_lines if 'LysR' in l), None)
            
            if lysr_line:
                # Check if LysR is adjacent to the id_pirin (either before or after)
                lysr_index = lines.index(lysr_line)
                lysr_adjacent_pfp = abs(i - lysr_index) == 2

                # Check if "isochorismatase" or "hydrolase" is in the lines next to the LysR line
                start_lysr = max(0, lysr_index - 2)
                end_lysr = min(len(lines), lysr_index + 3)
                adjacent_lines = lines[start_lysr:end_lysr]
                lysr_adjacent_ifp = any('isochorismatase' in l or 'hydrolase' in l for l in adjacent_lines)

                return surrounding_lines, lysr_line, lysr_adjacent_pfp, lysr_adjacent_ifp

            else:
                return surrounding_lines, None, False, False  # LysR not found

    return None, None, False, False  # Target ID not found



# Change Proteobacteria classes here!!!
strains = 'operon_wo_lysR/gamma.txt'  # 
macsyfinder_dir = 'macsyfinder_operon_v2/gamma'
gff_dir = 'ncbi_dataset_gamma/'

macsyfinder_file = 'all_systems.txt'
gff_file = 'genomic.gff'


num_strains = 0
num_lysr = 0
num_frameshift = 0
num_lysR_adjacent_pfp = 0
num_lysR_adjacent_ifp = 0

with open(strains, 'r') as file_strains:
    for strain in file_strains:
        num_strains +=1
        strain = strain.strip()
        print(f"\n\n############# {strain} #############\n")
        macsy_strain_dir = os.path.join(macsyfinder_dir, strain, macsyfinder_file)
        gff_strain_dir = os.path.join(gff_dir, strain, gff_file)
        
        if os.path.isfile(macsy_strain_dir):
            id_pirin = read_all_systems(macsy_strain_dir)
            if id_pirin:
                surrounding_lines, lysr_line, lysr_adjacent_pfp, lysr_adjacent_ifp = extract_context_from_gff(gff_strain_dir, id_pirin)
                if surrounding_lines:
                    print("".join(surrounding_lines))
                    if lysr_line:
                        print(f"LysR found in the following line:\n{lysr_line}")
                        num_lysr +=1
                        words = ["frameshift", "incomplete", "internal stop", "pseudo", "pseudogene", "disrupted",
                                 "partial", "truncated", "stop codon", "nonsense", "split gene", "degenerat"]
                        if any(x in lysr_line for x in words):
                            num_frameshift +=1
                        if lysr_adjacent_pfp:
                            num_lysR_adjacent_pfp +=1
                        if lysr_adjacent_ifp:
                            num_lysR_adjacent_ifp +=1
                    else:
                        print("LysR not found.")
                else:
                    print(f"{id_pirin} not found in {gff_strain_dir}.")
            else:
                print(f"{id_pirin} not found in {macsy_strain_dir}.")
        else:
            print(f"File {macsy_strain_dir} does not exist.")
            

print("\n\n############################################")
print(f"Total number of strains = {num_strains}")
print(f"{num_lysr} strains encode a LysR within +-2 genes distance from pfp ({num_lysR_adjacent_pfp} adjacent to pfp, {num_lysR_adjacent_ifp} adjacent to ifp)")
print(f"Of these LysRs, {num_frameshift} are knocked-out or frameshifted.")


