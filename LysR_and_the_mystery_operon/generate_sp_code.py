#!/usr/bin/python3

import pandas as pd

# Read the table
df = pd.read_csv("data_summary.tsv", sep='\t')

# Define a function to generate a unique number for each combination of genus and species
def generate_unique_number(group):
    group['unique_number'] = group.groupby(['genus', 'species']).cumcount() + 1
    return group

# Extract the genus and species names
df['genus'] = df['Organism Scientific Name'].str.split().str[0].str[:3].str.upper()
df['species'] = df['Organism Scientific Name'].str.split().str[1].str[:3].str.upper()

# Apply the generate_unique_number function to create the unique number column
df = df.groupby(['genus', 'species']).apply(generate_unique_number)

# Create the new column with the desired format
df['Species Code'] = df['genus'] + df['species'] + '.0723.' + df['unique_number'].astype(str).str.zfill(5)

# Edit the final table
df.drop(columns=['genus', 'species', 'unique_number'], inplace=True)
cols_order = ['Species Code', 'Organism Scientific Name', 'Organism Common Name', 'Organism Qualifier', 'Taxonomy id', 'Assembly Name', 'Assembly Accession', 'Source', 'Annotation', 'Level', 'Contig N50', 'Size', 'Submission Date', 'Gene Count', 'BioProject', 'BioSample']
df = df[cols_order]

# Save the updated table to a new CSV file named "output.csv"
output_file = "data_information.tsv"
df.to_csv(output_file, index=False, sep='\t')

