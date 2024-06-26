#! /usr/bin/python3

import pandas as pd
import glob


files_main = glob.glob('plot_DEGs_by_BP/*main.tsv')
files_main_prod = glob.glob('plot_DEGs_by_BP/*main_prod.tsv')
files_new_prod = glob.glob('plot_DEGs_by_BP/*_pOXA-48-vs-DlysR_prod.tsv')

main_df = []
main_df_prod = []
new_df_prod = []

for file in files_main:
    strain_name = file.replace('plot_DEGs_by_BP/', '').replace('_main.tsv', '')
    df = pd.read_csv(file, sep='\t', names=['log2FC', 'GO_BP_list'])
    df['Strain'] = strain_name
    main_df.append(df)

for file in files_main_prod:
    strain_name = file.replace('plot_DEGs_by_BP/', '').replace('_main_prod.tsv', '')
    df = pd.read_csv(file, sep='\t', names=['log2FC', 'Gene', 'Product'])
    df['Strain'] = strain_name
    main_df_prod.append(df)

for file in files_new_prod:
    strain_name = file.replace('plot_DEGs_by_BP/', '').replace('_prod.tsv', '')
    df = pd.read_csv(file, sep='\t', names=['log2FC', 'Gene', 'Product'])
    df['Strain'] = strain_name
    new_df_prod.append(df)


main_df = pd.concat(main_df, axis=0)
main_df.reset_index(drop=True, inplace=True)

main_df_prod = pd.concat(main_df_prod, axis=0)
main_df_prod.reset_index(drop=True, inplace=True)

new_df_prod = pd.concat(new_df_prod, axis=0)
new_df_prod.reset_index(drop=True, inplace=True)


# Split GO terms by ; and make them new rows

split_df = main_df.iloc[:,1].str.split(';', expand=True).stack().reset_index(level=1, drop=True).rename('GO_BP')
table = pd.merge(split_df, main_df, left_index=True, right_index=True).reset_index(drop=True)


# Reformat tables

table = table.drop('GO_BP_list', axis=1)
table.iloc[:,0] = table.iloc[:,0].str.replace(r'\[GO:[0-9]{2,10}\]', '', regex=True)
table.iloc[:,0] = table.iloc[:,0].str.replace(r'^ ', '', regex=True)
table['Direction'] = table.iloc[:,1] > 0
table['Direction'] = table['Direction'].replace({True: 'upregulated', False: 'downregulated'})

mask_main_df_prod = main_df_prod['Gene'] != '-'
main_df_prod.loc[mask_main_df_prod, 'Product'] += ' (' + main_df_prod.loc[mask_main_df_prod, 'Gene'] + ')'
main_df_prod = main_df_prod.drop('Gene', axis=1)
main_df_prod['Direction'] = main_df_prod.iloc[:,0] > 0
main_df_prod['Direction'] = main_df_prod['Direction'].replace({True: 'upregulated', False: 'downregulated'})

mask_new_df_prod = new_df_prod['Gene'] != '-'
new_df_prod.loc[mask_new_df_prod, 'Product'] += ' (' + new_df_prod.loc[mask_new_df_prod, 'Gene'] + ')'
new_df_prod = new_df_prod.drop('Gene', axis=1)
new_df_prod['Direction'] = new_df_prod.iloc[:,0] > 0
new_df_prod['Direction'] = new_df_prod['Direction'].replace({True: 'upregulated', False: 'downregulated'})


# Get the count of each BP
table_full = table.groupby(['GO_BP', 'Direction']).size().reset_index(name='Count')

# Get the count of each BP by strain
table_full_by_strain = table.groupby(['GO_BP', 'Strain', 'Direction']).size().reset_index(name='Count')


# Get BPs related to iron
filter_by = 'Fe|iron|enterobactin|ferrous|heme|siderophore|siderophore-iron'

table_iron = table_full_by_strain[table_full_by_strain['GO_BP'].str.contains(filter_by, case=True)]
table_iron_copy = table_iron.copy()
table_iron_copy.loc[table_iron_copy['Direction'] == 'downregulated', 'Count'] *= -1


# Get the count of each gene product with functions related to iron

table_main_prod = main_df_prod.groupby(['Product', 'Strain', 'Direction']).size().reset_index(name='Count')
table_main_prod_copy = table_main_prod.copy()
table_main_prod_copy.loc[table_main_prod_copy['Direction'] == 'downregulated', 'Count'] *= -1

table_new_prod = new_df_prod.groupby(['Product', 'Strain', 'Direction']).size().reset_index(name='Count')
table_new_prod_copy = table_new_prod.copy()
table_new_prod_copy.loc[table_new_prod_copy['Direction'] == 'downregulated', 'Count'] *= -1


# Saving tables

table_iron_copy.to_csv('plot_DEGs_by_BP/table_count_BPs_iron.tsv', sep='\t', index=False)
table_main_prod_copy.to_csv('plot_DEGs_by_BP/table_count_Products_iron.tsv', sep='\t', index=False)
table_new_prod_copy.to_csv('plot_DEGs_by_BP/table_count_Products_iron_new.tsv', sep='\t', index=False)

writer = pd.ExcelWriter('plot_DEGs_by_BP/table_DEGs_by_BP.xlsx', engine='xlsxwriter')
table_full.to_excel(writer, index=False, sheet_name='Full')
table_full_by_strain.to_excel(writer, index=False, sheet_name='Full_by_strain')
writer.save()

