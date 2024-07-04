# Exploratory Data Analysis, Stats and Figures

## Genelarized Additive Model (GAM)

Data from Fig 6A and S13 Fig was used to build a GAM (S14 Fig) with the script **GAM.R**.


## BP and DEGs involved in iron uptake

Code to generate Figs 6CD and S15 Fig. Run from directory RNA-Seq_1st_dataset:

```sh
mkdir plot_DEGs_by_BP/

# To plot the number of DEGs by BP (Fig 6C, First RNA-Seq dataset),
# first, extract the BP terms of each DEG with the log2FC:
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv <(cat RNAseq_C325/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/C325_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv <(cat RNAseq_CF13/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f21,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/CF13_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv <(cat RNAseq_H53/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/H53_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv <(cat RNAseq_J57/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/J57_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv <(cat RNAseq_K147/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/K147_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv <(cat RNAseq_EC10/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/EC10_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv <(cat RNAseq_KPN04/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/KPN04_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv <(cat RNAseq_KPN07/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/KPN07_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $17 in a {$20=a[$17]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv <(cat RNAseq_KPN10/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f20,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/KPN10_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv <(cat RNAseq_KPN15/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f21,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/KPN15_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv <(cat RNAseq_KPN18/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f22,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/KPN18_main.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$3; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv <(cat RNAseq_MG1655/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') | cut -f21,3 | awk -F'\t' '$2 != ""' > plot_DEGs_by_BP/MG1655_main.tsv

# Problem: some genes related to iron do not have the keywords of iron in the GO annotation,
# so extract genes that have those keywords in the product name (Fig S15, First RNA-Seq dataset):
cut -f3,20,21 RNAseq_C325/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/C325_main_prod.tsv
cut -f3,19,20 RNAseq_CF13/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/CF13_main_prod.tsv
cut -f3,20,21 RNAseq_H53/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/H53_main_prod.tsv
cut -f3,20,21 RNAseq_J57/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/J57_main_prod.tsv
cut -f3,20,21 RNAseq_K147/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/K147_main_prod.tsv
cut -f3,20,21 RNAseq_EC10/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/EC10_main_prod.tsv
cut -f3,20,21 RNAseq_KPN04/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN04_main_prod.tsv
cut -f3,20,21 RNAseq_KPN07/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN07_main_prod.tsv
cut -f3,18,19 RNAseq_KPN10/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN10_main_prod.tsv
cut -f3,19,20 RNAseq_KPN15/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN15_main_prod.tsv
cut -f3,20,21 RNAseq_KPN18/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN18_main_prod.tsv
cut -f3,19,20 RNAseq_MG1655/DE_results_filtered_chromosome_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/MG1655_main_prod.tsv

# and with the Second RNA-Seq dataset:
cut -f3,20,21 ../RNA-Seq_DeltaLysR/RNAseq_CF13/DE_results_pOXA-48-vs-pOXA-48DlysR_filtered_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/CF13_pOXA-48-vs-DlysR_prod.tsv
cut -f3,20,21 ../RNA-Seq_DeltaLysR/RNAseq_KPN15/DE_results_pOXA-48-vs-pOXA-48DlysR_filtered_padj.tsv | grep -E "Fe|iron|enterobactin|ferrous|heme|siderophore|ferritin" | sed 's/"//g' > plot_DEGs_by_BP/KPN15_pOXA-48-vs-DlysR_prod.tsv

# Then run:
./plot_DEGs_by_BP/get_table_DEGs_by_BP_and_Product.py
```
The outputted tables table_count_BPs_iron.tsv, table_count_Products_iron.tsv and table_count_Products_iron_new.tsv are imported into the R script **DEGs_by_BP_plot.R** for plotting. Note that BP and genes not related to iron uptake were removed after plotting.


## Phylogeny of enterobacterial strains

**Mashtree v1.2.0** was used to generate a phylogenetic tree of mash distances between the closed genomes of the 12 pOXA-48-carrying strains, using a bootstrap of 100 (S1A Fig). Run from directory RNA-Seq_1st_dataset:

```sh
mkdir mashtree
mashtree_bootstrap.pl --reps 100 ../../Closed_sequences/C325.fasta ../../Closed_sequences/CF13.fasta ../../Closed_sequences/H53.fasta ../../Closed_sequences/J57.fasta ../../Closed_sequences/K147.fasta ../../Closed_sequences/TC_EC10.fasta ../../Closed_sequences/TC_KPN04.fasta ../../Closed_sequences/TC_KPN07.fasta ../../Closed_sequences/TC_KPN10.fasta ../../Closed_sequences/TC_KPN15.fasta ../../Closed_sequences/TC_KPN18.fasta ../../Closed_sequences/MG1655p.fasta > mashtree/mashtree.dnd
```


## Dimensionality reduction and hierarchical clustering

To further analyze common transcriptomic responses between samples, 2488 single copy orthologs (SCO; orthogroups which contain only one gene per strain) were identified with **OrthoFinder v2.5.4**.

```sh
# Run from directory RNA-Seq_1st_dataset
mkdir plot_SCO_PCA_HC/   # Copy here the .faa files of the 12 strains
~/Software/OrthoFinder/orthofinder -f plot_SCO_PCA_HC/
```

From the R script **plot_PCA_HC.R**, first run the normalization steps (VST method) of the raw counts. Then run:

```sh
cd plot_SCO_PCA_HC/
# Create a table of the single copy orthologs between the 12 strains:
./get_SCO_mappings.py
# Map those SCOs to the GeneIDs and get the normalized counts of the SCOs (pOXA-48 carriers):
./get_SCO_norm_counts_PC.py
# Map those SCOs to the GeneIDs and get the normalized counts of the SCOs (pOXA-48-free):
./get_SCO_norm_counts_PF.py
# Create a table of the log2FC and padj values of the SCOs:
./get_SCO_log2FC.py

# Percentage of SCOs that are significant (13.6%):
cat SCO_log2FC_full.tsv | cut -f5 | sort | uniq -c
# Here are the products of the SCOs:
cat SCO_log2FC_full.tsv | cut -f7 | sort | uniq -c
cd ..
```

Continue running **plot_PCA_HC.R** to perform and plot the dimensionality reduction of the RNA-Seq normalized count data (**PCA** and **t-SNE**) and to plot a heatmap with **hierarchical clustering** of the differential expression of the SCOs (S3 Fig).


## Heatmap of DEGs with BP GO terms enriched

Run the script **table4heatmaps_1.py** to generate a table with all the DEGs that have a BP GO term enriched. The table (DEGs4heatmap_full.tsv) is manually edited so all genes have a gene name for later plotting. When the gene name is missing, it is infered from the Product column or by retrieving the gene name from UniProt SwissProt. If the Product description does not contain a protein name for the gene name and the UniProt gene name is not available, the refseq IDs and the Product names of all 12 strains are compared; if the 12 strains do not share the same refseq ID for the same Product, the rows for those genes are deleted, since we cannot know if they are the same gene or not.

Then run **table4heatmaps_2.py**. The script will only retain genes that are present in all strains for comparison, and that have at least the gene significantly DE in at least one strain, to remove non-significant genes. Finally, it also adds a column with the Parental description for grouping.

```sh
# Run from directory RNA-Seq_1st_dataset
mkdir plot_heatmap/
plot_heatmap/table4heatmaps_1.py
plot_heatmap/table4heatmaps_2.py
```

Plot the heatmap with the **heatmaps.R** script (Fig S5). Notice that figure is edited on InkScape to reorder rows based on the GSEA figure (Fig 2). Note: This R script also generates the heatmap of the *pfp-ifp* operon (Fig 3A).
