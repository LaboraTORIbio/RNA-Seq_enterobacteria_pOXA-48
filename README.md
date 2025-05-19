# RNA-Seq_enterobacteria_pOXA-48

</br>

**Code and extended Bioinformatic methods for the manuscripts:**

</br>


---

Laura Toribio-Celestino, Alicia Calvo-Villamañán, Cristina Herencias, Aida Alonso-del Valle, Jorge Sastre-Dominguez, Susana Quesada, Didier Mazel, Eduardo PC Rocha, Ariadna Fernández-Calvet, Alvaro San Millan (2024) **A plasmid-chromosome crosstalk in multidrug resistant clinical enterobacteria.** *Nat Commun* 15, 10859 doi: [https://doi.org/10.1101/2024.08.08.607126](https://doi.org/10.1038/s41467-024-55169-y).

---


### 📂 Directory [`Genome_assemblies`](./Genome_assemblies/)
Includes a MarkDown file describing the bioinformatics workflow for generating closed reference genomes: [`closing_reference_genomes.md`](./Genome_assemblies/closing_reference_genomes.md).


### 📂 Directory [`RNA-Seq_1st_dataset`](./RNA-Seq_1st_dataset/)
Contains bioinformatics workflows and scripts for:

**RNA-Seq read quality control.** Described in the MarkDown file [`read_quality_control.md`](./RNA-Seq_1st_dataset/read_quality_control.md), with associated files in subdirectory `reads_RNAseq`.

**RNA-Seq data analysis of the 1st dataset.** The general RNA-Seq analysis workflow is described in [`RNA-Seq_1st_dataset.md`](./RNA-Seq_1st_dataset/RNA-Seq_1st_dataset.md). Each `RNAseq_<strain>` subdirectory includes:
* The `quant_diffexpr.Rmd` script to perform read quantification and differential expression analysis
* The HTML file `quant_diffexpr.html` to visualize the outputs of the Rmd script (see the general workflow for links to the rendered HTMLs)
* The `TPM_calculation.Rmd` script to calcultate the TPM values for each strain

**Gene set enrichment analysis (GSEA)**. Workflow is also described in [`RNA-Seq_1st_dataset.md`](./RNA-Seq_1st_dataset/RNA-Seq_1st_dataset.md). Analysis rendered [here](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/GSEA/GSEA.html).

**Data visualization.** TPM plots, heatmaps of gene expression, plots of DEGs by BP, dimensionality reduction and hierarchical clustering.


### 📂 Directory [`RNA-Seq_DeltaLysR`](./RNA-Seq_DeltaLysR/)
Contains bioinformatics workflows and scripts for:

**RNA-Seq read quality control and data analysis.** The general RNA-Seq analysis workflow is described in [`RNA-Seq_2nd_dataset.md`](./RNA-Seq_DeltaLysR/RNA-Seq_2nd_dataset.md). Each `RNAseq_<strain>` subdirectory includes:
* The `quant_diffexpr.Rmd` script to perform read quantification and differential expression analysis
* The HTML file `quant_diffexpr.html` to visualize the outputs of the Rmd script (see the general workflow for links to the rendered HTMLs)
* Other relevant files for the analysis (see the general workflow)

**Data visualization.** Plots of the differential expression of genes regulated by LysRs.


### 📂 Directory [`LysR_and_the_mystery_operon`](./LysR_and_the_mystery_operon/)
Contains bioinformatics workflows and scripts for the study of the *lysR-pfp-ifp* cluster and the *pfp-ifp* operon. The general workflow is described in [`LysR_and_the_mystery_operon.md`](./LysR_and_the_mystery_operon/LysR_and_the_mystery_operon.md).
* Construction of MysteryOperonFinder and MysteryOperonFinder_v2 - MacSyFinder models
* Identification of the *lysR-pfp-ifp* cluster in a large database of Proteobacteria genomes
* Identification of *pfp-ifp* clusters without the associated *lysR*
* Association analysis of *lysR-pfp-ifp* clusters and pOXA-48 carriage
* Construction of phylogenetic trees of the *lysR-pfp-ifp* cluster and Proteobacteria
* Ancestral character reconstruction of the *lysR-pfp-ifp* cluster
* Genomic neighborhood of the *lysR-pfp-ifp* cluster
* GC content and CAI calculation
* Analysis of the conservation of LysR<sub>pOXA-48</sub>


### 📂 Directory [`EDA_Stats_Figs`](./EDA_Stats_Figs/)
Includes workflows for exploratory data analyses, generating plots and performing statistical tests, including a generalized additive model (GAM).


</br>

---

Jorge Sastre-Dominguez, Javier DelaFuente, Laura Toribio-Celestino, Cristina Herencias, Pedro Herrador Gómez, Coloma Costas Romero, Rafael Cantón, Jerónimo Rodríguez-Beltrán, Alfonso Santos-Lopez, Alvaro San Millan (2024) **Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria.** *Nat Ecol Evol* 8, 2097–2112 doi: [https://doi.org/10.1101/2024.03.01.582297](https://doi.org/10.1038/s41559-024-02523-4).

---

**Generation of closed reference genomes** described in [`closing_reference_genomes.md`](./Genome_assemblies/closing_reference_genomes.md).

**RNA-Seq read quality control** described in [`read_quality_control.md`](./RNA-Seq_1st_dataset/read_quality_control.md).

### 📂 Directory [`RNA-Seq_experimental_evolution`](./RNA-Seq_experimental_evolution/)
Contains bioinformatics workflows and scripts for:

**RNA-Seq read quality control and data analysis.** The general RNA-Seq analysis workflow is described in [`RNA-Seq_experimental_evolution.md`](./RNA-Seq_experimental_evolution/RNA-Seq_experimental_evolution.md). Each `RNAseq_<strain>` subdirectory includes:
* The `quant_diffexpr.Rmd` script to perform read quantification and differential expression analysis (also available in R format)
* The HTML file `quant_diffexpr.html` to visualize the outputs of the Rmd script (see the general workflow for links to the rendered HTMLs)
* Other relevant files for the analysis (see the general workflow)
