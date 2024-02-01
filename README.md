# RNA-Seq_enterobacteria_pOXA-48

</br>

Code and extended Bioinformatic methods for the manuscripts:

---

Jorge Sastre-Dominguez, Javier DelaFuente, Laura Toribio-Celestino, Cristina Herencias, Pedro Herrador Gómez, Coloma Costas Romero, Rafael Cantón, Jerónimo Rodríguez-Beltrán, Alfonso Santos-Lopez, Alvaro San Millan (2024) **Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria.** *bioRxiv*

---

&

---

*In process...*

---

</br>

### Directory `Genome_assemblies`
Includes a MarkDown file describing the bioinformatics workflow for generating closed reference genomes: `closing_reference_genomes.md`, common to both manuscripts.

### Directory `RNA-Seq`
Contains bioinformatics workflows for:
1. RNA-Seq read quality control, described in the MarkDown file `read_quality_control.md`, with associated files in subdirectory `reads_RNAseq`, common to both manuscripts.
2. Main RNA-Seq data analysis (*for manuscript in process...*)

### Directory `RNA-Seq_experimental_evolution`
Contains bioinformatics workflows specific to the manuscript *Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria*:
1. The MarkDown file `RNA-Seq_experimental_evolution.md` describes the general RNA-Seq analysis workflow
2. Each `RNAseq_<strain>` directory includes:
   * The `quant_diffexpr.Rmd` script to perform read quantification and differential expression analysis
   * The rendered HTML file `quant_diffexpr.html` to visualize the outputs of the Rmd script
   * Other relevant files for the analysis (see the general workflow)
