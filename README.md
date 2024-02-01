# RNA-Seq_enterobacteria_pOXA-48

Code and extended Bioinformatic methods for the manuscripts:

---

Jorge Sastre-Dominguez, Javier DelaFuente, Laura Toribio-Celestino, Cristina Herencias, Pedro Herrador Gómez, Coloma Costas Romero, Rafael Cantón, Jerónimo Rodríguez-Beltrán, Alfonso Santos-Lopez, Alvaro San Millan (2024) **Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria.** *bioRxiv*

---

&

---

*In process...*

---

1. Directory `Genome_assemblies` includes a MarkDown file describing the bioinformatics workflow for generating closed reference genomes: `closing_reference_genomes.md`
2. Directory `RNA-Seq` contains bioinformatics workflows for:
   1. RNA-Seq read quality control, described in the MarkDown file `read_quality_control.md`, with associated files in subdirectory `reads_RNAseq`
   2. Main RNA-Seq data analysis (*in process...*)
3. Directory `RNA-Seq_experimental_evolution` contains bioinformatics workflows specific to the manuscript *Plasmid-encoded insertion sequences promote rapid adaptation in clinical enterobacteria*:
   1. The MarkDown file `RNA-Seq_experimental_evolution.md` describes the general RNA-Seq analysis workflow
   2. Each `RNAseq_<strain>` directory includes:
      * The `quant_diffexpr.Rmd` script to perform read quantification and differential expression analysis in R
      * The rendered HTML file `quant_diffexpr.html` to visualize the outputs of the Rmd script
      * Other relevant files for the analysis (see the general workflow)
