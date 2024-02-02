# RNA-Seq Analyses for the Experimental Evolution Project

## Project setup

Before starting, check how reference genomes were assembled and closed ([closing_reference_genomes.md](Genome_assemblies/closing_reference_genomes.md)) and how RNA-Seq reads were processed ([read_quality_control.md](RNA-Seq/read_quality_control.md)).

In this project, we will not use the RNA-Seq data of the *E. coli* strain MG1655 since we are interested in clinical strains. Also, we will use the alternative strain nomenclature for the pOXA-48-naive strains, to match that of the main manuscript. Therefore: EC10 = C063, KPN04 = K091, KPN07 = K141, KPN10 = K209, KPN15 = K249 and KPN18 = K275. Here, the suffix *p* denotes the pOXA-48-carrying strain.

Create the project directory from which run all commands:

```sh
mkdir -p RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution
cd RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution
```

Here, each strain will have its own directory, denoted with the prefix *RNAseq* (e.g. RNAseq_C325), that will contain BAM mapping files and the Rproject and necessary files for conducting differential expression (DE) analysis in R.


## Identifying insertion sequences (ISs)

We analyzed our RNA-Seq data to study the transcriptional regulation of the IS1. ISs typically bear multiple copies with  high identity nucleotide sequences between them, and thus, short Illumina reads can map ambiguously to different IS copies in different locations of the genome, reducing the statistical power to detect differential expression of ISs between conditions (carrying/not-carrying pOXA-48). To avoid this, we masked IS sequences from the reference genomes, leaving only one IS copy per IS family unmasked. This way, reads belonging to an IS family (like IS1) can map with low ambiguity to one region of the genome.

We **inspected the annotations of the reference genomes to find all IS copies of each IS family**, as annotated by PGAP. Since strain J57 does not encode IS1 copies in the chromosome or in other plasmids, only in pOXA-48, the differential expression analysis could show an acute over-expression of IS1 in the pOXA-48-carrying strains, possibly biasing the correlation (see main manuscript). Therefore, this strain was excluded from further analyses.

## Masking IS regions

The files **masked_regions.tsv** were manually generated and placed at each strain's *RNAseq* directory. These are BED format tables in which the columns are the chromosome, start and end positions of the ISs that will be masked (adding Ns). Note that one copy per IS family is not included in the table and thus will not be masked. These unmasked copies per IS family are preferentially located at the chromosome. Masking is performed with the **bedtools maskfasta v2.27.1** tool:

```sh
bedtools maskfasta -fi ../../Closed_sequences/C325.fasta -fo RNAseq_C325/C325_maskedIS.fasta -bed RNAseq_C325/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/CF13.fasta -fo RNAseq_CF13/CF13_maskedIS.fasta -bed RNAseq_CF13/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/H53.fasta -fo RNAseq_H53/H53_maskedIS.fasta -bed RNAseq_H53/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/K147.fasta -fo RNAseq_K147/K147_maskedIS.fasta -bed RNAseq_K147/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_EC10.fasta -fo RNAseq_C063/C063p_maskedIS.fasta -bed RNAseq_C063/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_KPN04.fasta -fo RNAseq_K091/K091p_maskedIS.fasta -bed RNAseq_K091/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_KPN07.fasta -fo RNAseq_K141/K141p_maskedIS.fasta -bed RNAseq_K141/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_KPN10.fasta -fo RNAseq_K209/K209p_maskedIS.fasta -bed RNAseq_K209/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_KPN15.fasta -fo RNAseq_K249/K249p_maskedIS.fasta -bed RNAseq_K249/masked_regions.tsv
bedtools maskfasta -fi ../../Closed_sequences/TC_KPN18.fasta -fo RNAseq_K275/K275p_maskedIS.fasta -bed RNAseq_K275/masked_regions.tsv
```

Then, the masked FASTAs are reannotated with **PGAP v2021-07-01.build5508** and the GFF3 annotation file is placed at each strain's *RNAseq* directory.


## Mapping RNA-Seq reads

Since splicing is rare in bacteria, we can use a regular aligner to map the RNA-Seq reads to their respective masked reference genomes. Here, we used **BWA-MEM v0.7.17**:

```sh
bwa index RNAseq_C325/C325_maskedIS.fasta
bwa index RNAseq_CF13/CF13_maskedIS.fasta
bwa index RNAseq_H53/H53_maskedIS.fasta
bwa index RNAseq_K147/K147_maskedIS.fasta
bwa index RNAseq_C063/C063p_maskedIS.fasta
bwa index RNAseq_K091/K091p_maskedIS.fasta
bwa index RNAseq_K141/K141p_maskedIS.fasta
bwa index RNAseq_K209/K209p_maskedIS.fasta
bwa index RNAseq_K249/K249p_maskedIS.fasta
bwa index RNAseq_K275/K275p_maskedIS.fasta

# C325
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325.1_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325.1.bam
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325.2_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325.2.bam
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325.3_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325.3.bam
# C325c1
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.1_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325c1.1.bam
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.2_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325c1.2.bam
bwa mem RNAseq_C325/C325_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C325c1.3_val_2.fq.gz | samtools sort -o ./RNAseq_C325/C325c1.3.bam

# CF13
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13.1_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13.1.bam
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13.2_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13.2.bam
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13.3_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13.3.bam
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13.4_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13.4_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13.4.bam
# CF13c1
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13c1.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13c1.1_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13c1.1.bam
bwa mem RNAseq_CF13/CF13_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_CF13c1.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_CF13c1.2_val_2.fq.gz | samtools sort -o ./RNAseq_CF13/CF13c1.2.bam

# H53
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53.1_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53.1.bam
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53.2_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53.2.bam
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53.3_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53.3.bam
# H53c1
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.1_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53c1.1.bam
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.2_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53c1.2.bam
bwa mem RNAseq_H53/H53_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_H53c1.3_val_2.fq.gz | samtools sort -o ./RNAseq_H53/H53c1.3.bam

# K147
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147.1_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147.1.bam
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147.2_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147.2.bam
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147.3_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147.3.bam
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147.4_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147.4_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147.4.bam
# K147c1
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147c1.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147c1.1_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147c1.1.bam
bwa mem RNAseq_K147/K147_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K147c1.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K147c1.2_val_2.fq.gz | samtools sort -o ./RNAseq_K147/K147c1.2.bam

# C063
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063.1_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063.1.bam
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063.2_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063.2.bam
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063.3_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063.3.bam
# C063p
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063p.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063p.1_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063p.1.bam
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063p.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063p.2_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063p.2.bam
bwa mem RNAseq_C063/C063p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_C063p.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_C063p.3_val_2.fq.gz | samtools sort -o ./RNAseq_C063/C063p.3.bam

# K091
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091.1_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091.1.bam
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091.2_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091.2.bam
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091.3_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091.3.bam
# K091p
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091p.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091p.1_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091p.1.bam
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091p.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091p.2_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091p.2.bam
bwa mem RNAseq_K091/K091p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_K091p.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_K091p.3_val_2.fq.gz | samtools sort -o ./RNAseq_K091/K091p.3.bam

# K141
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.1_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141.1.bam
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.2_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141.2.bam
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN07.3_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141.3.bam
# K141p
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.1_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141p.1.bam
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.2_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141p.2.bam
bwa mem RNAseq_K141/K141p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN07.3_val_2.fq.gz | samtools sort -o ./RNAseq_K141/K141p.3.bam

# K209
bwa mem RNAseq_K209/K209p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN10.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN10.1_val_2.fq.gz | samtools sort -o ./RNAseq_K209/K209.1.bam
bwa mem RNAseq_K209/K209p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN10.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN10.2_val_2.fq.gz | samtools sort -o ./RNAseq_K209/K209.2.bam
# K209p
bwa mem RNAseq_K209/K209p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN10.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN10.1_val_2.fq.gz | samtools sort -o ./RNAseq_K209/K209p.1.bam
bwa mem RNAseq_K209/K209p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN10.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN10.2_val_2.fq.gz | samtools sort -o ./RNAseq_K209/K209p.2.bam

# K249
bwa mem RNAseq_K249/K249p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN15.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN15.1_val_2.fq.gz | samtools sort -o ./RNAseq_K249/K249.1.bam
bwa mem RNAseq_K249/K249p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN15.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN15.2_val_2.fq.gz | samtools sort -o ./RNAseq_K249/K249.2.bam
# K249p
bwa mem RNAseq_K249/K249p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.1_val_2.fq.gz | samtools sort -o ./RNAseq_K249/K249p.1.bam
bwa mem RNAseq_K249/K249p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.2_val_2.fq.gz | samtools sort -o ./RNAseq_K249/K249p.2.bam
bwa mem RNAseq_K249/K249p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN15.3_val_2.fq.gz | samtools sort -o ./RNAseq_K249/K249p.3.bam

# K275
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.1_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275.1.bam
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.2_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275.2.bam
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_PF_KPN18.3_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275.3.bam
# K275p
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.1_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.1_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275p.1.bam
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.2_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.2_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275p.2.bam
bwa mem RNAseq_K275/K275p_maskedIS.fasta ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.3_val_1.fq.gz ../RNA-Seq/reads_RNAseq/WTCHG_TC_KPN18.3_val_2.fq.gz | samtools sort -o ./RNAseq_K275/K275p.3.bam
```

The percentage of mapped reads was >99% for all samples (**samtools v1.12**):

```sh
for bam in ./RNAseq_*/*bam
do
	echo "####################################################"
	echo "MAPPING STATS OF $bam"
	samtools flagstats $bam
done
```


## Differential expression analysis

The first step in DE analysis is counting the number of reads mapped to each genomic feature, which is performed here with featureCounts. This program requires that the annotations of the reference genomes are in SAF format. The script **gff3_to_saf.py**, located at `../RNA-Seq/`, parses the GFF3 lines annotated as CDS, ncRNA, tmRNA (*ssrA*), RNase P (*rnpB*), tRNA, antisense RNA and SRP RNA and outputs them in SAF format:

```sh
../RNA-Seq/gff3_to_saf.py ./RNAseq_C325/C325_maskedIS.gff > ./RNAseq_C325/C325_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_CF13/CF13_maskedIS.gff > ./RNAseq_CF13/CF13_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_H53/H53_maskedIS.gff > ./RNAseq_H53/H53_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K147/K147_maskedIS.gff > ./RNAseq_K147/K147_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_C063/C063p_maskedIS.gff > ./RNAseq_C063/C063p_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K091/K091p_maskedIS.gff > ./RNAseq_K091/K091p_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K141/K141p_maskedIS.gff > ./RNAseq_K141/K141p_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K209/K209p_maskedIS.gff > ./RNAseq_K209/K209p_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K249/K249p_maskedIS.gff > ./RNAseq_K249/K249p_maskedIS.saf
../RNA-Seq/gff3_to_saf.py ./RNAseq_K275/K275p_maskedIS.gff > ./RNAseq_K275/K275p_maskedIS.saf
```

Finally, the script **quant_diffexpr.Rmd** (located at each strain's *RNAseq* directory) was run for each strain  to perform read quantification with **featureCounts v2.14.2**, data exploratory analyses and DE with **DESeq2 v1.40.1**. DE analyses are always performed comparing the pOXA-48-carrying versions of the strains (transconjugants or WT pOXA-48-carrying strains) against the pOXA-48-free strains (WT naive strains or cured from pOXA-48). In general, the percentage of successfully assigned alignments to features was high (>90%), and replicates between conditions separated well by PC1, except in K275/KPN18, possibly reflecting low impact of pOXA-48 carriage in this strain. See the **quant_diffexpr.html** file of each analysis:

* [Read quantification and DE analysis of strain C325](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_C325/quant_diffexpr.html)
* [Read quantification and DE analysis of strain CF13](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_CF13/quant_diffexpr.html)
* [Read quantification and DE analysis of strain H53](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_H53/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K147](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K147/quant_diffexpr.html)
* [Read quantification and DE analysis of strain C063](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_C063/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K091](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K091/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K141](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K141/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K209](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K209/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K249](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K249/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K275](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_experimental_evolution/RNAseq_K275/quant_diffexpr.html)

