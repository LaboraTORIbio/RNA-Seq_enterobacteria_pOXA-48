# Closing Reference Genomes

## Project setup

The transcriptomes of 12 enterobacterial strains, with and without plasmid pOXA-48, were sequenced:
* **5 clinical strains naturally carrying pOXA-48:** 1 *E. coli* (C325), 1 *C. freundii*  (CF13) and 3 *K. pneumoniae* (H53, J57 and K147) strains. pOXA-48 was cured from these strains; cured strains are denoted with the suffix *c1*. The reference genomes of these strains (carrying pOXA-48) were already closed for [Fernandez-Calvet *et al.* 2023](https://doi.org/10.1099/mic.0.001369) (check [GitHub](https://github.com/LaboraTORIbio/CRISPR_cured_pOXA-48)) and are available under the BioProject **PRJNA626430**.
* **6 clinical strains naive for pOXA-48 but ecologically compatible with it:** 1 *E. coli* (PF_EC10/C063) and 5 *K. pneumoniae* (PF_KPN04/K091, PF_KPN07/K141, PF_KPN10/K209, PF_KPN15/K249 and PF_KPN18/K275) strains. pOXA-48 was introduced into these strains by conjugation; transconjugants are denoted with the prefix *TC* instead of *PF*. The genomes of the transconjugant strains will be closed here.
* **1 laboratory *E. coli* MG1655 strain.** pOXA-48 was introduced into this strain by conjugation; transconjugant is denoted with the suffix *p*. The genome of the transconjugant strain will be closed here.

Create the project directory from which run all commands:

```sh
mkdir -p RNA-Seq_enterobacteria_pOXA-48/Genome_assemblies
cd RNA-Seq_enterobacteria_pOXA-48/Genome_assemblies
```

## Read trimming and genome assembly

The raw Illumina reads of the 6 TC clinical strains (**PRJNA641166**; strain nomenclature equivalence: TC_EC10 = Ec010, TC_KPN04 = Kpn04, TC_KPN07 = Kpn08, TC_KPN10 = Kpn13, TC_KPN15 = Kpn20, TC_KPN18 = Kpn23), were located at `../../Reads_PF_TC/raw_reads_originals/`. Illumina reads for the *E. coli* MG1655 strain (**PRJNA1071971**) were located at `../../Reads_others/`. Long Nanopore reads generated at MiGS (**PRJNA1071971**) were located at `../../Reads_long/` (note the suffix *nanopore*).

Illumina reads (generated at the Wellcome Trust Centre for Human Genetics) were trimmed with **Trim Galore v0.6.4** (Cutadapt v2.8). Note that the base name is the prefix *WTCHG* (sequencing center) followed by the strain name (e.g. WTCHG_TC_EC10). The reads of MG1655 came already trimmed by the MicrobesNG sequencing center (prefix *MicNG*).

```sh
trim_galore --quality 20 --length 50 --fastqc --basename WTCHG_<strain_name> --output_dir ../../Reads_PF_TC/trimmed_originals --paired ../../Reads_PF_TC/raw_reads_originals/<fq1> ../../Reads_PF_TC/raw_reads_originals/<fq2>
```

Hybrid assemblies were generated with **Unicycler v0.4.9**:

```sh
# TC_EC10
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_EC10_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_EC10_val_2.fq.gz -l ../../Reads_long/TC_EC10_nanopore.fastq.gz -o ./assemblies_unicycler/TC_EC10
# TC_KPN04 (discarded assembly)
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN04_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN04_val_2.fq.gz -l ../../Reads_long/TC_KPN04_nanopore.fastq.gz -o ./assemblies_unicycler/TC_KPN04
# TC_KPN07 (discarded assembly)
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN07_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN07_val_2.fq.gz -l ../../Reads_long/TC_KPN07_nanopore.fastq.gz -o ./assemblies_unicycler/TC_KPN07
# TC_KPN10
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN10_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN10_val_2.fq.gz -l ../../Reads_long/TC_KPN10_nanopore.fastq.gz -o ./assemblies_unicycler/TC_KPN10
# TC_KPN15
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN15_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN15_val_2.fq.gz -l ../../Reads_long/TC_KPN15_nanopore.fastq.gz -o ./assemblies_unicycler/TC_KPN15
# TC_KPN18 (discarded assembly)
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN18_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_originals/WTCHG_TC_KPN18_val_2.fq.gz -l ../../Reads_long/TC_KPN18_nanopore.fastq.gz -o ./assemblies_unicycler/TC_KPN18
# MG1655p
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_others/MicNG_MG1655p_1_trimmed.fastq.gz -2 ../../Reads_others/MicNG_MG1655p_2_trimmed.fastq.gz -l ../../Reads_long/MG1655p_nanopore.fastq.gz -o ./assemblies_unicycler/MG1655p
```

The assemblies of strains TC_KPN04, TC_KPN07 and TC_KPN18 were fragmented, possibly because short and long reads came from different DNA extractions. Thus, these strains were sequenced again from a new DNA extraction with Illumina at MiGS and with minION at the lab (reads available at **PRJNA1071971**).

Illumina reads were trimmed with **Trim Galore v0.6.4** as before; trimmed reads were outputted to `../../Reads_PF_TC/trimmed_MiGS/` (note the prefix *MiGS* on the trimmed reads). The minION reads were placed at `../../Reads_long/` and identified with the suffix *minion*. Then, hybrid assemblies for these three strains were obtained with **Unicycler v0.4.9**:

```sh
# TC_KPN04
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_2.fq.gz -l ../../Reads_long/TC_KPN04_minion.fastq.gz -o ./assemblies_unicycler/TC_KPN04_minion
# TC_KPN07
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN07_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN07_val_2.fq.gz -l ../../Reads_long/TC_KPN07_minion.fastq.gz -o ./assemblies_unicycler/TC_KPN07_minion
# TC_KPN18
unicycler --spades_path ~/Software/SPAdes-3.13.0-Linux/bin/spades.py -1 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_1.fq.gz -2 ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_2.fq.gz -l ../../Reads_long/TC_KPN18_minion.fastq.gz -o ./assemblies_unicycler/TC_KPN18_minion
```

These assemblies were not fragmented, but plasmid 2 of strains TC_KPN04 and TC_KPN18 could not be closed. We attempted to close these plasmids by assembling only the long reads with Flye v2.9. First, we compared the quality metrics of the Nanopore (MiGS) and minION reads:

```sh
NanoStat --fastq ../../Reads_long/TC_KPN04_minion.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN04_minion
NanoStat --fastq ../../Reads_long/TC_KPN04_nanopore.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN04_nanopore
NanoStat --fastq ../../Reads_long/TC_KPN07_minion.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN07_minion
NanoStat --fastq ../../Reads_long/TC_KPN07_nanopore.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN07_nanopore
NanoStat --fastq ../../Reads_long/TC_KPN18_minion.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN18_minion
NanoStat --fastq ../../Reads_long/TC_KPN18_nanopore.fastq.gz --outdir ../../Reads_long/statreports -n TC_KPN18_nanopore
grep "Total bases:" ../../Reads_long/statreports/*
grep "Mean read quality" ../../Reads_long/statreports/*
grep "Median read quality" ../../Reads_long/statreports/*
```

Since Nanopore reads (MiGS) had higher coverage and mean/median read quality compared to minION reads, they were chosen for the assemblies. First, Nanopore reads were filtered with **filtlong v0.2.1** to obtain a subset of high identity reads with a minimum read depth of 85x. Filtered long reads were assembled with **Flye v2.9**, and resulting assemblies were polished with **Medaka v1.4.3** and several rounds of **Pilon v1.24**, mapping the trimmed Illumina reads (MiGS), until no more changes were made. Finally, polished assemblies were recircularized with **circlator v1.5.5**. This way, plasmid 2 was closed in both strains.

```sh
### TC_KPN04

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/TC_KPN04_nanopore.fastq.gz | gzip > ../../Reads_long/TC_KPN04_nanopore_filt.fastq.gz
flye --nano-raw ../../Reads_long/TC_KPN04_nanopore_filt.fastq.gz --out-dir ./assemblies_flye/TC_KPN04
medaka_consensus -i ../../Reads_long/TC_KPN04_nanopore_filt.fastq.gz -d ./assemblies_flye/TC_KPN04/assembly.fasta -o ./assemblies_flye/TC_KPN04/medaka
cp ./assemblies_flye/TC_KPN04/medaka/consensus.fasta ./assemblies_flye/TC_KPN04_medaka.fasta

# 4 rounds of pilon + circlator
# round 1
bwa index ./assemblies_flye/TC_KPN04_medaka.fasta
bwa mem ./assemblies_flye/TC_KPN04_medaka.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN04_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN04_medaka_align.bam
pilon --changes --genome ./assemblies_flye/TC_KPN04_medaka.fasta --frags ./assemblies_flye/TC_KPN04_medaka_align.bam --output TC_KPN04_pilon --outdir ./assemblies_flye/pilon
# round 2
bwa index ./assemblies_flye/pilon/TC_KPN04_pilon.fasta
bwa mem ./assemblies_flye/pilon/TC_KPN04_pilon.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN04_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN04_medaka_align.bam
pilon --changes --genome ./assemblies_flye/pilon/TC_KPN04_pilon.fasta --frags ./assemblies_flye/TC_KPN04_medaka_align.bam --output TC_KPN04_pilon_2 --outdir ./assemblies_flye/pilon
# round 3
bwa index ./assemblies_flye/pilon/TC_KPN04_pilon_2.fasta
bwa mem ./assemblies_flye/pilon/TC_KPN04_pilon_2.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN04_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN04_medaka_align.bam
pilon --changes --genome ./assemblies_flye/pilon/TC_KPN04_pilon_2.fasta --frags ./assemblies_flye/TC_KPN04_medaka_align.bam --output TC_KPN04_pilon_3 --outdir ./assemblies_flye/pilon
# round 4
bwa index ./assemblies_flye/pilon/TC_KPN04_pilon_3.fasta
bwa mem ./assemblies_flye/pilon/TC_KPN04_pilon_3.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN04_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN04_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN04_medaka_align.bam
pilon --changes --genome ./assemblies_flye/pilon/TC_KPN04_pilon_3.fasta --frags ./assemblies_flye/TC_KPN04_medaka_align.bam --output TC_KPN04_pilon_4 --outdir ./assemblies_flye/pilon

circlator fixstart ./assemblies_flye/pilon/TC_KPN04_pilon_4.fasta TC_KPN04_circlator

### TC_KPN18

filtlong --min_length 1000 --keep_percent 90 --target_bases 500000000 ../../Reads_long/TC_KPN18_nanopore.fastq.gz | gzip > ../../Reads_long/TC_KPN18_nanopore_filt.fastq.gz
flye --nano-raw ../../Reads_long/TC_KPN18_nanopore_filt.fastq.gz --out-dir ./assemblies_flye/TC_KPN18
medaka_consensus -i ../../Reads_long/TC_KPN18_nanopore_filt.fastq.gz -d ./assemblies_flye/TC_KPN18/assembly.fasta -o ./assemblies_flye/TC_KPN18/medaka
cp ./assemblies_flye/TC_KPN18/medaka/consensus.fasta ./assemblies_flye/TC_KPN18_medaka.fasta

# 2 rounds of pilon + circlator
# round 1
bwa index ./assemblies_flye/TC_KPN18_medaka.fasta
bwa mem ./assemblies_flye/TC_KPN18_medaka.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN18_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN18_medaka_align.bam
pilon --changes --genome ./assemblies_flye/TC_KPN18_medaka.fasta --frags ./assemblies_flye/TC_KPN18_medaka_align.bam --output TC_KPN18_pilon --outdir ./assemblies_flye/pilon
# round 2
bwa index ./assemblies_flye/pilon/TC_KPN18_pilon.fasta
bwa mem ./assemblies_flye/pilon/TC_KPN18_pilon.fasta ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_1.fq.gz ../../Reads_PF_TC/trimmed_MiGS/MiGS_TC_KPN18_val_2.fq.gz | samtools sort -o ./assemblies_flye/TC_KPN18_medaka_align.bam
samtools index ./assemblies_flye/TC_KPN18_medaka_align.bam
pilon --changes --genome ./assemblies_flye/pilon/TC_KPN18_pilon.fasta --frags ./assemblies_flye/TC_KPN18_medaka_align.bam --output TC_KPN18_pilon_2 --outdir ./assemblies_flye/pilon

circlator fixstart ./assemblies_flye/pilon/TC_KPN18_pilon_2.fasta TC_KPN18_circlator
```


## Final selection of reference genomes

For strains TC_EC10, TC_KPN10, TC_KPN15 and MG1655p, we selected the first Unicycler assemblies. For strains TC_KPN04, TC_KPN07 and TC_KPN18, we selected the second Unicycler assemblies, but replacing the sequence of plasmid 2 by the sequence of the closed plasmid obtained from the Flye assemblies. These final closed genomes were then coppied to `../../Closed_sequences/` and annotated with **PGAP v2021-07-01.build5508**.
