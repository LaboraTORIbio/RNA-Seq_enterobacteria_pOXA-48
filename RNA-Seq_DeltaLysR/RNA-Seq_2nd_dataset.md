# Transcriptomic Analyses for the Second RNA-Seq Dataset

## Project setup

The transcriptomes of strains KPN15 and CF13, carrying pOXA-48 or pOXA-48Δ*lysR*, were sequenced (3 biological replicates for each strain and condition). Reads are available at BioProject **PRJNA1071971**:

|strain_replicate|fastq code|
|:----|:----|
|CF13_pOXA-48_1|4328_10_S22|
|CF13_pOXA-48_2|4328_11_S23|
|CF13_pOXA-48_3|4328_12_S24|
|CF13_pOXA-48DlysR_1|4328_16_S28|
|CF13_pOXA-48DlysR_2|4328_17_S29|
|CF13_pOXA-48DlysR_3|4328_18_S30|
|KPN15_pOXA-48_1|4328_1_S13|
|KPN15_pOXA-48_2|4328_2_S14|
|KPN15_pOXA-48_3|4328_3_S15|
|KPN15pOXA-48DlysR_1|4328_7_S19|
|KPN15pOXA-48DlysR_2|4328_8_S20|
|KPN15pOXA-48DlysR_3|4328_9_S21|

Create the project directory from which run all commands:

```sh
mkdir -p RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_DeltaLysR
cd RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_DeltaLysR
```

Here, each strain will have its own directory, denoted with the prefix *RNAseq* (e.g. RNAseq_C325), that will contain BAM mapping files, the Rproject and necessary files for conducting differential expression (DE) analysis in R.


## Read trimming and quality control

RNA-Seq raw reads (**PRJNA1071971**) were located at `./reads_RNAseq/raw_reads/`. Before trimming, we performed an initial analysis of read quality with **FastQC v0.11.9** and **MultiQC v1.11**. In general, all reads showed good quality metrics.

```sh
mkdir reads_RNAseq/raw_reads/fastQC_analysis/
fastqc reads_RNAseq/raw_reads/* -o reads_RNAseq/raw_reads/fastQC_analysis/
multiqc reads_RNAseq/raw_reads/fastQC_analysis/* -o reads_RNAseq/raw_reads/fastQC_analysis/multiqc
```

Reads were trimmed with **Trim Galore v0.6.4**, generating FastQC reports. Note that the base name is the prefix *SeqCoast* (sequencing center) followed by the strain name and replicate number (e.g. SeqCoast_CF13_pOXA-48_1):

```sh
trim_galore --paired --quality 20 --length 50 --fastqc --basename SeqCoast_<strain_name> --output_dir reads_RNAseq/ reads_RNAseq/raw_reads/<fastq1> reads_RNAseq/raw_reads/<fastq2>
```

FastQC reports were combined with **MultiQC v1.11**:

```sh
mkdir reads_RNAseq/fastQC_analysis
mv *html reads_RNAseq/fastQC_analysis/
mv *zip reads_RNAseq/fastQC_analysis/
multiqc reads_RNAseq/fastQC_analysis/* -o reads_RNAseq/fastQC_analysis/multiqc
```

#### General statistics

Samples have >6 million reads. The GC content is 51-53% for *C. freundii* and 54-555% for *K. pneumoniae* strains (note that MGEs and overrepresented sequences could be skewing the GC content). Finally, the percentage of duplicated sequences is 76-82% across samples, which probably represent overrepresented sequences.

### Analysis of overrepresented sequences

Most of the overrepresented sequences are commonly present in more than one sample, counted up to 17 times across all samples:

```sh
# get all the overrepresented sequences from the fastQC output (of trimmed reads) and find the unique ones, with count:
unzip -p reads_RNAseq/fastQC_analysis/\*zip \*/fastqc_data.txt | sed -n '/>>Overrepresented/,/>>END_MODULE/{/>>END_MODULE/!p}' | grep -E "^[A|T|G|C]" | cut -f1 | sort | uniq -c | sort -rn
# and storing the sequences in fasta
unzip -p reads_RNAseq/fastQC_analysis/\*zip \*/fastqc_data.txt | sed -n '/>>Overrepresented/,/>>END_MODULE/{/>>END_MODULE/!p}' | grep -E "^[A|T|G|C]" | cut -f1 | sort | uniq -c | sort -rn | awk '{print $2}' | awk '{printf("%s%s\n", (++num==1 ? "" : ">"num"\n"), $0)}' | sed '1s/^/>1\n/' > reads_RNAseq/overrepresented_sequences.fasta
```

We analyzed to what genetic features the overrepresented sequences correspond. For this, they were blasted (**BLASTn v2.11.0**) against the reference genomes, storing the hits in BED format:

```sh
mkdir blastn_overrepresented_seqs

makeblastdb -in ../../Closed_sequences/CF13.fasta -dbtype nucl
makeblastdb -in ../../Closed_sequences/TC_KPN15.fasta -dbtype nucl

# the blast table is filtered by 100% identity and full alignment length (50), then only the chr, start and end of the match is stored (bed format)
blastn -query reads_RNAseq/overrepresented_sequences.fasta -db ../../Closed_sequences/CF13.fasta -outfmt 6 | awk '{ if ($3 == "100.000" && $4 == "50") { print $2"\t"$9"\t"$10 } }' | awk '{ if ($3-$2 < 0) {print $1"\t"$3"\t"$2} else {print $1"\t"$2"\t"$3} }' > blastn_overrepresented_seqs/CF13.bed
blastn -query reads_RNAseq/overrepresented_sequences.fasta -db ../../Closed_sequences/TC_KPN15.fasta -outfmt 6 | awk '{ if ($3 == "100.000" && $4 == "50") { print $2"\t"$9"\t"$10 } }' | awk '{ if ($3-$2 < 0) {print $1"\t"$3"\t"$2} else {print $1"\t"$2"\t"$3} }' > blastn_overrepresented_seqs/TC_KPN15.bed
```

Next, **bedtools v2.27.1** is used to find the gene annotation matching the overrepresented regions. The number of matches per gene, per strain, is counted:

```sh
bedtools intersect -wa -wb -a blastn_overrepresented_seqs/CF13.bed -b ../../Closed_sequences/CF13.gff | awk '{ if ($6 == "gene") {print $0} }' | cut -d";" -f2 | cut -d"=" -f2 | sort | uniq -c | sort -rn
bedtools intersect -wa -wb -a blastn_overrepresented_seqs/TC_KPN15.bed -b ../../Closed_sequences/TC_KPN15.gff | awk '{ if ($6 == "gene") {print $0} }' | cut -d";" -f2 | cut -d"=" -f2 | sort | uniq -c | sort -rn
```

In all strains, *ssrA* (or tmRNA) is the most overrepresented gene (as with the first RNA-Seq dataset). The rest of overrepresented sequences correspond to 23S and 16S ribosomal RNA, but the percentage of overrepresented sequences is <6% in all samples and rRNA will not be counted at the read count step, so we continued with the analyses.



## Mapping RNA-Seq reads

Since splicing is rare in bacteria, we can use a regular aligner to map the RNA-Seq reads to their respective reference genomes. Here, we used **BWA-MEM v0.7.17**:

```sh
mkdir ./RNAseq_CF13/
mkdir ./RNAseq_KPN15/

bwa index ../../Closed_sequences/CF13.fasta
bwa index ../../Closed_sequences/TC_KPN15.fasta

for fq1 in reads_RNAseq/SeqCoast_CF13*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:22}
	name=${name::-12}
	bwa mem ../../Closed_sequences/CF13.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_CF13/$name".bam"
done

for fq1 in reads_RNAseq/SeqCoast_KPN15*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:22}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN15.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN15/$name".bam"
done

# Mapping stats (>99% mapped reads):
for bam in ./RNAseq_*/*bam
do
	echo "####################################################"
	echo "MAPPING STATS OF $bam"
	samtools flagstats $bam
done
```


## Differential expression analysis

The first step in DE analysis is counting the number of reads mapped to each genomic feature, which is performed here with featureCounts. This program requires that the annotations of the reference genomes are in SAF format. The **SAF** files generated for the first RNA-Seq dataset are reused for these DE analyses (copied to the respective *RNAseq* directories).

The script **quant_diffexpr.Rmd** (located at each strain's *RNAseq* directory) was run for each strain  to perform read quantification with **featureCounts v2.14.2**, data exploratory analyses and DE with **DESeq2 v1.40.1**. DE analyses are performed comparing the pOXA-48-carrying versions of the strains against the pOXA-48Δ*lysR*-carrying strains. In general, the percentage of successfully assigned alignments to features was high (>85%), and replicates between conditions separated well by PC1. See the **quant_diffexpr.html** file of each analysis:

* [Read quantification and DE analysis of strain CF13](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_DeltaLysR/RNAseq_CF13/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN15](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_DeltaLysR/RNAseq_KPN15/quant_diffexpr.html)

The Chr number of pOXA-48 is changed to pOXA-48 in the filtered tables to merge the results with the **DEGs_table_parser_v2.py** script, generating a summary table of DEGs. Only includes CDSs, summarized by "Chr", "Type", "RefSeq", "Gene" and "Product". Other features like tRNAs and CDSs with no RefSeq and annotated as hypothetical proteins are excluded because we cannot differenciate them.

```sh
./DEGs_table_parser_v2.py -d RNAseq_CF13/ RNAseq_KPN15/ -n DEGs_summary
```

Genes putatively controlled by other chromosomal LysRs were searched manually in the DE raw results files and saved to the **first_genes_reg_by_lysRs.txt** files located at each strains *RNAseq* directory. The script **plot_genes_reg_by_lysRs.R** (in the homonymous directory) plots the DE of these genes (Fig 5).


## Variant calling

To assess whether mutations accumulated in the chromosome or other plasmids during growth or construction of pOXA-48-carrying or pOXA-48-free samples affected gene expression, variant calling was performed on the transcribed regions with **Snippy v4.6.0**.

```sh
# CF13
snippy --outdir variants_snippy/CF13_pOXA-48_1 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48_1.bam
snippy --outdir variants_snippy/CF13_pOXA-48_2 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48_2.bam
snippy --outdir variants_snippy/CF13_pOXA-48_3 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48_3.bam
snippy --outdir variants_snippy/CF13_pOXA-48DlysR_1 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48DlysR_1.bam
snippy --outdir variants_snippy/CF13_pOXA-48DlysR_2 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48DlysR_2.bam
snippy --outdir variants_snippy/CF13_pOXA-48DlysR_3 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13_pOXA-48DlysR_3.bam

# KPN15
snippy --outdir variants_snippy/KPN15_pOXA-48_1 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48_1.bam
snippy --outdir variants_snippy/KPN15_pOXA-48_2 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48_2.bam
snippy --outdir variants_snippy/KPN15_pOXA-48_3 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48_3.bam
snippy --outdir variants_snippy/KPN15_pOXA-48DlysR_1 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48DlysR_1.bam
snippy --outdir variants_snippy/KPN15_pOXA-48DlysR_2 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48DlysR_2.bam
snippy --outdir variants_snippy/KPN15_pOXA-48DlysR_3 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/KPN15_pOXA-48DlysR_3.bam
```

To confirm or discard suspected SNPs in plasmid pOXA-48, the RNA-Seq reads were mapped to the sequence of the pOXA-48_K8 variant (accession number MT441554):

```sh
# CF13
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48_1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48_1_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48_1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48_2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48_2_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48_2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48_3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48_3_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48_3_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48DlysR_1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_1_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48DlysR_2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_2_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13_pOXA-48DlysR_3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_3_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_CF13_pOXA-48DlysR_3_val_2.fq.gz

# KPN15
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48_1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48_1_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48_1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48_2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48_2_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48_2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48_3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48_3_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48_3_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48DlysR_1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_1_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48DlysR_2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_2_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_KPN15_pOXA-48DlysR_3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_3_val_1.fq.gz --R2 reads_RNAseq/SeqCoast_KPN15_pOXA-48DlysR_3_val_2.fq.gz
```
