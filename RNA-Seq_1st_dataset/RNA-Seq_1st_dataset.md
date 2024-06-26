# Transcriptomic Analyses for the First RNA-Seq Dataset

## Project setup

Before starting, check how reference genomes were assembled and closed ([closing_reference_genomes.md](RNA-Seq_enterobacteria_pOXA-48/Genome_assemblies/closing_reference_genomes.md)) and how RNA-Seq reads were processed ([read_quality_control.md](RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/read_quality_control.md)).

Create the project directory from which run all commands:

```sh
mkdir -p RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset
cd RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset
```

Here, each strain will have its own directory, denoted with the prefix *RNAseq* (e.g. RNAseq_C325), that will contain BAM mapping files, the Rproject and necessary files for conducting differential expression (DE) analysis in R.


## Mapping RNA-Seq reads

Since splicing is rare in bacteria, we can use a regular aligner to map the RNA-Seq reads to their respective reference genomes. Here, we used **BWA-MEM v0.7.17**:

```sh
mkdir ./RNAseq_C325/
mkdir ./RNAseq_CF13/
mkdir ./RNAseq_H53/
mkdir ./RNAseq_J57/
mkdir ./RNAseq_K147/
mkdir ./RNAseq_EC10/
mkdir ./RNAseq_KPN04/
mkdir ./RNAseq_KPN07/
mkdir ./RNAseq_KPN10/
mkdir ./RNAseq_KPN15/
mkdir ./RNAseq_KPN18/
mkdir ./RNAseq_MG1655/

bwa index ../../Closed_sequences/C325.fasta
bwa index ../../Closed_sequences/CF13.fasta
bwa index ../../Closed_sequences/H53.fasta
bwa index ../../Closed_sequences/J57.fasta
bwa index ../../Closed_sequences/K147.fasta
bwa index ../../Closed_sequences/TC_EC10.fasta
bwa index ../../Closed_sequences/TC_KPN04.fasta
bwa index ../../Closed_sequences/TC_KPN07.fasta
bwa index ../../Closed_sequences/TC_KPN10.fasta
bwa index ../../Closed_sequences/TC_KPN15.fasta
bwa index ../../Closed_sequences/TC_KPN18.fasta
bwa index ../../Closed_sequences/MG1655p.fasta

for fq1 in reads_RNAseq/WTCHG_C325*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/C325.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_C325/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_CF13*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/CF13.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_CF13/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_H53*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/H53.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_H53/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_J57*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/J57.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_J57/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_K147*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/K147.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_K147/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*EC10*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_EC10.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_EC10/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*KPN04*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN04.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN04/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*KPN07*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN07.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN07/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*KPN10*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN10.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN10/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*KPN15*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN15.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN15/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_*KPN18*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/TC_KPN18.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_KPN18/$name".bam"
done

for fq1 in reads_RNAseq/WTCHG_MG1655*_val_1.fq.gz
do
	fq2=${fq1%%1.fq.gz}"2.fq.gz"
	name=${fq1:19}
	name=${name::-12}
	bwa mem ../../Closed_sequences/MG1655p.fasta $fq1 $fq2 | samtools sort -o ./RNAseq_MG1655/$name".bam"
done
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

The first step in DE analysis is counting the number of reads mapped to each genomic feature, which is performed here with featureCounts. This program requires that the annotations of the reference genomes are in SAF format. The script **gff3_to_saf.py**, located at this directory, parses the GFF3 lines annotated as CDS, ncRNA, tmRNA (*ssrA*), RNase P (*rnpB*), tRNA, antisense RNA and SRP RNA and outputs them in SAF format:

```sh
./gff3_to_saf.py ../../Closed_sequences/C325.gff > ./RNAseq_C325/C325.saf
./gff3_to_saf.py ../../Closed_sequences/CF13.gff > ./RNAseq_CF13/CF13.saf
./gff3_to_saf.py ../../Closed_sequences/H53.gff > ./RNAseq_H53/H53.saf
./gff3_to_saf.py ../../Closed_sequences/J57.gff > ./RNAseq_J57/J57.saf
./gff3_to_saf.py ../../Closed_sequences/K147.gff > ./RNAseq_K147/K147.saf
./gff3_to_saf.py ../../Closed_sequences/TC_EC10.gff > ./RNAseq_EC10/TC_EC10.saf
./gff3_to_saf.py ../../Closed_sequences/TC_KPN04.gff > ./RNAseq_KPN04/TC_KPN04.saf
./gff3_to_saf.py ../../Closed_sequences/TC_KPN07.gff > ./RNAseq_KPN07/TC_KPN07.saf
./gff3_to_saf.py ../../Closed_sequences/TC_KPN10.gff > ./RNAseq_KPN10/TC_KPN10.saf
./gff3_to_saf.py ../../Closed_sequences/TC_KPN15.gff > ./RNAseq_KPN15/TC_KPN15.saf
./gff3_to_saf.py ../../Closed_sequences/TC_KPN18.gff > ./RNAseq_KPN18/TC_KPN18.saf
./gff3_to_saf.py ../../Closed_sequences/MG1655p.gff > ./RNAseq_MG1655/MG1655p.saf
```

Finally, the script **quant_diffexpr.Rmd** (located at each strain's *RNAseq* directory) was run for each strain  to perform read quantification with **featureCounts v2.14.2**, data exploratory analyses and DE with **DESeq2 v1.40.1**. DE analyses are always performed comparing the pOXA-48-carrying versions of the strains (transconjugants or WT pOXA-48-carrying strains) against the pOXA-48-free strains (WT naive strains or cured from pOXA-48). In general, the percentage of successfully assigned alignments to features was high (>85%), and replicates between conditions separated well by PC1, except in KPN18, possibly reflecting low impact of pOXA-48 carriage in this strain. See the **quant_diffexpr.html** file of each analysis:

* [Read quantification and DE analysis of strain C325](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_C325/quant_diffexpr.html)
* [Read quantification and DE analysis of strain CF13](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_CF13/quant_diffexpr.html)
* [Read quantification and DE analysis of strain H53](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_H53/quant_diffexpr.html)
* [Read quantification and DE analysis of strain J57](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_J57/quant_diffexpr.html)
* [Read quantification and DE analysis of strain K147](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_K147/quant_diffexpr.html)
* [Read quantification and DE analysis of strain EC10](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_EC10/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN04](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_KPN04/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN07](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_KPN07/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN10](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_KPN10/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN15](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_KPN15/quant_diffexpr.html)
* [Read quantification and DE analysis of strain KPN18](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_KPN18/quant_diffexpr.html)
* [Read quantification and DE analysis of strain MG1655](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/RNAseq_MG1655/quant_diffexpr.html)

The DE results of chromosomes and plasmids (excluding pOXA-48) are separated in different .tsv files for following analyses:

```sh
# raw
cat RNAseq_C325/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_C325/DE_results_raw_chromosome.tsv
cat RNAseq_C325/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 4' > RNAseq_C325/DE_results_raw_plasmids.tsv
cat RNAseq_CF13/DE_results_raw.tsv | awk 'NR == 1 || $12 == 1' > RNAseq_CF13/DE_results_raw_chromosome.tsv
cat RNAseq_CF13/DE_results_raw.tsv | awk 'NR == 1 || $12 != 1 && $12 != 2' > RNAseq_CF13/DE_results_raw_plasmids.tsv
cat RNAseq_H53/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_H53/DE_results_raw_chromosome.tsv
cat RNAseq_H53/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 3' > RNAseq_H53/DE_results_raw_plasmids.tsv
cat RNAseq_J57/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_J57/DE_results_raw_chromosome.tsv
cat RNAseq_J57/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 4' > RNAseq_J57/DE_results_raw_plasmids.tsv
cat RNAseq_K147/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_K147/DE_results_raw_chromosome.tsv
cat RNAseq_K147/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 3' > RNAseq_K147/DE_results_raw_plasmids.tsv
cat RNAseq_EC10/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_EC10/DE_results_raw_chromosome.tsv
cat RNAseq_EC10/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 3' > RNAseq_EC10/DE_results_raw_plasmids.tsv
cat RNAseq_KPN04/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_KPN04/DE_results_raw_chromosome.tsv
cat RNAseq_KPN04/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 3' > RNAseq_KPN04/DE_results_raw_plasmids.tsv
cat RNAseq_KPN07/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_KPN07/DE_results_raw_chromosome.tsv
cat RNAseq_KPN07/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 4' > RNAseq_KPN07/DE_results_raw_plasmids.tsv
cat RNAseq_KPN10/DE_results_raw.tsv | awk 'NR == 1 || $11 == 1' > RNAseq_KPN10/DE_results_raw_chromosome.tsv
cat RNAseq_KPN10/DE_results_raw.tsv | awk 'NR == 1 || $11 != 1 && $11 != 4' > RNAseq_KPN10/DE_results_raw_plasmids.tsv
cat RNAseq_KPN15/DE_results_raw.tsv | awk 'NR == 1 || $12 == 1' > RNAseq_KPN15/DE_results_raw_chromosome.tsv
cat RNAseq_KPN15/DE_results_raw.tsv | awk 'NR == 1 || $12 != 1 && $12 != 4' > RNAseq_KPN15/DE_results_raw_plasmids.tsv
cat RNAseq_KPN18/DE_results_raw.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_KPN18/DE_results_raw_chromosome.tsv
cat RNAseq_KPN18/DE_results_raw.tsv | awk 'NR == 1 || $13 != 1 && $13 != 4' > RNAseq_KPN18/DE_results_raw_plasmids.tsv
cat RNAseq_MG1655/DE_results_raw.tsv | awk 'NR == 1 || $12 == 1' > RNAseq_MG1655/DE_results_raw_chromosome.tsv
# filtered
cat RNAseq_C325/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_C325/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_C325/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 4' > RNAseq_C325/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_CF13/DE_results_filtered_padj.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_CF13/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_CF13/DE_results_filtered_padj.tsv | awk 'NR == 1 || $13 != 1 && $13 != 2' > RNAseq_CF13/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_H53/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_H53/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_H53/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 3' > RNAseq_H53/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_J57/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_J57/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_J57/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 4' > RNAseq_J57/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_K147/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_K147/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_K147/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 3' > RNAseq_K147/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_EC10/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_EC10/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_EC10/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 3' > RNAseq_EC10/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_KPN04/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_KPN04/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_KPN04/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 3' > RNAseq_KPN04/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_KPN07/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_KPN07/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_KPN07/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 4' > RNAseq_KPN07/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_KPN10/DE_results_filtered_padj.tsv | awk 'NR == 1 || $12 == 1' > RNAseq_KPN10/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_KPN10/DE_results_filtered_padj.tsv | awk 'NR == 1 || $12 != 1 && $12 != 4' > RNAseq_KPN10/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_KPN15/DE_results_filtered_padj.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_KPN15/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_KPN15/DE_results_filtered_padj.tsv | awk 'NR == 1 || $13 != 1 && $13 != 4' > RNAseq_KPN15/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_KPN18/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 == 1' > RNAseq_KPN18/DE_results_filtered_chromosome_padj.tsv
cat RNAseq_KPN18/DE_results_filtered_padj.tsv | awk 'NR == 1 || $14 != 1 && $14 != 4' > RNAseq_KPN18/DE_results_filtered_plasmids_padj.tsv
cat RNAseq_MG1655/DE_results_filtered_padj.tsv | awk 'NR == 1 || $13 == 1' > RNAseq_MG1655/DE_results_filtered_chromosome_padj.tsv

```


## Exploring DE results

Number of features with DE data (not significant, which excludes pOXA-48 genes and 0 count genes):

```sh
echo "C325: of" $(expr $(cat RNAseq_C325/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_C325/DE_results_raw.tsv | awk '$13 != 4' | wc -l) - 1) "have DE data"
# C325: of 5479 features, 5424 have DE data
echo "CF13: of" $(expr $(cat RNAseq_CF13/*saf | awk '$2 != 2' |  wc -l) - 1) "features," $(expr $(cat RNAseq_CF13/DE_results_raw.tsv | awk '$12 != 2' | wc -l) - 1) "have DE data"
# CF13: of 4928 features, 4908 have DE data
echo "H53: of" $(expr $(cat RNAseq_H53/*saf | awk '$2 != 3' |  wc -l) - 1) "features," $(expr $(cat RNAseq_H53/DE_results_raw.tsv | awk '$13 != 3' | wc -l) - 1) "have DE data"
# H53: of 5336 features, 5302 have DE data
echo "J57: of" $(expr $(cat RNAseq_J57/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_J57/DE_results_raw.tsv | awk '$13 != 4' | wc -l) - 1) "have DE data"
# J57: of 5710 features, 5685 have DE data
echo "K147: of" $(expr $(cat RNAseq_K147/*saf | awk '$2 != 3' |  wc -l) - 1) "features," $(expr $(cat RNAseq_K147/DE_results_raw.tsv | awk '$13 != 3' | wc -l) - 1) "have DE data"
# K147: of 5283 features, 5248 have DE data
echo "EC10: of" $(expr $(cat RNAseq_EC10/*saf | awk '$2 != 3' |  wc -l) - 1) "features," $(expr $(cat RNAseq_EC10/DE_results_raw.tsv | awk '$13 != 3' | wc -l) - 1) "have DE data"
# EC10: of 4999 features, 4963 have DE data
echo "KPN04: of" $(expr $(cat RNAseq_KPN04/*saf | awk '$2 != 3' |  wc -l) - 1) "features," $(expr $(cat RNAseq_KPN04/DE_results_raw.tsv | awk '$13 != 3' | wc -l) - 1) "have DE data"
# KPN04: of 5395 features, 5373 have DE data
echo "KPN07: of" $(expr $(cat RNAseq_KPN07/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_KPN07/DE_results_raw.tsv | awk '$13 != 4' | wc -l) - 1) "have DE data"
# KPN07: of 5526 features, 5489 have DE data
echo "KPN10: of" $(expr $(cat RNAseq_KPN10/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_KPN10/DE_results_raw.tsv | awk '$11 != 4' | wc -l) - 1) "have DE data"
# KPN10: of 5521 features, 5471 have DE data
echo "KPN15: of" $(expr $(cat RNAseq_KPN15/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_KPN15/DE_results_raw.tsv | awk '$12 != 4' | wc -l) - 1) "have DE data"
# KPN15: of 5532 features, 5483 have DE data
echo "KPN18: of" $(expr $(cat RNAseq_KPN18/*saf | awk '$2 != 4' |  wc -l) - 1) "features," $(expr $(cat RNAseq_KPN18/DE_results_raw.tsv | awk '$13 != 4' | wc -l) - 1) "have DE data"
# KPN18: of 5635 features, 5629 have DE data
echo "MG1655: of" $(expr $(cat RNAseq_MG1655/*saf | awk '$2 != 2' |  wc -l) - 1) "features," $(expr $(cat RNAseq_MG1655/DE_results_raw.tsv | awk '$12 != 2' | wc -l) - 1) "have DE data"
# MG1655: of 4474 features, 4413 have DE data
```

Percentage of significant chromosomal DEGs:


```sh
echo "C325:" $(expr $(expr $(cat RNAseq_C325/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_C325/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# C325: 15 % of chromosomal genes are DE
echo "CF13:" $(expr $(expr $(cat RNAseq_CF13/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_CF13/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# CF13: 2 % of chromosomal genes are DE
echo "H53:" $(expr $(expr $(cat RNAseq_H53/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_H53/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# H53: 0 % of chromosomal genes are DE
echo "J57:" $(expr $(expr $(cat RNAseq_J57/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_J57/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# J57: 2 % of chromosomal genes are DE
echo "K147:" $(expr $(expr $(cat RNAseq_K147/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_K147/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# K147: 2 % of chromosomal genes are DE
echo "EC10:" $(expr $(expr $(cat RNAseq_EC10/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_EC10/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# EC10: 23 % of chromosomal genes are DE
echo "KPN04:" $(expr $(expr $(cat RNAseq_KPN04/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_KPN04/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# KPN04: 6 % of chromosomal genes are DE
echo "KPN07:" $(expr $(expr $(cat RNAseq_KPN07/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_KPN07/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# KPN07: 37 % of chromosomal genes are DE
echo "KPN10:" $(expr $(expr $(cat RNAseq_KPN10/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_KPN10/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# KPN10: 21 % of chromosomal genes are DE
echo "KPN15:" $(expr $(expr $(cat RNAseq_KPN15/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_KPN15/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# KPN15: 2 % of chromosomal genes are DE
echo "KPN18:" $(expr $(expr $(cat RNAseq_KPN18/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_KPN18/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# KPN18: 2 % of chromosomal genes are DE
echo "MG1655:" $(expr $(expr $(cat RNAseq_MG1655/DE_results_filtered_chromosome_padj.tsv | wc -l) - 1) \* 100 / $(expr $(cat RNAseq_MG1655/*saf | awk '$2 == 1' | wc -l) - 1)) "% of chromosomal genes are DE"
# MG1655: 1 % of chromosomal genes are DE
```
Number of total DEGs by strain:

```sh
for file in RNAseq_*/DE_results_filtered_chromosome_padj.tsv; do echo ${file::-40}":" $(expr $(cat $file | wc -l) - 1); done
# RNAseq_C325: 820
# RNAseq_CF13: 127
# RNAseq_EC10: 1145
# RNAseq_H53: 6
# RNAseq_J57: 136
# RNAseq_K147: 132
# RNAseq_KPN04: 330
# RNAseq_KPN07: 1919
# RNAseq_KPN10: 1076
# RNAseq_KPN15: 116
# RNAseq_KPN18: 142
# RNAseq_MG1655: 63
```


## Retrieving GO annotations

The RefSeq IDs of genes with DE data are retrieved and mapped to **UniProtKB IDs** using the [UniProt database](https://www.uniprot.org/id-mapping). Then, the **GO annotations** are downloaded from the UniProt database, filtering columns to include only the From, Entry, Gene Ontology (biological process), Gene Ontology (cellular component), Gene Ontology (molecular function) and Gene Ontology (GO) columns (May 2023). Tables are saved as **refseq2uniprot_\<strain>_\<chr>.tsv** in the `GSEA/refseq2uniprot/` directory. Since a RefSeq ID can be mapped to multiple UniProt IDs, that generally share the same GO annotations, only one entry is left.

```sh
mkdir -p GSEA/refseq_IDs/

# Retrieving RefSeq IDs:
# chromosome
awk -F"\t" '{print $(NF-2)}' RNAseq_C325/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/C325_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_CF13/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/CF13_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_H53/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/H53_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_J57/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/J57_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_K147/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/K147_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_EC10/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_EC10_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN04/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN04_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN07/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN07_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN10/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN10_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN15/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN15_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN18/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN18_chromosome.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_MG1655/DE_results_raw_chromosome.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/MG1655_chromosome.txt
# plasmids
awk -F"\t" '{print $(NF-2)}' RNAseq_C325/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/C325_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_CF13/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/CF13_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_H53/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/H53_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_J57/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/J57_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_K147/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/K147_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_EC10/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_EC10_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN04/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN04_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN07/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN07_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN10/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN10_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN15/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN15_plasmids.txt
awk -F"\t" '{print $(NF-2)}' RNAseq_KPN18/DE_results_raw_plasmids.tsv | sed 's/"//g' | grep -Ev "RefSeq|-" | sort | uniq | sed 's/\.[0-9]//g' > GSEA/refseq_IDs/TC_KPN18_plasmids.txt

mkdir GSEA/refseq2uniprot/

# Removing duplicate entries:
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_C325_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_CF13_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_H53_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_J57_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_K147_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_EC10_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN04_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN07_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN10_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN15_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN18_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_MG1655_chromosome.tsv > GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_C325_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_CF13_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_H53_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_J57_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_K147_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_EC10_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN04_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN07_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN10_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN15_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv
awk '!seen[$1]++' GSEA/refseq2uniprot/refseq2uniprot_KPN18_plasmids.tsv > GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv
```

More than 65% and 63% of RefSeq IDs and DE features were annotated with GO terms (chromosome+plasmids), respectively:

```sh
echo "C325: of" $(expr $(cat RNAseq_C325/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/C325_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_C325_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_C325_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_C325_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/C325_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_C325_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_C325/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# C325: of 5425 DE features, 5027 have a RefSeq ID, of which 4488 were mapped to UniProt IDs. Finally, 3640 ( 72 % ) RefSeq IDs ( 67 % DE features ) can be annotated with GO terms

echo "CF13: of" $(expr $(cat RNAseq_CF13/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/CF13_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_CF13_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_CF13_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_CF13_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/CF13_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_CF13_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_CF13/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# CF13: of 4909 DE features, 4701 have a RefSeq ID, of which 3658 were mapped to UniProt IDs. Finally, 3099 ( 65 % ) RefSeq IDs ( 63 % DE features ) can be annotated with GO terms

echo "H53: of" $(expr $(cat RNAseq_H53/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/H53_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_H53_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_H53_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_H53_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/H53_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_H53_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_H53/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# H53: of 5303 DE features, 5083 have a RefSeq ID, of which 4614 were mapped to UniProt IDs. Finally, 3718 ( 73 % ) RefSeq IDs ( 70 % DE features ) can be annotated with GO terms

echo "J57: of" $(expr $(cat RNAseq_J57/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/J57_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_J57_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_J57_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_J57_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/J57_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_J57_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_J57/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# J57: of 5686 DE features, 5357 have a RefSeq ID, of which 4668 were mapped to UniProt IDs. Finally, 3793 ( 70 % ) RefSeq IDs ( 66 % DE features ) can be annotated with GO terms

echo "K147: of" $(expr $(cat RNAseq_K147/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/K147_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_K147_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_K147_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_K147_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/K147_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_K147_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_K147/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# K147: of 5249 DE features, 5050 have a RefSeq ID, of which 4729 were mapped to UniProt IDs. Finally, 3769 ( 74 % ) RefSeq IDs ( 71 % DE features ) can be annotated with GO terms

echo "EC10: of" $(expr $(cat RNAseq_EC10/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_EC10_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_EC10_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_EC10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_EC10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_EC10_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_EC10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_EC10/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# EC10: of 4964 DE features, 4757 have a RefSeq ID, of which 4286 were mapped to UniProt IDs. Finally, 3585 ( 75 % ) RefSeq IDs ( 72 % DE features ) can be annotated with GO terms

)echo "KPN04: of" $(expr $(cat RNAseq_KPN04/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_KPN04_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_KPN04_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN04_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN04_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_KPN04_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN04_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_KPN04/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# KPN04: of 5374 DE features, 5078 have a RefSeq ID, of which 4634 were mapped to UniProt IDs. Finally, 3732 ( 73 % ) RefSeq IDs ( 69 % DE features ) can be annotated with GO terms

echo "KPN07: of" $(expr $(cat RNAseq_KPN07/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_KPN07_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_KPN07_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN07_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN07_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_KPN07_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN07_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_KPN07/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# KPN07: of 5490 DE features, 5213 have a RefSeq ID, of which 4724 were mapped to UniProt IDs. Finally, 3782 ( 72 % ) RefSeq IDs ( 68 % DE features ) can be annotated with GO terms

echo "KPN10: of" $(expr $(cat RNAseq_KPN10/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_KPN10_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_KPN10_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_KPN10_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN10_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_KPN10/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# KPN10: of 5472 DE features, 5198 have a RefSeq ID, of which 4713 were mapped to UniProt IDs. Finally, 3779 ( 72 % ) RefSeq IDs ( 69 % DE features ) can be annotated with GO terms

echo "KPN15: of" $(expr $(cat RNAseq_KPN15/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_KPN15_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_KPN15_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN15_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN15_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_KPN15_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN15_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_KPN15/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# KPN15: of 5484 DE features, 5237 have a RefSeq ID, of which 4872 were mapped to UniProt IDs. Finally, 3860 ( 73 % ) RefSeq IDs ( 70 % DE features ) can be annotated with GO terms

echo "KPN18: of" $(expr $(cat RNAseq_KPN18/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/TC_KPN18_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_KPN18_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN18_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN18_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/TC_KPN18_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_KPN18_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_KPN18/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# KPN18: of 5630 DE features, 5326 have a RefSeq ID, of which 4919 were mapped to UniProt IDs. Finally, 3875 ( 72 % ) RefSeq IDs ( 68 % DE features ) can be annotated with GO terms

echo "MG1655: of" $(expr $(cat RNAseq_MG1655/DE_results_raw_* | wc -l) - 1) "DE features," $(cat GSEA/refseq_IDs/MG1655_*.txt | wc -l) "have a RefSeq ID, of which" $(expr $(cat GSEA/refseq2uniprot/refseq2uniprot_MG1655_*_filt.tsv | wc -l) - 1) "were mapped to UniProt IDs. Finally," $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_MG1655_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) "(" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_MG1655_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(cat GSEA/refseq_IDs/MG1655_*.txt | wc -l)) "% ) RefSeq IDs (" $(expr $(expr $(cut -f6 GSEA/refseq2uniprot/refseq2uniprot_MG1655_*_filt.tsv | grep -P "[a-z]" | grep -v "Gene Ontology" | wc -l) \* 100) / $(expr $(cat RNAseq_MG1655/DE_results_raw_* | wc -l) - 1)) "% DE features ) can be annotated with GO terms"
# MG1655: of 4413 DE features, 4252 have a RefSeq ID, of which 4011 were mapped to UniProt IDs. Finally, 3627 ( 85 % ) RefSeq IDs ( 82 % DE features ) can be annotated with GO terms
```

Adding GO annotations to the significant DEGs table:

```sh
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv <(cat RNAseq_C325/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_C325/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_C325/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv <(cat RNAseq_C325/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_C325/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_C325/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv <(cat RNAseq_CF13/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_CF13/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_CF13/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv <(cat RNAseq_CF13/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_CF13/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_CF13/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv <(cat RNAseq_H53/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_H53/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_H53/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv <(cat RNAseq_H53/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_H53/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_H53/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv <(cat RNAseq_J57/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_J57/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_J57/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv <(cat RNAseq_J57/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_J57/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_J57/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv <(cat RNAseq_K147/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_K147/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_K147/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv <(cat RNAseq_K147/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_K147/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_K147/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv <(cat RNAseq_EC10/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_EC10/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_EC10/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv <(cat RNAseq_EC10/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_EC10/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_EC10/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv <(cat RNAseq_KPN04/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN04/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN04/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv <(cat RNAseq_KPN04/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN04/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN04/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv <(cat RNAseq_KPN07/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN07/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN07/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv <(cat RNAseq_KPN07/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN07/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN07/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $17 in a {$20=a[$17]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv <(cat RNAseq_KPN10/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN10/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN10/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $17 in a {$20=a[$17]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv <(cat RNAseq_KPN10/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN10/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN10/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv <(cat RNAseq_KPN15/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN15/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN15/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv <(cat RNAseq_KPN15/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN15/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN15/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv <(cat RNAseq_KPN18/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN18/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN18/DE_results_filtered_chromosome_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $19 in a {$22=a[$19]} 1' GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv <(cat RNAseq_KPN18/DE_results_filtered_plasmids_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_KPN18/DE_results_filtered_plasmids_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_KPN18/DE_results_filtered_plasmids_padj_GOannot.tsv
awk -F'\t' 'BEGIN {OFS=FS} NR==FNR {a[$1]=$6; next} $18 in a {$21=a[$18]} 1' GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv <(cat RNAseq_MG1655/DE_results_filtered_chromosome_padj.tsv | sed 's/"//g' | sed 's/\.1\t/\t/g') > RNAseq_MG1655/DE_results_filtered_chromosome_padj_GOannot.tsv
sed -i 's/Product/Product\tGO/g' RNAseq_MG1655/DE_results_filtered_chromosome_padj_GOannot.tsv
```

Run the **DEGs_table_parser.py** script to generate a summary table of DEGs. Only includes CDSs, summarized by "Type", "RefSeq", "Gene" and "Product". Other features like tRNAs and CDSs with no RefSeq and annotated as hypothetical proteins are excluded because we cannot differenciate them.

```sh
./DEGs_table_parser.py -d RNAseq_C325/ RNAseq_EC10/ RNAseq_MG1655/ RNAseq_CF13/ RNAseq_J57/ RNAseq_H53/ RNAseq_K147/ RNAseq_KPN04/ RNAseq_KPN07/ RNAseq_KPN10/ RNAseq_KPN15/ RNAseq_KPN18/ -n DEGs_summary
```


## Gene Set Enrichment Analysis (GSEA)

GSEA will be performed with clusterProfiler, which requieres the geneList, TERM2NAME and TERM2GENE mapping files. For **geneList**, we generate a table with gene IDs (CDSs) as first column and log<sub>2</sub>FC (shrunk with DESeq2’s *apeglm* method) as second column. Note that CDSs annotated as pseudogenes are shown in two rows with identical GeneIDs and values, so we will only store the first occurrence:

```sh
mkdir GSEA/geneLists/

grep -P "\t\"CDS\"\t" RNAseq_C325/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/C325_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_C325/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/C325_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_CF13/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/CF13_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_CF13/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/CF13_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_H53/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/H53_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_H53/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/H53_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_J57/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/J57_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_J57/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/J57_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_K147/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/K147_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_K147/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/K147_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_EC10/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/EC10_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_EC10/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/EC10_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN04/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN04_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN04/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN04_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN07/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN07_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN07/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN07_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN10/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN10_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN10/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN10_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN15/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN15_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN15/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN15_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN18/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN18_chromosome.tsv
grep -P "\t\"CDS\"\t" RNAseq_KPN18/DE_results_raw_plasmids.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/KPN18_plasmids.tsv
grep -P "\t\"CDS\"\t" RNAseq_MG1655/DE_results_raw_chromosome.tsv | cut -f1,3 | sed 's/"//g' | awk '!seen[$0]++' > GSEA/geneLists/MG1655_chromosome.tsv
```

In the **TERM2NAME** table, the first column contains GO IDs and the second column the GO annotations:

```sh
mkdir GSEA/TERM2NAME/

# Biological Process
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_chromosome_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_plasmids_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/MG1655_chromosome_BP.tsv

# Cellular Component
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_chromosome_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_plasmids_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/MG1655_chromosome_CC.tsv

# Molecular Function
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/C325_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/CF13_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/H53_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/J57_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/K147_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/EC10_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN04_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN07_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN10_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN15_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_chromosome_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/KPN18_plasmids_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep "GO" | sed 's/; /\n/g' | sed 's/ \[GO/\tGO/g' | awk -F"\t" '{print $2"\t"$1}' | sed 's/\]//g' | sort | uniq > GSEA/TERM2NAME/MG1655_chromosome_MF.tsv
```

In the **TERM2GENE** table, the first column contains the GO IDs and the second column the gene IDs. First, get Gene IDs - RefSeq mappings. As before, only one row of the pseudogene CDSs is saved:

```sh
mkdir GSEA/TERM2GENE/

awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_C325/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_C325_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_C325/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_C325_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_CF13/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_CF13_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_CF13/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_CF13_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_H53/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_H53_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_H53/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_H53_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_J57/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_J57_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_J57/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_J57_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_K147/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_K147_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_K147/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_K147_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_EC10/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_EC10_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_EC10/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_EC10_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN04/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN04_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN04/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN04_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN07/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN07_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN07/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN07_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN10/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN10_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN10/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN10_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN15/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN15_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN15/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN15_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN18/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN18_chr.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_KPN18/DE_results_raw_plasmids.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_KPN18_plas.tsv
awk -F"\t" '{print $1"\t"$(NF-2)}' RNAseq_MG1655/DE_results_raw_chromosome.tsv | grep -v "\"-\"" | sed 's/"//g' | sed 's/\.[0-9]//g' | awk '!seen[$0]++' > GSEA/TERM2GENE/geneID_refseq_MG1655_chr.tsv
```

Now, get a list of GO terms (unique) by sub-ontology:

```sh
# Biological Process
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_chr_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_plas_BP.tsv
cut -f3 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_MG1655_chr_BP.tsv

# Cellular Component
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_chr_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_plas_CC.tsv
cut -f4 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_MG1655_chr_CC.tsv

# Molecular Function
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_C325_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_CF13_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_H53_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_J57_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_K147_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_EC10_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN04_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN07_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN10_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN15_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_chr_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_KPN18_plas_MF.tsv
cut -f5 GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv | grep -oP "GO\:[0-9]{7}" | sort | uniq > GSEA/TERM2GENE/GOs_MG1655_chr_MF.tsv
```

Now, run the provided **generate_TERM2GENE.py** script to obtain the TERM2GENE table, providing the list of GO terms, the annotation file of RefSeq IDs with GOs and the file mapping GeneIDs to RefSeq IDs:

```sh
# Biological Process
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_chr_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_plas_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_plas_BP.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_MG1655_chr_BP.tsv -a GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_MG1655_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_MG1655_chr_BP.tsv

# Cellular Component
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_chr_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_plas_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_plas_CC.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_MG1655_chr_CC.tsv -a GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_MG1655_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_MG1655_chr_CC.tsv

# Molecular Function
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_C325_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_C325_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_C325_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_C325_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_CF13_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_CF13_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_CF13_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_CF13_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_H53_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_H53_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_H53_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_H53_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_J57_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_J57_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_J57_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_J57_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_K147_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_K147_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_K147_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_K147_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_EC10_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_EC10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_EC10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_EC10_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN04_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN04_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN04_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN04_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN07_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN07_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN07_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN07_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN10_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN10_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN10_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN10_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN15_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN15_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN15_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN15_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_chr_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_KPN18_plas_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_KPN18_plas_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_KPN18_plas.tsv -o GSEA/TERM2GENE/TERM2GENE_KPN18_plas_MF.tsv
./GSEA/generate_TERM2GENE.py -g GSEA/TERM2GENE/GOs_MG1655_chr_MF.tsv -a GSEA/refseq2uniprot/refseq2uniprot_MG1655_chr_filt.tsv -m GSEA/TERM2GENE/geneID_refseq_MG1655_chr.tsv -o GSEA/TERM2GENE/TERM2GENE_MG1655_chr_MF.tsv
```

Everything is set to run GSEA with **clusterProfiler v4.8.1**, correcting for multiple tests with the Benjamini-Hochberg procedure, running the R script **GSEA.Rmd** (check the HTML-formatted output [here](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/GSEA/GSEA.html)). For representing Fig 2 and Fig S6, GO terms were grouped by parental terms based on observing the GO tree in QuickGO (June 2023).


## Transcript Per Million (TPM)

In each strain's *RNAseq_* directory, the Rmd script **TPM_calculation.Rmd** calculates the TPM values of pOXA-48-carrying strains from FPKM values computed with DESeq2 from raw counts. Then, the script **TPM_calculation.Rmd** inside the directory `plot_TPM/` imports the previously generated TPM values, normalizes them by median chromosomal TPM and plots a heatmap of pOXA-48 gene expression (Fig 1). Check the HTML-formatted output [here](https://laboratoribio.github.io/RNA-Seq_enterobacteria_pOXA-48/RNA-Seq_1st_dataset/plot_TPM/TPM_calculation.html).


## Variant calling

To assess whether mutations accumulated in the chromosome or other plasmids during growth or construction of pOXA-48-carrying or pOXA-48-free samples affected gene expression, variant calling was performed on the transcribed regions with **Snippy v4.6.0**.

```sh
# C325
snippy --outdir variants_snippy/C325.1 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325.1.bam
snippy --outdir variants_snippy/C325.2 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325.2.bam
snippy --outdir variants_snippy/C325.3 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325.3.bam
snippy --outdir variants_snippy/C325c1.1 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325c1.1.bam
snippy --outdir variants_snippy/C325c1.2 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325c1.2.bam
snippy --outdir variants_snippy/C325c1.3 --ref ../../Closed_sequences/C325.gbk --bam RNAseq_C325/C325c1.3.bam

#CF13
snippy --outdir variants_snippy/CF13.1 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13.1.bam
snippy --outdir variants_snippy/CF13.2 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13.2.bam
snippy --outdir variants_snippy/CF13.3 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13.3.bam
snippy --outdir variants_snippy/CF13c1.1 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13c1.1.bam
snippy --outdir variants_snippy/CF13c1.2 --ref ../../Closed_sequences/CF13.gbk --bam RNAseq_CF13/CF13c1.2.bam

# H53
snippy --outdir variants_snippy/H53.1 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53.1.bam
snippy --outdir variants_snippy/H53.2 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53.2.bam
snippy --outdir variants_snippy/H53.3 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53.3.bam
snippy --outdir variants_snippy/H53c1.1 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53c1.1.bam
snippy --outdir variants_snippy/H53c1.2 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53c1.2.bam
snippy --outdir variants_snippy/H53c1.3 --ref ../../Closed_sequences/H53.gbk --bam RNAseq_H53/H53c1.3.bam

# J57
snippy --outdir variants_snippy/J57.1 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57.1.bam
snippy --outdir variants_snippy/J57.2 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57.2.bam
snippy --outdir variants_snippy/J57.3 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57.3.bam
snippy --outdir variants_snippy/J57c1.1 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57c1.1.bam
snippy --outdir variants_snippy/J57c1.2 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57c1.2.bam
snippy --outdir variants_snippy/J57c1.3 --ref ../../Closed_sequences/J57.gbk --bam RNAseq_J57/J57c1.3.bam

# K147
snippy --outdir variants_snippy/K147.1 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147.1.bam
snippy --outdir variants_snippy/K147.2 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147.2.bam
snippy --outdir variants_snippy/K147.3 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147.3.bam
snippy --outdir variants_snippy/K147.4 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147.4.bam
snippy --outdir variants_snippy/K147c1.1 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147c1.1.bam
snippy --outdir variants_snippy/K147c1.2 --ref ../../Closed_sequences/K147.gbk --bam RNAseq_K147/K147c1.2.bam

# EC10
snippy --outdir variants_snippy/PF_EC10.1 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/PF_EC10.1.bam
snippy --outdir variants_snippy/PF_EC10.2 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/PF_EC10.2.bam
snippy --outdir variants_snippy/PF_EC10.3 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/PF_EC10.3.bam
snippy --outdir variants_snippy/TC_EC10.1 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/TC_EC10.1.bam
snippy --outdir variants_snippy/TC_EC10.2 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/TC_EC10.2.bam
snippy --outdir variants_snippy/TC_EC10.3 --ref ../../Closed_sequences/TC_EC10.gbk --bam RNAseq_EC10/TC_EC10.3.bam

# KPN04
snippy --outdir variants_snippy/PF_KPN04.1 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/PF_KPN04.1.bam
snippy --outdir variants_snippy/PF_KPN04.2 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/PF_KPN04.2.bam
snippy --outdir variants_snippy/PF_KPN04.3 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/PF_KPN04.3.bam
snippy --outdir variants_snippy/TC_KPN04.1 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/TC_KPN04.1.bam
snippy --outdir variants_snippy/TC_KPN04.2 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/TC_KPN04.2.bam
snippy --outdir variants_snippy/TC_KPN04.3 --ref ../../Closed_sequences/TC_KPN04.gbk --bam RNAseq_KPN04/TC_KPN04.3.bam

# KPN07
snippy --outdir variants_snippy/PF_KPN07.1 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/PF_KPN07.1.bam
snippy --outdir variants_snippy/PF_KPN07.2 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/PF_KPN07.2.bam
snippy --outdir variants_snippy/PF_KPN07.3 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/PF_KPN07.3.bam
snippy --outdir variants_snippy/TC_KPN07.1 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/TC_KPN07.1.bam
snippy --outdir variants_snippy/TC_KPN07.2 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/TC_KPN07.2.bam
snippy --outdir variants_snippy/TC_KPN07.3 --ref ../../Closed_sequences/TC_KPN07.gbk --bam RNAseq_KPN07/TC_KPN07.3.bam

# KPN10
snippy --outdir variants_snippy/PF_KPN10.1 --ref ../../Closed_sequences/TC_KPN10.gbk --bam RNAseq_KPN10/PF_KPN10.1.bam
snippy --outdir variants_snippy/PF_KPN10.2 --ref ../../Closed_sequences/TC_KPN10.gbk --bam RNAseq_KPN10/PF_KPN10.2.bam
snippy --outdir variants_snippy/TC_KPN10.1 --ref ../../Closed_sequences/TC_KPN10.gbk --bam RNAseq_KPN10/TC_KPN10.1.bam
snippy --outdir variants_snippy/TC_KPN10.2 --ref ../../Closed_sequences/TC_KPN10.gbk --bam RNAseq_KPN10/TC_KPN10.2.bam

# KPN15
snippy --outdir variants_snippy/PF_KPN15.1 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/PF_KPN15.1.bam
snippy --outdir variants_snippy/PF_KPN15.2 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/PF_KPN15.2.bam
snippy --outdir variants_snippy/TC_KPN15.1 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/TC_KPN15.1.bam
snippy --outdir variants_snippy/TC_KPN15.2 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/TC_KPN15.2.bam
snippy --outdir variants_snippy/TC_KPN15.3 --ref ../../Closed_sequences/TC_KPN15.gbk --bam RNAseq_KPN15/TC_KPN15.3.bam

# KPN18
snippy --outdir variants_snippy/PF_KPN18.1 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/PF_KPN18.1.bam
snippy --outdir variants_snippy/PF_KPN18.2 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/PF_KPN18.2.bam
snippy --outdir variants_snippy/PF_KPN18.3 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/PF_KPN18.3.bam
snippy --outdir variants_snippy/TC_KPN18.1 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/TC_KPN18.1.bam
snippy --outdir variants_snippy/TC_KPN18.2 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/TC_KPN18.2.bam
snippy --outdir variants_snippy/TC_KPN18.3 --ref ../../Closed_sequences/TC_KPN18.gbk --bam RNAseq_KPN18/TC_KPN18.3.bam

# MG1655
snippy --outdir variants_snippy/MG1655.1 --ref ../../Closed_sequences/MG1655p.gbk --bam RNAseq_MG1655/MG1655.1.bam
snippy --outdir variants_snippy/MG1655.2 --ref ../../Closed_sequences/MG1655p.gbk --bam RNAseq_MG1655/MG1655.2.bam
snippy --outdir variants_snippy/MG1655p.1 --ref ../../Closed_sequences/MG1655p.gbk --bam RNAseq_MG1655/MG1655p.1.bam
snippy --outdir variants_snippy/MG1655p.2 --ref ../../Closed_sequences/MG1655p.gbk --bam RNAseq_MG1655/MG1655p.2.bam
snippy --outdir variants_snippy/MG1655p.3 --ref ../../Closed_sequences/MG1655p.gbk --bam RNAseq_MG1655/MG1655p.3.bam
```

To confirm or discard suspected SNPs in plasmid pOXA-48, the RNA-Seq reads were mapped to the sequence of the pOXA-48_K8 variant (accession number MT441554):

```sh
# C325
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_C325.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_C325.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_C325.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_C325.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_C325.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_C325.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_C325.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_C325.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_C325.3_val_2.fq.gz

# CF13
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_CF13.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_CF13.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_CF13.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_CF13.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_CF13.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_CF13.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_CF13.3_val_2.fq.gz

# H53
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_H53.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_H53.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_H53.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_H53.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_H53.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_H53.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_H53.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_H53.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_H53.3_val_2.fq.gz

# J57
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_J57.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_J57.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_J57.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_J57.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_J57.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_J57.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_J57.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_J57.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_J57.3_val_2.fq.gz

# K147
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_K147.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_K147.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_K147.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_K147.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_K147.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_K147.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_K147.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_K147.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_K147.3_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_K147.4 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_K147.4_val_1.fq.gz --R2 reads_RNAseq/WTCHG_K147.4_val_2.fq.gz

# EC10
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_EC10.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_EC10.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_EC10.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_EC10.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_EC10.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_EC10.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_EC10.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_EC10.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_EC10.3_val_2.fq.gz

# KPN04
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN04.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN04.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN04.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN04.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN04.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN04.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN04.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN04.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN04.3_val_2.fq.gz

# KPN07
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN07.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN07.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN07.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN07.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN07.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN07.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN07.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN07.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN07.3_val_2.fq.gz

# KPN10
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN10.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN10.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN10.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN10.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN10.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN10.2_val_2.fq.gz

# KPN15
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN15.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN15.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN15.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN15.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN15.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN15.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN15.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN15.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN15.3_val_2.fq.gz

# KPN18
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN18.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN18.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN18.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN18.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN18.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN18.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_TC_KPN18.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_TC_KPN18.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_TC_KPN18.3_val_2.fq.gz

# MG1655
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_MG1655p.1 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_MG1655p.1_val_1.fq.gz --R2 reads_RNAseq/WTCHG_MG1655p.1_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_MG1655p.2 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_MG1655p.2_val_1.fq.gz --R2 reads_RNAseq/WTCHG_MG1655p.2_val_2.fq.gz
snippy --outdir variants_snippy/ref_pOXA-48_K8_map_MG1655p.3 --ref ../../Closed_sequences/plasmids/pOXA-48_K8.gb --R1 reads_RNAseq/WTCHG_MG1655p.3_val_1.fq.gz --R2 reads_RNAseq/WTCHG_MG1655p.3_val_2.fq.gz
```


## Checking LysR<sub>pOXA-48</sub> identity with other LysRs

The protein sequence of LysR<sub>pOXA-48</sub> (place **lysR_pOXA-48.faa** in `./blastp_LysR_pOXA-48_to_refs/`) was blasted (**BLASTp v2.11.0**) against the reference genomes to assess the percentage of identity of LysR<sub>pOXA-48</sub> to other chromosomally-encoded LysRs:

```sh
mkdir ./blastp_LysR_pOXA-48_to_refs

makeblastdb -in ../../Closed_sequences/C325.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/C325.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_C325.tsv
makeblastdb -in ../../Closed_sequences/TC_EC10.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_EC10.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_EC10.tsv
makeblastdb -in ../../Closed_sequences/MG1655p.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/MG1655p.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_MG1655p.tsv
makeblastdb -in ../../Closed_sequences/CF13.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/CF13.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_CF13.tsv
makeblastdb -in ../../Closed_sequences/H53.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/H53.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_H53.tsv
makeblastdb -in ../../Closed_sequences/J57.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/J57.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_J57.tsv
makeblastdb -in ../../Closed_sequences/K147.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/K147.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_K147.tsv
makeblastdb -in ../../Closed_sequences/TC_KPN04.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_KPN04.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_KPN04.tsv
makeblastdb -in ../../Closed_sequences/TC_KPN07.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_KPN07.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_KPN07.tsv
makeblastdb -in ../../Closed_sequences/TC_KPN10.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_KPN10.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_KPN10.tsv
makeblastdb -in ../../Closed_sequences/TC_KPN15.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_KPN15.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_KPN15.tsv
makeblastdb -in ../../Closed_sequences/TC_KPN18.faa -dbtype prot
blastp -query blastp_LysR_pOXA-48_to_refs/lysR_pOXA-48.faa -db ../../Closed_sequences/TC_KPN18.faa -outfmt 6 > blastp_LysR_pOXA-48_to_refs/LysR_pOXA-48_to_TC_KPN18.tsv
```
