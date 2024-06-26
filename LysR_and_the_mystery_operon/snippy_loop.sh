#!/bin/bash

for contigs in /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/dnds_analysis/input_refseq_sequences/*.fasta
do
	sample=$( echo ${contigs%%.fasta})
	sample=$( basename $sample )
	echo snippy --outdir /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/dnds_analysis/variant_calling_snippy/$sample --prefix $sample --ref reference/pOXA48K8.gbk --ctgs $contigs --cpus 17 --force
done


