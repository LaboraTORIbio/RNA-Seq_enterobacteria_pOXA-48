#!/bin/bash

# Extract genomic regions of operon and surrounding areas (+-10kbp) from the coords given by macsyfinder

samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/AGRTUM.0723.00003/GCF_002591665.3_ASM259166v3_genomic.fna NZ_CP042274.1:2077513-2099945 > ./clinker_analysis/input_fasta/AGRTUM.0723.00003_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/AGRTUM.0723.00040/GCF_025560485.1_ASM2556048v1_genomic.fna NZ_CP048551.1:2272898-2295607 > ./clinker_analysis/input_fasta/AGRTUM.0723.00040_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ALCFAE.0723.00002/GCF_001641975.2_ASM164197v2_genomic.fna NZ_CP021079.1:3771195-3793867 > ./clinker_analysis/input_fasta/ALCFAE.0723.00002_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ALCFAE.0723.00005/GCF_002443155.1_ASM244315v1_genomic.fna NZ_CP023667.1:3515389-3538034 > ./clinker_analysis/input_fasta/ALCFAE.0723.00005_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/BORTRE.0723.00001/GCF_002860045.1_ASM286004v1_genomic.fna NZ_CP018898.1:1945858-1968470 > ./clinker_analysis/input_fasta/BORTRE.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/BRUSP..0723.00013/GCF_025908475.1_ASM2590847v1_genomic.fna NZ_CP109781.1:1443853-1466559 > ./clinker_analysis/input_fasta/BRUSP..0723.00013_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/BURCEP.0723.00020/GCF_029962485.1_ASM2996248v1_genomic.fna NZ_CP073638.1:321358-344188 > ./clinker_analysis/input_fasta/BURCEP.0723.00020_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/CF13/CF13.fasta 1:1670069-1692643 > ./clinker_analysis/input_fasta/CF13_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/CITFRE.0723.00011/GCF_002903215.1_ASM290321v1_genomic.fna NZ_CP026235.1:1682239-1704813 > ./clinker_analysis/input_fasta/CITFRE.0723.00011_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/CITSP..0723.00032/GCF_013889015.1_ASM1388901v1_genomic.fna NZ_CP057611.1:2831853-2854425 > ./clinker_analysis/input_fasta/CITSP..0723.00032_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ENTCLO.0723.00041/GCF_014169315.1_ASM1416931v1_genomic.fna NZ_AP022133.1:1657955-1680528 > ./clinker_analysis/input_fasta/ENTCLO.0723.00041_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ENTSP..0723.00014/GCF_010692925.1_ASM1069292v1_genomic.fna NZ_CP048736.1:1700728-1723300 > ./clinker_analysis/input_fasta/ENTSP..0723.00014_op10k.fasta
### Second Enterobacter spp contains operon in an IncFII plasmid, needed to circularize and extract regions from end and start of the plasmid
# Start of plasmid
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ENTSP..0723.00039/GCF_019968805.1_ASM1996880v1_genomic.fna NZ_CP074160.1:1-15487 > ./clinker_analysis/input_fasta/ENTSP..0723.00039_op10k_1.fasta
# End of plasmid
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ENTSP..0723.00039/GCF_019968805.1_ASM1996880v1_genomic.fna NZ_CP074160.1:184116-191037 > ./clinker_analysis/input_fasta/ENTSP..0723.00039_op10k_2.fasta
# Then merged using https://www.bioinformatics.org/sms2/combine_fasta.html

samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ESCCOL.0723.02099/GCF_022493275.1_ASM2249327v1_genomic.fna NZ_CP028742.1:4894420-4916938 > ./clinker_analysis/input_fasta/ESCCOL.0723.02099_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/ESCCOL.0723.02103/GCF_022493355.1_ASM2249335v1_genomic.fna NZ_CP028747.1:859126-881650 > ./clinker_analysis/input_fasta/ESCCOL.0723.02103_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/J57/J57.fasta 1:3524120-3546696 > ./clinker_analysis/input_fasta/J57_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/KLEAER.0723.00045/GCF_029027985.1_ASM2902798v1_genomic.fna NZ_CP119076.1:2998529-3021107 > ./clinker_analysis/input_fasta/KLEAER.0723.00045_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/KLEOXY.0723.00013/GCF_009707385.1_ASM970738v1_genomic.fna NZ_CP046115.1:2147543-2170121 > ./clinker_analysis/input_fasta/KLEOXY.0723.00013_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/KLEPNE.0723.00003/GCF_000220485.1_ASM22048v1_genomic.fna NC_017540.1:1753066-1775642 > ./clinker_analysis/input_fasta/KLEPNE.0723.00003_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/KLUASC.0723.00001/GCF_023195735.1_ASM2319573v1_genomic.fna NZ_CP096201.1:3111979-3134556 > ./clinker_analysis/input_fasta/KLUASC.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/PSEAER.0723.00007/GCF_000223945.1_ASM22394v2_genomic.fna NZ_AFXJ01000001.1:2740078-2762756 > ./clinker_analysis/input_fasta/PSEAER.0723.00007_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/RAOORN.0723.00004/GCF_001723565.1_ASM172356v1_genomic.fna NZ_CP012555.1:1773023-1795493 > ./clinker_analysis/input_fasta/RAOORN.0723.00004_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/RAOPLA.0723.00001/GCF_000783935.2_ASM78393v2_genomic.fna NZ_CP026047.1:3740387-3762964 > ./clinker_analysis/input_fasta/RAOPLA.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SALENT.0723.00028/GCF_000240905.2_ASM24090v3_genomic.fna NZ_CP019186.1:2758837-2781429 > ./clinker_analysis/input_fasta/SALENT.0723.00028_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SALENT.0723.00459/GCF_002211965.1_ASM221196v1_genomic.fna NZ_CP022142.1:3384240-3406817 > ./clinker_analysis/input_fasta/SALENT.0723.00459_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SERMAR.0723.00023/GCF_003031545.1_ASM303154v1_genomic.fna NZ_CP020503.1:983982-1006642 > ./clinker_analysis/input_fasta/SERMAR.0723.00023_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SHEDEN.0723.00001/GCF_000013765.1_ASM1376v1_genomic.fna NC_007954.1:4501528-4524108 > ./clinker_analysis/input_fasta/SHEDEN.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SHEPSY.0723.00001/GCF_002005305.1_ASM200530v1_genomic.fna NZ_CP014782.1:5634217-5656612 > ./clinker_analysis/input_fasta/SHEPSY.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SHIBLA.0723.00001/GCF_000262305.1_ASM26230v1_genomic.fna NC_017910.1:2212066-2234577 > ./clinker_analysis/input_fasta/SHIBLA.0723.00001_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/SHIBLA.0723.00002/GCF_900635135.1_36603_E01_genomic.fna NZ_LR133996.1:2210456-2232733 > ./clinker_analysis/input_fasta/SHIBLA.0723.00002_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/TC_KPN15/TC_KPN15.fasta 1:3237851-3260427 > ./clinker_analysis/input_fasta/TC_KPN15_op10k.fasta
samtools faidx /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/VIBCHO.0723.00072/GCF_013357745.1_ASM1335774v1_genomic.fna NZ_CP053807.1:267378-289889 > ./clinker_analysis/input_fasta/VIBCHO.0723.00072_op10k.fasta

# Detect IS elements in surrounding regions

for genome in /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/clinker_analysis/input_fasta/*.fasta
do
	strain=$( basename $genome )
	strain=${strain::-6}
	#echo $genome
	echo isescan.py --seqfile $genome --nthread 16 --output /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/isescan_synteny/$strain
done

# Annotate with prokka as gbk for clinker input

for genome in /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/clinker_analysis/input_fasta/*_op10k.fasta
do
	sample=$( basename $genome )
	samplename=${sample%%_op10k.fasta}
	#echo $samplename
	echo prokka --outdir /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/clinker_analysis/prokka_annots/$samplename --proteins /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/synteny/$samplename/$samplename.faa --cpus 16 $genome
done

# Execute clinker

clinker /media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/clinker_analysis/prokka_annots/*.gbk -p

