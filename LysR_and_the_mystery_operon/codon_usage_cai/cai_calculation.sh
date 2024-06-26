#########################################################################
#### Compute the CAI of chromosomal genes and the operon system genes ###
#########################################################################

# The Codon Adaptation Index (CAI) was used to assess the similarity between the synonymous codon usage of the operon system genes and that of the chromosomal ribosomal proteins as reference, since they are highly expressed, and highly expressed genes have an optimized codon usage (Sharp and Li 1987).

# Run from directory ./codon_usage_cai


# Get seqs of lysR, pirin and isochorismatase (in that order, in forward strand)
samtools faidx ../../../Closed_sequences/CF13.fasta 1:1681732-1682643 > operon_seqs_CF13.fasta
samtools faidx ../../../Closed_sequences/CF13.fasta 1:1680763-1681623 -i >> operon_seqs_CF13.fasta
samtools faidx ../../../Closed_sequences/CF13.fasta 1:1680069-1680746 -i >> operon_seqs_CF13.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta 1:3247851-3248756 -i > operon_seqs_TC_KPN15.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta 1:3248871-3249731 >> operon_seqs_KPN15.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta 1:3249750-3250427 >> operon_seqs_KPN15.fasta
samtools faidx ../../../Closed_sequences/J57.fasta 1:3534120-3535025 -i > operon_seqs_J57.fasta
samtools faidx ../../../Closed_sequences/J57.fasta 1:3535140-3536000 >> operon_seqs_J57.fasta
samtools faidx ../../../Closed_sequences/J57.fasta 1:3536019-3536696 >> operon_seqs_J57.fasta
samtools faidx ../GC_content/SHEDEN.0723.00001.fna NC_007954.1:4513203-4514108 > operon_seqs_SHEDEN.0723.00001.fasta
samtools faidx ../GC_content/SHEDEN.0723.00001.fna NC_007954.1:4512218-4513081 -i >> operon_seqs_SHEDEN.0723.00001.fasta
samtools faidx ../GC_content/SHEDEN.0723.00001.fna NC_007954.1:4511528-4512154 -i >> operon_seqs_SHEDEN.0723.00001.fasta
samtools faidx ../GC_content/SERMAR.0723.00023.fna NZ_CP020503.1:993982-994890 -i > operon_seqs_SERMAR.0723.00023.fasta
samtools faidx ../GC_content/SERMAR.0723.00023.fna NZ_CP020503.1:995016-995897 >> operon_seqs_SERMAR.0723.00023.fasta
samtools faidx ../GC_content/SERMAR.0723.00023.fna NZ_CP020503.1:996016-996642 >> operon_seqs_SERMAR.0723.00023.fasta
samtools faidx ../GC_content/ESCCOL.0723.02103.fna NZ_CP028747.1:870721-871650 > operon_seqs_ESCCOL.0723.02103.fasta
samtools faidx ../GC_content/ESCCOL.0723.02103.fna NZ_CP028747.1:869126-869818 -i >> operon_seqs_ESCCOL.0723.02103.fasta
samtools faidx ../GC_content/ESCCOL.0723.02103.fna NZ_CP028747.1:869859-870584 -i >> operon_seqs_ESCCOL.0723.02103.fasta
samtools faidx ../GC_content/SALENT.0723.00028.fna NZ_CP019186.1:2768837-2769763 -i > operon_seqs_SALENT.0723.00028.fasta
samtools faidx ../GC_content/SALENT.0723.00028.fna NZ_CP019186.1:2769872-2770732 >> operon_seqs_SALENT.0723.00028.fasta
samtools faidx ../GC_content/SALENT.0723.00028.fna NZ_CP019186.1:2770752-2771429 >> operon_seqs_SALENT.0723.00028.fasta
samtools faidx ../GC_content/ENTCLO.0723.00001.fna NC_014121.1:2647683-2648588 -i > operon_seqs_ENTCLO.0723.00001.fasta
samtools faidx ../GC_content/ENTCLO.0723.00001.fna NC_014121.1:2648698-2649558 >> operon_seqs_ENTCLO.0723.00001.fasta
samtools faidx ../GC_content/ENTCLO.0723.00001.fna NC_014121.1:2649579-2650256 >> operon_seqs_ENTCLO.0723.00001.fasta


# Get regions of ribosomal proteins
grep "ribosomal protein" ../RNAseq_CF13/*saf | awk '{if ($5 == "+") { print $2":"$3"-"$4} }' > ribosomal_proteins_plus_CF13.txt
grep "ribosomal protein" ../RNAseq_CF13/*saf | awk '{if ($5 == "-") { print $2":"$3"-"$4} }' > ribosomal_proteins_minus_CF13.txt
grep "ribosomal protein" ../RNAseq_KPN15/*saf | awk '{if ($5 == "+") { print $2":"$3"-"$4} }' > ribosomal_proteins_plus_KPN15.txt
grep "ribosomal protein" ../RNAseq_KPN15/*saf | awk '{if ($5 == "-") { print $2":"$3"-"$4} }' > ribosomal_proteins_minus_KPN15.txt
grep "ribosomal protein" ../RNAseq_J57/*saf | awk '{if ($5 == "+") { print $2":"$3"-"$4} }' > ribosomal_proteins_plus_J57.txt
grep "ribosomal protein" ../RNAseq_J57/*saf | awk '{if ($5 == "-") { print $2":"$3"-"$4} }' > ribosomal_proteins_minus_J57.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SHEDEN.0723.00001/*gff | grep "ribosomal protein" | awk '{if ($8 == "+") { print $1":"$5"-"$6} }' > ribosomal_proteins_plus_SHEDEN.0723.00001.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SHEDEN.0723.00001/*gff | grep "ribosomal protein" | awk '{if ($8 == "-") { print $1":"$5"-"$6} }' > ribosomal_proteins_minus_SHEDEN.0723.00001.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SERMAR.0723.00023/*gff | grep "ribosomal protein" | awk '{if ($8 == "+") { print $1":"$5"-"$6} }' > ribosomal_proteins_plus_SERMAR.0723.00023.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SERMAR.0723.00023/*gff | grep "ribosomal protein" | awk '{if ($8 == "-") { print $1":"$5"-"$6} }' > ribosomal_proteins_minus_SERMAR.0723.00023.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/ESCCOL.0723.02103/*gff | grep "ribosomal protein" | awk '{if ($8 == "+") { print $1":"$5"-"$6} }' > ribosomal_proteins_plus_ESCCOL.0723.02103.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/ESCCOL.0723.02103/*gff | grep "ribosomal protein" | awk '{if ($8 == "-") { print $1":"$5"-"$6} }' > ribosomal_proteins_minus_ESCCOL.0723.02103.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SALENT.0723.00028/*gff | grep "ribosomal protein" | awk '{if ($8 == "+") { print $1":"$5"-"$6} }' > ribosomal_proteins_plus_SALENT.0723.00028.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/SALENT.0723.00028/*gff | grep "ribosomal protein" | awk '{if ($8 == "-") { print $1":"$5"-"$6} }' > ribosomal_proteins_minus_SALENT.0723.00028.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/ENTCLO.0723.00001/*gff | grep "ribosomal protein" | awk '{if ($8 == "+") { print $1":"$5"-"$6} }' > ribosomal_proteins_plus_ENTCLO.0723.00001.txt
grep -P "\tCDS\t" ../../LysR_and_the_mystery_operon/ncbi_dataset_*/ENTCLO.0723.00001/*gff | grep "ribosomal protein" | awk '{if ($8 == "-") { print $1":"$5"-"$6} }' > ribosomal_proteins_minus_ENTCLO.0723.00001.txt

# Get seqs of ribosomal proteins
samtools faidx ../../../Closed_sequences/CF13.fasta -r ribosomal_proteins_plus_CF13.txt > ribosomal_proteins_seqs_CF13.fasta
samtools faidx ../../../Closed_sequences/CF13.fasta -r ribosomal_proteins_minus_CF13.txt -i >> ribosomal_proteins_seqs_CF13.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta -r ribosomal_proteins_plus_KPN15.txt > ribosomal_proteins_seqs_KPN15.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta -r ribosomal_proteins_minus_KPN15.txt -i >> ribosomal_proteins_seqs_KPN15.fasta
samtools faidx ../../../Closed_sequences/J57.fasta -r ribosomal_proteins_plus_J57.txt > ribosomal_proteins_seqs_J57.fasta
samtools faidx ../../../Closed_sequences/J57.fasta -r ribosomal_proteins_minus_J57.txt -i >> ribosomal_proteins_seqs_J57.fasta
samtools faidx ../GC_content/SHEDEN.0723.00001.fna -r ribosomal_proteins_plus_SHEDEN.0723.00001.txt > ribosomal_proteins_seqs_SHEDEN.0723.00001.fasta
samtools faidx ../GC_content/SHEDEN.0723.00001.fna -r ribosomal_proteins_minus_SHEDEN.0723.00001.txt -i >> ribosomal_proteins_seqs_SHEDEN.0723.00001.fasta
samtools faidx ../GC_content/SERMAR.0723.00023.fna -r ribosomal_proteins_plus_SERMAR.0723.00023.txt > ribosomal_proteins_seqs_SERMAR.0723.00023.fasta
samtools faidx ../GC_content/SERMAR.0723.00023.fna -r ribosomal_proteins_minus_SERMAR.0723.00023.txt -i >> ribosomal_proteins_seqs_SERMAR.0723.00023.fasta
samtools faidx ../GC_content/ESCCOL.0723.02103.fna -r ribosomal_proteins_plus_ESCCOL.0723.02103.txt > ribosomal_proteins_seqs_ESCCOL.0723.02103.fasta
samtools faidx ../GC_content/ESCCOL.0723.02103.fna -r ribosomal_proteins_minus_ESCCOL.0723.02103.txt -i >> ribosomal_proteins_seqs_ESCCOL.0723.02103.fasta
samtools faidx ../GC_content/SALENT.0723.00028.fna -r ribosomal_proteins_plus_SALENT.0723.00028.txt > ribosomal_proteins_seqs_SALENT.0723.00028.fasta
samtools faidx ../GC_content/SALENT.0723.00028.fna -r ribosomal_proteins_minus_SALENT.0723.00028.txt -i >> ribosomal_proteins_seqs_SALENT.0723.00028.fasta
samtools faidx ../GC_content/ENTCLO.0723.00001.fna -r ribosomal_proteins_plus_ENTCLO.0723.00001.txt > ribosomal_proteins_seqs_ENTCLO.0723.00001.fasta
samtools faidx ../GC_content/ENTCLO.0723.00001.fna -r ribosomal_proteins_minus_ENTCLO.0723.00001.txt -i >> ribosomal_proteins_seqs_ENTCLO.0723.00001.fasta



# Make codon usage table
cusp -sequence ribosomal_proteins_seqs_CF13.fasta -outfile ribosomal_proteins_cusp_CF13.cut
cusp -sequence ribosomal_proteins_seqs_KPN15.fasta -outfile ribosomal_proteins_cusp_KPN15.cut
cusp -sequence ribosomal_proteins_seqs_J57.fasta -outfile ribosomal_proteins_cusp_J57.cut
cusp -sequence ribosomal_proteins_seqs_SHEDEN.0723.00001.fasta -outfile ribosomal_proteins_cusp_SHEDEN.0723.00001.cut
cusp -sequence ribosomal_proteins_seqs_SERMAR.0723.00023.fasta -outfile ribosomal_proteins_cusp_SERMAR.0723.00023.cut
cusp -sequence ribosomal_proteins_seqs_ESCCOL.0723.02103.fasta -outfile ribosomal_proteins_cusp_ESCCOL.0723.02103.cut
cusp -sequence ribosomal_proteins_seqs_SALENT.0723.00028.fasta -outfile ribosomal_proteins_cusp_SALENT.0723.00028.cut
cusp -sequence ribosomal_proteins_seqs_ENTCLO.0723.00001.fasta -outfile ribosomal_proteins_cusp_ENTCLO.0723.00001.cut



# Compute CAI of ribosomal proteins, for reference
cai ribosomal_proteins_seqs_CF13.fasta -cfile ribosomal_proteins_cusp_CF13.cut -outfile cai_ribosomal_proteins_CF13.txt
cai ribosomal_proteins_seqs_KPN15.fasta -cfile ribosomal_proteins_cusp_KPN15.cut -outfile cai_ribosomal_proteins_KPN15.txt
cai ribosomal_proteins_seqs_J57.fasta -cfile ribosomal_proteins_cusp_J57.cut -outfile cai_ribosomal_proteins_J57.txt
cai ribosomal_proteins_seqs_SHEDEN.0723.00001.fasta -cfile ribosomal_proteins_cusp_SHEDEN.0723.00001.cut -outfile cai_ribosomal_proteins_SHEDEN.0723.00001.txt
cai ribosomal_proteins_seqs_SERMAR.0723.00023.fasta -cfile ribosomal_proteins_cusp_SERMAR.0723.00023.cut -outfile cai_ribosomal_proteins_SERMAR.0723.00023.txt
cai ribosomal_proteins_seqs_ESCCOL.0723.02103.fasta -cfile ribosomal_proteins_cusp_ESCCOL.0723.02103.cut -outfile cai_ribosomal_proteins_ESCCOL.0723.02103.txt
cai ribosomal_proteins_seqs_SALENT.0723.00028.fasta -cfile ribosomal_proteins_cusp_SALENT.0723.00028.cut -outfile cai_ribosomal_proteins_SALENT.0723.00028.txt
cai ribosomal_proteins_seqs_ENTCLO.0723.00001.fasta -cfile ribosomal_proteins_cusp_ENTCLO.0723.00001.cut -outfile cai_ribosomal_proteins_ENTCLO.0723.00001.txt

# Compute CAI of lysR, pirin and isochorismatase
cai operon_seqs_CF13.fasta -cfile ribosomal_proteins_cusp_CF13.cut -outfile cai_operon_CF13.txt
cai operon_seqs_KPN15.fasta -cfile ribosomal_proteins_cusp_KPN15.cut -outfile cai_operon_KPN15.txt
cai operon_seqs_J57.fasta -cfile ribosomal_proteins_cusp_J57.cut -outfile cai_operon_J57.txt
cai operon_seqs_SHEDEN.0723.00001.fasta -cfile ribosomal_proteins_cusp_SHEDEN.0723.00001.cut -outfile cai_operon_SHEDEN.0723.00001.txt
cai operon_seqs_SERMAR.0723.00023.fasta -cfile ribosomal_proteins_cusp_SERMAR.0723.00023.cut -outfile cai_operon_SERMAR.0723.00023.txt
cai operon_seqs_ESCCOL.0723.02103.fasta -cfile ribosomal_proteins_cusp_ESCCOL.0723.02103.cut -outfile cai_operon_ESCCOL.0723.02103.txt
cai operon_seqs_SALENT.0723.00028.fasta -cfile ribosomal_proteins_cusp_SALENT.0723.00028.cut -outfile cai_operon_SALENT.0723.00028.txt
cai operon_seqs_ENTCLO.0723.00001.fasta -cfile ribosomal_proteins_cusp_ENTCLO.0723.00001.cut -outfile cai_operon_ENTCLO.0723.00001.txt



# Generate table
echo -e "Strain\tFeature\t(median)CAI"

cut -d" " -f4 cai_ribosomal_proteins_CF13.txt | sort | awk ' { a[i++]=$1; } END { print "CF13\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "CF13\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "CF13\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "CF13\tifp\t"$4}'
    fi
done < cai_operon_CF13.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_KPN15.txt | sort | awk ' { a[i++]=$1; } END { print "KPN15\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "KPN15\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "KPN15\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "KPN15\tifp\t"$4}'
    fi
done < cai_operon_KPN15.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_J57.txt | sort | awk ' { a[i++]=$1; } END { print "J57\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "J57\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "J57\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "J57\tifp\t"$4}'
    fi
done < cai_operon_J57.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_SHEDEN.0723.00001.txt | sort | awk ' { a[i++]=$1; } END { print "SHEDEN.0723.00001\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "SHEDEN.0723.00001\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "SHEDEN.0723.00001\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "SHEDEN.0723.00001\tifp\t"$4}'
    fi
done < cai_operon_SHEDEN.0723.00001.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_SERMAR.0723.00023.txt | sort | awk ' { a[i++]=$1; } END { print "SERMAR.0723.00023\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "SERMAR.0723.00023\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "SERMAR.0723.00023\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "SERMAR.0723.00023\tifp\t"$4}'
    fi
done < cai_operon_SERMAR.0723.00023.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_ESCCOL.0723.02103.txt | sort | awk ' { a[i++]=$1; } END { print "ESCCOL.0723.02103\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "ESCCOL.0723.02103\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "ESCCOL.0723.02103\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "ESCCOL.0723.02103\tifp\t"$4}'
    fi
done < cai_operon_ESCCOL.0723.02103.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_SALENT.0723.00028.txt | sort | awk ' { a[i++]=$1; } END { print "SALENT.0723.00028\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "SALENT.0723.00028\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "SALENT.0723.00028\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "SALENT.0723.00028\tifp\t"$4}'
    fi
done < cai_operon_SALENT.0723.00028.txt
line_number=0

cut -d" " -f4 cai_ribosomal_proteins_ENTCLO.0723.00001.txt | sort | awk ' { a[i++]=$1; } END { print "ENTCLO.0723.00001\tRibosomal proteins\t"a[int(i/2)]; }'

line_number=0
while IFS= read -r line; do
    ((line_number++))

    if [ $line_number -eq 1 ]; then
        echo "$line" | awk '{print "ENTCLO.0723.00001\tlysR\t"$4}'
    elif [ $line_number -eq 2 ]; then
        echo "$line" | awk '{print "ENTCLO.0723.00001\tpfp\t"$4}'
    elif [ $line_number -eq 3 ]; then
        echo "$line" | awk '{print "ENTCLO.0723.00001\tifp\t"$4}'
    fi
done < cai_operon_ENTCLO.0723.00001.txt
line_number=0


