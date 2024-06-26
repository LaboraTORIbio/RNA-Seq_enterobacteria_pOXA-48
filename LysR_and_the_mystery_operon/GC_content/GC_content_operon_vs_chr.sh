########################################################################################################################
#### Calculate the GC content of the chromosome and the operon system (from start of lysR to end of isochorismatase) ###
########################################################################################################################

# Code adapted from: https://www.biostars.org/p/70167/

# Copy the fasta sequences of the RefSeq genomes to ./GC_content

# Run from directory ./GC_content

# Get operon seqs:
samtools faidx ../../../Closed_sequences/CF13.fasta 1:1680069-1682643 > seq_operon_CF13.fasta
samtools faidx ../../../Closed_sequences/TC_KPN15.fasta 1:3247851-3250427 > seq_operon_TC_KPN15.fasta
samtools faidx ../../../Closed_sequences/J57.fasta 1:3534120-3536696 > seq_operon_J57.fasta
samtools faidx SHEDEN.0723.00001.fna NC_007954.1:4511528-4514108 > seq_operon_SHEDEN.0723.00001.fasta
samtools faidx SERMAR.0723.00023.fna NZ_CP020503.1:993982-996642 > seq_operon_SERMAR.0723.00023.fasta
samtools faidx ESCCOL.0723.02103.fna NZ_CP028747.1:869126-871650 > seq_operon_ESCCOL.0723.02103.fasta
samtools faidx SALENT.0723.00028.fna NZ_CP019186.1:2768837-2771429 > seq_operon_SALENT.0723.00028.fasta

samtools faidx ENTCLO.0723.00001.fna NC_014121.1:2647683-2650256 > seq_operon_ENTCLO.0723.00001.fasta



### Calculate GC content:

echo -e "Strain\tHeader\tFeature\tGC content"

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "CF13\t"h"\tOperon system\t"(cg/t) }' seq_operon_CF13.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "CF13\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "CF13\t"h"\tGenome\t"(cg/t); }' ../../../Closed_sequences/CF13.fasta

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "KPN15\t"h"\tOperon system\t"(cg/t) }' seq_operon_TC_KPN15.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "KPN15\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "KPN15\t"h"\tGenome\t"(cg/t); }' ../../../Closed_sequences/TC_KPN15.fasta

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "J57\t"h"\tOperon system\t"(cg/t) }' seq_operon_J57.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "J57\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "J57\t"h"\tGenome\t"(cg/t); }' ../../../Closed_sequences/J57.fasta

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SHEDEN.0723.00001\t"h"\tOperon system\t"(cg/t) }' seq_operon_SHEDEN.0723.00001.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "SHEDEN.0723.00001\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SHEDEN.0723.00001\t"h"\tGenome\t"(cg/t); }' SHEDEN.0723.00001.fna

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SERMAR.0723.00023\t"h"\tOperon system\t"(cg/t) }' seq_operon_SERMAR.0723.00023.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "SERMAR.0723.00023\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SERMAR.0723.00023\t"h"\tGenome\t"(cg/t); }' SERMAR.0723.00023.fna

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "ESCCOL.0723.02103\t"h"\tOperon system\t"(cg/t) }' seq_operon_ESCCOL.0723.02103.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "ESCCOL.0723.02103\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "ESCCOL.0723.02103\t"h"\tGenome\t"(cg/t); }' ESCCOL.0723.02103.fna

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SALENT.0723.00028\t"h"\tOperon system\t"(cg/t) }' seq_operon_SALENT.0723.00028.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "SALENT.0723.00028\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "SALENT.0723.00028\t"h"\tGenome\t"(cg/t); }' SALENT.0723.00028.fna

#-------------------------------------------

awk 'BEGIN { FS=""; h="" }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (h != "") print h"\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "ENTCLO.0723.00001\t"h"\tOperon system\t"(cg/t) }' seq_operon_ENTCLO.0723.00001.fasta

awk 'BEGIN { FS=""; }
{
    if ($1 != ">") {
        for (i = 1; i <= NF; i++) {
            if ($i ~ /[ACTGactg]/) t++;
            if ($i ~ /[CGcg]/) cg++;
        }
    } else {
        if (t > 0) print "ENTCLO.0723.00001\t"h"\tGenome\t"(cg/t);
        h = substr($0, 2);
        cg = t = 0;
    }
}
END { print "ENTCLO.0723.00001\t"h"\tGenome\t"(cg/t); }' ENTCLO.0723.00001.fna

