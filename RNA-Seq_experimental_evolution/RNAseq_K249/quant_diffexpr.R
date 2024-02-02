library(Rsubread)
library(DESeq2)
library(apeglm)
library(dplyr)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)
library(EnhancedVolcano)



# Quantification with featureCounts

bamfiles <- dir(".", ".bam$")
saffile <- dir(".", "saf$")

fc <- featureCounts(bamfiles, annot.ext=saffile, strandSpecific=2, isPairedEnd=TRUE, nthreads=15)

write.table(x=data.frame(fc$annotation[,c("GeneID","Chr", "Length")], fc$counts, stringsAsFactors=FALSE), file="featureCounts.tsv", quote=FALSE, sep="\t", row.names=FALSE)
head(fc$counts)
summary(fc$counts)
fc$stat

countdata <- fc$counts
colnames(countdata) <- gsub("\\.bam$", "", colnames(countdata))
nrow(countdata)

par(mar=c(8,4,4,1)+0.1)
barplot(colSums(countdata), las=3, main="Counts")



# DESeqDataSet object and experimental design

condition <- c('wt', 'wt', 'tc', 'tc', 'tc')
replicate <- c('1', '2', '1', '2', '3')
coldata <- as.data.frame(cbind(colnames(countdata), condition, replicate))
coldata

dds <- DESeqDataSetFromMatrix(countData = countdata, colData = coldata, design = ~ condition)
head(dds)
nrow(dds)

dds$condition <- relevel(dds$condition, ref = "wt")

# Pre-filtering the dataset

keep <- rowSums(counts(dds)) > 1
dds <- dds[keep,]
nrow(dds)




# Analysis of count data

rld <- rlog(dds, blind = FALSE)
head(assay(rld))

# Heatmap of distances between samples

sampleDists <- dist(t(assay(rld)))
sampleDistMatrix <- as.matrix(sampleDists)
rownames(sampleDistMatrix) <- paste(rld$condition, sep="-")
colnames(sampleDistMatrix) <- NULL
colors <- colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
pheatmap(sampleDistMatrix,
         clustering_distance_rows=sampleDists,
         clustering_distance_cols=sampleDists,
         col=colors,
         main="Heatmap of Euclidean distances")

# PCA plot

pcaData <- plotPCA(rld, intgroup = c("condition", "replicate"), returnData = TRUE)
percentVar <- round(100 * attr(pcaData, "percentVar"))
ggplot(pcaData, aes(x = PC1, y = PC2, color = condition, shape = replicate)) +
  geom_point(size =3) +
  xlab(paste0("PC1: ", percentVar[1], "% variance")) +
  ylab(paste0("PC2: ", percentVar[2], "% variance")) +
  coord_fixed() +
  theme_bw()




# Differential expression analysis

dds <- DESeq(dds)
resultsNames(dds)

res <- results(dds)
res <- lfcShrink(dds, coef="condition_tc_vs_wt", type="apeglm")
res
summary(res)

saftable <- read.csv(saffile, header=TRUE, sep="\t")
saftable <- as.data.frame(saftable)
resdata <- merge(as.data.frame(res), as.data.frame(counts(dds, normalized=TRUE)), by="row.names", sort=FALSE)
names(resdata)[1] <- "GeneID"
resdata <- merge(x=as.data.frame(resdata), saftable, by="GeneID")
head(resdata)
write.table(resdata, file="DE_results_raw.tsv", sep="\t", row.names=FALSE)




# MA plot

DESeq2::plotMA(res, ylim=c(-5,5), main="MA plot")

# Volcano plot

EnhancedVolcano(res, lab = rownames(res), x = 'log2FoldChange', y = 'padj', pCutoff=0.05, FCcutoff=0, subtitleLabSize=0, axisLabSize=10, labSize=3, legendLabSize=10, legendIconSize=3, col=c('black', 'black', 'gold2', 'blue'), colAlpha = 0.3, drawConnectors = TRUE)

