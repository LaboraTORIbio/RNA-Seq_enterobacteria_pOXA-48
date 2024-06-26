library(DESeq2)
library(ggfortify)
library(Rtsne)
library(umap)
library(tidyr)
library(dplyr)
library(pheatmap)


### Importing raw counts

counts_C325 <- read.table("../RNAseq_C325/featureCounts.tsv", header = TRUE)
counts_CF13 <- read.table("../RNAseq_CF13/featureCounts.tsv", header = TRUE)
counts_H53 <- read.table("../RNAseq_H53/featureCounts.tsv", header = TRUE)
counts_J57 <- read.table("../RNAseq_J57/featureCounts.tsv", header = TRUE)
counts_K147 <- read.table("../RNAseq_K147/featureCounts.tsv", header = TRUE)
counts_EC10 <- read.table("../RNAseq_EC10/featureCounts.tsv", header = TRUE)
counts_KPN04 <- read.table("../RNAseq_KPN04/featureCounts.tsv", header = TRUE)
counts_KPN07 <- read.table("../RNAseq_KPN07/featureCounts.tsv", header = TRUE)
counts_KPN10 <- read.table("../RNAseq_KPN10/featureCounts.tsv", header = TRUE)
counts_KPN15 <- read.table("../RNAseq_KPN15/featureCounts.tsv", header = TRUE)
counts_KPN18 <- read.table("../RNAseq_KPN18/featureCounts.tsv", header = TRUE)
counts_MG1655 <- read.table("../RNAseq_MG1655/featureCounts.tsv", header = TRUE)

colnames(counts_C325) <- gsub("\\.bam$", "", colnames(counts_C325))
colnames(counts_CF13) <- gsub("\\.bam$", "", colnames(counts_CF13))
colnames(counts_H53) <- gsub("\\.bam$", "", colnames(counts_H53))
colnames(counts_J57) <- gsub("\\.bam$", "", colnames(counts_J57))
colnames(counts_K147) <- gsub("\\.bam$", "", colnames(counts_K147))
colnames(counts_EC10) <- gsub("\\.bam$", "", colnames(counts_EC10))
colnames(counts_KPN04) <- gsub("\\.bam$", "", colnames(counts_KPN04))
colnames(counts_KPN07) <- gsub("\\.bam$", "", colnames(counts_KPN07))
colnames(counts_KPN10) <- gsub("\\.bam$", "", colnames(counts_KPN10))
colnames(counts_KPN15) <- gsub("\\.bam$", "", colnames(counts_KPN15))
colnames(counts_KPN18) <- gsub("\\.bam$", "", colnames(counts_KPN18))
colnames(counts_MG1655) <- gsub("\\.bam$", "", colnames(counts_MG1655))

rownames(counts_C325) <- counts_C325[,1]
counts_C325[,1] <- NULL
counts_C325$Chr <- NULL
counts_C325$Length <- NULL
rownames(counts_CF13) <- counts_CF13[,1]
counts_CF13[,1] <- NULL
counts_CF13$Chr <- NULL
counts_CF13$Length <- NULL
rownames(counts_H53) <- counts_H53[,1]
counts_H53[,1] <- NULL
counts_H53$Chr <- NULL
counts_H53$Length <- NULL
rownames(counts_J57) <- counts_J57[,1]
counts_J57[,1] <- NULL
counts_J57$Chr <- NULL
counts_J57$Length <- NULL
rownames(counts_K147) <- counts_K147[,1]
counts_K147[,1] <- NULL
counts_K147$Chr <- NULL
counts_K147$Length <- NULL
rownames(counts_EC10) <- counts_EC10[,1]
counts_EC10[,1] <- NULL
counts_EC10$Chr <- NULL
counts_EC10$Length <- NULL
rownames(counts_KPN04) <- counts_KPN04[,1]
counts_KPN04[,1] <- NULL
counts_KPN04$Chr <- NULL
counts_KPN04$Length <- NULL
rownames(counts_KPN07) <- counts_KPN07[,1]
counts_KPN07[,1] <- NULL
counts_KPN07$Chr <- NULL
counts_KPN07$Length <- NULL
rownames(counts_KPN10) <- counts_KPN10[,1]
counts_KPN10[,1] <- NULL
counts_KPN10$Chr <- NULL
counts_KPN10$Length <- NULL
rownames(counts_KPN15) <- counts_KPN15[,1]
counts_KPN15[,1] <- NULL
counts_KPN15$Chr <- NULL
counts_KPN15$Length <- NULL
rownames(counts_KPN18) <- counts_KPN18[,1]
counts_KPN18[,1] <- NULL
counts_KPN18$Chr <- NULL
counts_KPN18$Length <- NULL
rownames(counts_MG1655) <- counts_MG1655[,1]
counts_MG1655[,1] <- NULL
counts_MG1655$Chr <- NULL
counts_MG1655$Length <- NULL


###  Creating DESeqDataSet objects

# C325
condition_C325 <- c('wt', 'wt', 'wt', 'cur', 'cur', 'cur')
replicate_C325 <- c('1', '2', '3', '1', '2', '3')
coldata_C325 <- as.data.frame(cbind(colnames(counts_C325), condition_C325, replicate_C325))
dds_C325 <- DESeqDataSetFromMatrix(countData = counts_C325, colData = coldata_C325, design = ~ condition_C325)
dds_C325$condition_C325 <- relevel(dds_C325$condition_C325, ref = "cur")
# CF13
condition_CF13 <- c('wt', 'wt', 'wt', 'cur', 'cur')
replicate_CF13 <- c('1', '2', '3', '1', '2')
coldata_CF13 <- as.data.frame(cbind(colnames(counts_CF13), condition_CF13, replicate_CF13))
dds_CF13 <- DESeqDataSetFromMatrix(countData = counts_CF13, colData = coldata_CF13, design = ~ condition_CF13)
dds_CF13$condition_CF13 <- relevel(dds_CF13$condition_CF13, ref = "cur")
# H53
condition_H53 <- c('wt', 'wt', 'wt', 'cur', 'cur', 'cur')
replicate_H53 <- c('1', '2', '3', '1', '2', '3')
coldata_H53 <- as.data.frame(cbind(colnames(counts_H53), condition_H53, replicate_H53))
dds_H53 <- DESeqDataSetFromMatrix(countData = counts_H53, colData = coldata_H53, design = ~ condition_H53)
dds_H53$condition_H53 <- relevel(dds_H53$condition_H53, ref = "cur")
# J57
condition_J57 <- c('wt', 'wt', 'wt', 'cur', 'cur', 'cur')
replicate_J57 <- c('1', '2', '3', '1', '2', '3')
coldata_J57 <- as.data.frame(cbind(colnames(counts_J57), condition_J57, replicate_J57))
dds_J57 <- DESeqDataSetFromMatrix(countData = counts_J57, colData = coldata_J57, design = ~ condition_J57)
dds_J57$condition_J57 <- relevel(dds_J57$condition_J57, ref = "cur")
# K147
condition_K147 <- c('wt', 'wt', 'wt', 'wt', 'cur', 'cur')
replicate_K147 <- c('1', '2', '3', '4', '1', '2')
coldata_K147 <- as.data.frame(cbind(colnames(counts_K147), condition_K147, replicate_K147))
dds_K147 <- DESeqDataSetFromMatrix(countData = counts_K147, colData = coldata_K147, design = ~ condition_K147)
dds_K147$condition_K147 <- relevel(dds_K147$condition_K147, ref = "cur")
# EC10
condition_EC10 <- c('wt', 'wt', 'wt', 'tc', 'tc', 'tc')
replicate_EC10 <- c('1', '2', '3', '1', '2', '3')
coldata_EC10 <- as.data.frame(cbind(colnames(counts_EC10), condition_EC10, replicate_EC10))
dds_EC10 <- DESeqDataSetFromMatrix(countData = counts_EC10, colData = coldata_EC10, design = ~ condition_EC10)
dds_EC10$condition_EC10 <- relevel(dds_EC10$condition_EC10, ref = "tc")
# KPN04
condition_KPN04 <- c('wt', 'wt', 'wt', 'tc', 'tc', 'tc')
replicate_KPN04 <- c('1', '2', '3', '1', '2', '3')
coldata_KPN04 <- as.data.frame(cbind(colnames(counts_KPN04), condition_KPN04, replicate_KPN04))
dds_KPN04 <- DESeqDataSetFromMatrix(countData = counts_KPN04, colData = coldata_KPN04, design = ~ condition_KPN04)
dds_KPN04$condition_KPN04 <- relevel(dds_KPN04$condition_KPN04, ref = "tc")
# KPN07
condition_KPN07 <- c('wt', 'wt', 'wt', 'tc', 'tc', 'tc')
replicate_KPN07 <- c('1', '2', '3', '1', '2', '3')
coldata_KPN07 <- as.data.frame(cbind(colnames(counts_KPN07), condition_KPN07, replicate_KPN07))
dds_KPN07 <- DESeqDataSetFromMatrix(countData = counts_KPN07, colData = coldata_KPN07, design = ~ condition_KPN07)
dds_KPN07$condition_KPN07 <- relevel(dds_KPN07$condition_KPN07, ref = "tc")
# KPN10
condition_KPN10 <- c('wt', 'wt', 'tc', 'tc')
replicate_KPN10 <- c('1', '2', '1', '2')
coldata_KPN10 <- as.data.frame(cbind(colnames(counts_KPN10), condition_KPN10, replicate_KPN10))
dds_KPN10 <- DESeqDataSetFromMatrix(countData = counts_KPN10, colData = coldata_KPN10, design = ~ condition_KPN10)
dds_KPN10$condition_KPN10 <- relevel(dds_KPN10$condition_KPN10, ref = "tc")
# KPN15
condition_KPN15 <- c('wt', 'wt', 'tc', 'tc', 'tc')
replicate_KPN15 <- c('1', '2', '1', '2', '3')
coldata_KPN15 <- as.data.frame(cbind(colnames(counts_KPN15), condition_KPN15, replicate_KPN15))
dds_KPN15 <- DESeqDataSetFromMatrix(countData = counts_KPN15, colData = coldata_KPN15, design = ~ condition_KPN15)
dds_KPN15$condition_KPN15 <- relevel(dds_KPN15$condition_KPN15, ref = "tc")
# KPN18
condition_KPN18 <- c('wt', 'wt', 'wt', 'tc', 'tc', 'tc')
replicate_KPN18 <- c('1', '2', '3', '1', '2', '3')
coldata_KPN18 <- as.data.frame(cbind(colnames(counts_KPN18), condition_KPN18, replicate_KPN18))
dds_KPN18 <- DESeqDataSetFromMatrix(countData = counts_KPN18, colData = coldata_KPN18, design = ~ condition_KPN18)
dds_KPN18$condition_KPN18 <- relevel(dds_KPN18$condition_KPN18, ref = "tc")
# MG1655
condition_MG1655 <- c('wt', 'wt', 'tc', 'tc', 'tc')
replicate_MG1655 <- c('1', '2', '1', '2', '3')
coldata_MG1655 <- as.data.frame(cbind(colnames(counts_MG1655), condition_MG1655, replicate_MG1655))
dds_MG1655 <- DESeqDataSetFromMatrix(countData = counts_MG1655, colData = coldata_MG1655, design = ~ condition_MG1655)
dds_MG1655$condition_MG1655 <- relevel(dds_MG1655$condition_MG1655, ref = "tc")


# Normalization with VST

# C325
vst_C325 <- vst(dds_C325, blind = FALSE)
vst_C325_out <- assay(vst_C325)
vst_C325_out <- cbind(GeneID = rownames(vst_C325_out), vst_C325_out)
write.table(x=vst_C325_out, file="vst_C325.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# CF13
vst_CF13 <- vst(dds_CF13, blind = FALSE)
vst_CF13_out <- assay(vst_CF13)
vst_CF13_out <- cbind(GeneID = rownames(vst_CF13_out), vst_CF13_out)
write.table(x=vst_CF13_out, file="vst_CF13.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# H53
vst_H53 <- vst(dds_H53, blind = FALSE)
vst_H53_out <- assay(vst_H53)
vst_H53_out <- cbind(GeneID = rownames(vst_H53_out), vst_H53_out)
write.table(x=vst_H53_out, file="vst_H53.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# J57
vst_J57 <- vst(dds_J57, blind = FALSE)
vst_J57_out <- assay(vst_J57)
vst_J57_out <- cbind(GeneID = rownames(vst_J57_out), vst_J57_out)
write.table(x=vst_J57_out, file="vst_J57.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# K147
vst_K147 <- vst(dds_K147, blind = FALSE)
vst_K147_out <- assay(vst_K147)
vst_K147_out <- cbind(GeneID = rownames(vst_K147_out), vst_K147_out)
write.table(x=vst_K147_out, file="vst_K147.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# EC10
vst_EC10 <- vst(dds_EC10, blind = FALSE)
vst_EC10_out <- assay(vst_EC10)
vst_EC10_out <- cbind(GeneID = rownames(vst_EC10_out), vst_EC10_out)
write.table(x=vst_EC10_out, file="vst_EC10.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# KPN04
vst_KPN04 <- vst(dds_KPN04, blind = FALSE)
vst_KPN04_out <- assay(vst_KPN04)
vst_KPN04_out <- cbind(GeneID = rownames(vst_KPN04_out), vst_KPN04_out)
write.table(x=vst_KPN04_out, file="vst_KPN04.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# KPN07
vst_KPN07 <- vst(dds_KPN07, blind = FALSE)
vst_KPN07_out <- assay(vst_KPN07)
vst_KPN07_out <- cbind(GeneID = rownames(vst_KPN07_out), vst_KPN07_out)
write.table(x=vst_KPN07_out, file="vst_KPN07.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# KPN10
vst_KPN10 <- vst(dds_KPN10, blind = FALSE)
vst_KPN10_out <- assay(vst_KPN10)
vst_KPN10_out <- cbind(GeneID = rownames(vst_KPN10_out), vst_KPN10_out)
write.table(x=vst_KPN10_out, file="vst_KPN10.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# KPN15
vst_KPN15 <- vst(dds_KPN15, blind = FALSE)
vst_KPN15_out <- assay(vst_KPN15)
vst_KPN15_out <- cbind(GeneID = rownames(vst_KPN15_out), vst_KPN15_out)
write.table(x=vst_KPN15_out, file="vst_KPN15.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# KPN18
vst_KPN18 <- vst(dds_KPN18, blind = FALSE)
vst_KPN18_out <- assay(vst_KPN18)
vst_KPN18_out <- cbind(GeneID = rownames(vst_KPN18_out), vst_KPN18_out)
write.table(x=vst_KPN18_out, file="vst_KPN18.tsv", quote=FALSE, sep="\t", row.names=FALSE)
# MG1655
vst_MG1655 <- vst(dds_MG1655, blind = FALSE)
vst_MG1655_out <- assay(vst_MG1655)
vst_MG1655_out <- cbind(GeneID = rownames(vst_MG1655_out), vst_MG1655_out)
write.table(x=vst_MG1655_out, file="vst_MG1655.tsv", quote=FALSE, sep="\t", row.names=FALSE)


### Run the script get_SCO_mappings.py to get a table of the single copy orthologs between the 12 strains
### and get_SCO_norm_counts.py to map those SCOs to the GeneIDs and get the normalized counts of the SCOs

color_palette <- c("indianred2", "lightsalmon", "gray95", "darkgoldenrod1", "#A5D6A7",
                   "#80CBC4", "lightblue1", "lightskyblue", "skyblue3", "plum2", "#BA68C8", "mediumorchid4")
set.seed(123)


### pOXA-48-carrying strains

strain_order <- c("C325", "TC_EC10", "MG1655p", "CF13", "J57", "H53", "K147", "TC_KPN15", "TC_KPN18", "TC_KPN04", "TC_KPN07", "TC_KPN10")

SCO_norm_counts_PC <- read.table("SCO_norm_counts_PC.tsv", header = TRUE)
rownames(SCO_norm_counts_PC) <- SCO_norm_counts_PC$SCO
SCO_norm_counts_PC <- SCO_norm_counts_PC[, -which(names(SCO_norm_counts_PC) == "SCO")]
SCO_norm_counts_PC_t <- t(SCO_norm_counts_PC)
SCO_norm_counts_PC_t <- as.data.frame(SCO_norm_counts_PC_t)
SCO_norm_counts_PC_t$Strain <- c(rep("C325",3), rep("CF13",3), rep("TC_EC10",3), rep("H53",3), rep("J57",3), rep("K147",4), 
                            rep("TC_KPN04",3), rep("TC_KPN07",3), rep("TC_KPN10",2), rep("TC_KPN15",3), rep("TC_KPN18",3), rep("MG1655p",3))
SCO_norm_counts_PC_t$Species <- c(rep("E_coli",3), rep("C_freundii",3), rep("E_coli",3), rep("K_pneumoniae",3),
                               rep("K_variicola",3), rep("K_pneumoniae",18), rep("E_coli",3))
SCO_norm_counts_PC_t$Strain <- factor(SCO_norm_counts_PC_t$Strain, levels = strain_order)

# PCA
pca_SCO_PC <- prcomp(SCO_norm_counts_PC_t[,c(1:2488)])

autoplot(pca_SCO_PC, data = SCO_norm_counts_PC_t, colour = "Strain", shape = "Species") + 
  theme_bw() +
  geom_point(aes(color = Strain, shape = factor(Species)), size = 3) +
  scale_color_manual(values = color_palette)

# t-SNE
tsne_result_PC <- Rtsne(SCO_norm_counts_PC_t[,c(1:2488)], perplexity=5)
tsne_df_PC <- data.frame(tsne_result_PC$Y, Strain = SCO_norm_counts_PC_t$Strain, Species = SCO_norm_counts_PC_t$Species)

ggplot(tsne_df_PC, aes(x = X1, y = X2, color = Strain)) +
  geom_point(aes(shape=factor(Species)), size=3) + labs(x = "t-SNE 1", y = "t-SNE 2") + theme_bw() +
  scale_color_manual(values = color_palette)


### pOXA-48-free strains

strain_order2 <- c("C325c1", "PF_EC10", "MG1655", "CF13c1", "J57c1", "H53c1", "K147c1", "PF_KPN15", "PF_KPN18", "PF_KPN04", "PF_KPN07", "PF_KPN10")

SCO_norm_counts_PF <- read.table("SCO_norm_counts_PF.tsv", header = TRUE)
rownames(SCO_norm_counts_PF) <- SCO_norm_counts_PF$SCO
SCO_norm_counts_PF <- SCO_norm_counts_PF[, -which(names(SCO_norm_counts_PF) == "SCO")]
SCO_norm_counts_PF_t <- t(SCO_norm_counts_PF)
SCO_norm_counts_PF_t <- as.data.frame(SCO_norm_counts_PF_t)
SCO_norm_counts_PF_t$Strain <- c(rep("C325c1",3), rep("CF13c1",2), rep("PF_EC10",3), rep("H53c1",3), rep("J57c1",3), rep("K147c1",2), 
                                 rep("PF_KPN04",3), rep("PF_KPN07",3), rep("PF_KPN10",2), rep("PF_KPN15",2), rep("PF_KPN18",3), rep("MG1655",2))
SCO_norm_counts_PF_t$Species <- c(rep("E_coli",3), rep("C_freundii",2), rep("E_coli",3), rep("K_pneumoniae",3),
                                  rep("K_variicola",3), rep("K_pneumoniae",15), rep("E_coli",2))
SCO_norm_counts_PF_t$Strain <- factor(SCO_norm_counts_PF_t$Strain, levels = strain_order2)

# PCA
pca_SCO_PF <- prcomp(SCO_norm_counts_PF_t[,c(1:2488)])

autoplot(pca_SCO_PF, data = SCO_norm_counts_PF_t, colour = "Strain", shape = "Species") + 
  theme_bw() +
  geom_point(aes(color = Strain, shape = factor(Species)), size = 3) +
  scale_color_manual(values = color_palette)

# t-SNE
tsne_result_PF <- Rtsne(SCO_norm_counts_PF_t[,c(1:2488)], perplexity=5)
tsne_df_PF <- data.frame(tsne_result_PF$Y, Strain = SCO_norm_counts_PF_t$Strain, Species = SCO_norm_counts_PF_t$Species)

ggplot(tsne_df_PF, aes(x = X1, y = X2, color = Strain)) +
  geom_point(aes(shape=factor(Species)), size=3) + labs(x = "t-SNE 1", y = "t-SNE 2") + theme_bw() +
  scale_color_manual(values = color_palette)


### All together

color_palette2 <- c("indianred2", "lightsalmon", "gray95", "darkgoldenrod1", "#A5D6A7",
                   "#80CBC4", "lightblue1", "lightskyblue", "skyblue3", "plum2", "#BA68C8", "mediumorchid4",
                   "indianred2", "lightsalmon", "gray95", "darkgoldenrod1", "#A5D6A7",
                   "#80CBC4", "lightblue1", "lightskyblue", "skyblue3", "plum2", "#BA68C8", "mediumorchid4")

SCO_norm_counts_PC_t$Type <- c(rep("pOXA-48-carrying", 36))
SCO_norm_counts_PC_t$Strain_short <- c(rep("C325",3), rep("CF13",3), rep("EC10",3), rep("H53",3), rep("J57",3), rep("K147",4), 
                                  rep("KPN04",3), rep("KPN07",3), rep("KPN10",2), rep("KPN15",3), rep("KPN18",3), rep("MG1655",3))
SCO_norm_counts_PF_t$Type <- c(rep("pOXA-48-free", 31))
SCO_norm_counts_PF_t$Strain_short <- c(rep("C325",3), rep("CF13",2), rep("EC10",3), rep("H53",3), rep("J57",3), rep("K147",2), 
                                  rep("KPN04",3), rep("KPN07",3), rep("KPN10",2), rep("KPN15",2), rep("KPN18",3), rep("MG1655",2))

all_data <- rbind(SCO_norm_counts_PC_t, SCO_norm_counts_PF_t)

# PCA
pca_SCO_all <- prcomp(all_data[,c(1:2488)])

autoplot(pca_SCO_all, data = all_data, colour = "Strain", shape = "Species") + 
  theme_bw() +
  geom_point(aes(color = Strain, shape = factor(Species)), size = 3) +
  scale_color_manual(values = color_palette2) +
  geom_point(shape=21, size=5, aes(color=Strain, fill=Type)) +
  scale_fill_manual(values = c("black", "transparent"))

# t-SNE
tsne_result_all <- Rtsne(all_data[,c(1:2488)], perplexity=5)
tsne_df_all <- data.frame(tsne_result_all$Y, Strain = all_data$Strain, Species = all_data$Species, Type = all_data$Type)

ggplot(tsne_df_all, aes(x = X1, y = X2, color = Strain)) +
  geom_point(aes(shape=factor(Species)), size=3) + labs(x = "t-SNE 1", y = "t-SNE 2") + theme_bw() +
  scale_color_manual(values = color_palette2)+
  geom_point(data = subset(tsne_df_all, Type == "pOXA-48-carrying"), shape=21, size=5)




### Run the script get_SCO_log2FC.py to get a table of the log2FC and padj values of the SCOs

SCO_logs <- read.table("SCO_log2FC.tsv", header = TRUE)

count_per_sco <- SCO_logs %>%
  group_by(SCO) %>%
  summarise(unique_strains = n_distinct(Strain))

df_filtered <- SCO_logs %>%
  inner_join(count_per_sco, by = "SCO") %>%
  filter(unique_strains == 12) %>%
  select(-unique_strains)  # Remove the count column

df_subset1 <- df_filtered[, c("Strain", "SCO", "log2FoldChange")]
df_subset1 <- distinct(df_subset1)
log2FC_mat <- pivot_wider(df_subset1, names_from = Strain, values_from = log2FoldChange)
log2FC_mat <- log2FC_mat[, c(2:13)]

df_subset2 <- df_filtered[, c("Strain", "SCO", "padj")]
df_subset2 <- distinct(df_subset2)
padj_mat <- pivot_wider(df_subset2, names_from = Strain, values_from = padj)
padj_mat <- padj_mat[, c(2:13)]
padj_mat_logical <- sapply(padj_mat, as.logical)


### Heatmap with hierarchical clustering

heatmap_result <- pheatmap(log2FC_mat, 
         cluster_rows = TRUE,
         cluster_cols = TRUE,
         show_rownames = FALSE,
         color = colorRampPalette(c("#0000D5", "#FFFFFF", "#D50000"))(100),
         return_dendrogram = TRUE
)


### Heatmap of padj values

ordered_padj_mat_logical <- padj_mat_logical[heatmap_result$tree_row$order, heatmap_result$tree_col$order]
ordered_padj_mat_numeric <- 1*ordered_padj_mat_logical

pheatmap(ordered_padj_mat_numeric, 
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         show_rownames = FALSE,
         color = c("gray70", "black")
)

