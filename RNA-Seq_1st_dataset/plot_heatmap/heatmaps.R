
library(ggplot2)
library(cowplot)



# Load data

## DEGs with enriched biological processes
degs <- read.csv("DEGs4heatmap.tsv", sep="\t")
# Convert padj column to a factor with appropriate levels
degs$padj <- factor(degs$padj < 0.05, levels = c(FALSE, TRUE))
# Set colors for the heatmap
heatmap_colors <- c("gray70", "black")
## Operon genes
operon <- read.csv("operon_DEGs.csv", sep="\t")



# Heatmap of DEGs with enriched biological processes

degs$RefSeq <- factor(degs$RefSeq, levels = unique(degs$RefSeq[order(degs$Parental)]))
St_order <- c("C325", "EC10", "MG1655", "CF13", "H53", "J57", "K147", "KPN15", "KPN18", "KPN04", "KPN07", "KPN10")
gene_order <- c("matP", "mraY", "murA", "murC", "murE", "murF", "pal", "sdiA", "slmA", "ftsW", "mrdB", "mukB", "mukE", "sdiA", "bacA", "dacA", "dacD", "ftsI_WP_000642204", "ftsW", "ispU", "mltA", "mltD", "mrdA", "mtgA", "murI", "pbpC", "pbpG", "sltY", "yceG", "ypjD", "apbC", "iscA", "iscU", "sufA", "sufB", "sufE", "lpoB", "pbpC", "purC", "purD", "purE", "purF", "purL", "purM", "purN", "purT", "frdA", "glpC", "argA", "argE", "argG", "tyrR", "bglA", "glmM", "gph", "malZ", "nagA", "nagZ", "pgm", "aceE", "glk", "hisB", "hisC", "hisD", "hisF", "hisG", "hisIE", "ilvC", "ilvD", "ilvG", "arnA", "arnB", "arnT", "lpxC", "lpxH", "lpxK", "galF", "nirD", "citG", "cmk", "gmk", "pdxK", "pdxY", "pyk", "pykF", "thiD", "thrB", "fdnI", "fdoI", "hycE", "thiC", "thiI", "rplC", "rplD", "rplE", "rplF", "rplI", "rplJ", "rplK", "rplM", "rplN", "rplO", "rplP", "rplQ", "rplR", "rplS", "rplU", "rplW", "rplX", "rplY", "rpmA", "rpmB", "rpmC", "rpmD", "rpmE", "rpmG", "rpmI", "rpmJ", "rpsB", "rpsD", "rpsE", "rpsF", "rpsH", "rpsI", "rpsK", "rpsL", "rpsN", "rpsO", "rpsP", "rpsQ", "rpsR", "rpsS", "rpsT", "rpsU", "sra", "ybaK", "ykgO", "frdB", "mqo", "sdhB", "sucA", "sucD", "degP", "dsbB", "fklB", "fkpB", "grpE", "ppiA", "ppiB", "secB", "slyD", "asnC", "cra", "csiE", "cytR", "galS", "rapA", "rbsR", "rcsA", "sdiA", "yebC", "mglB", "lexA", "nfi", "nfo", "artM", "glnH", "glnP", "gltK", "hisQ", "plaP", "proY", "ydgI", "chbA", "chbB", "crr", "fruA", "fruB", "manX", "ptsH", "ptsI", "srlA", "srlB", "treB", "ulaB", "lptF", "metI", "plaP", "potB", "potC", "ugpA")
gene_order <- unique(gene_order)


# Log2FC heatmap
hm <- ggplot(degs, aes(x = factor(Strain, level=St_order), y = factor(Gene, level=rev(gene_order)), fill = log2FC)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#0000D5",
                       mid = "#FFFFFF",
                       high = "#D50000") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm')) +
  scale_y_discrete(labels = function(Gene) {
    italic_gene <- sapply(Gene, function(g) as.expression(substitute(italic(x), list(x = g))))
    return(italic_gene)
  }) +
  theme(axis.text.y = element_text(size = 3))


# Significance heatmap
padj <- ggplot(degs, aes(x = factor(Strain, level=St_order), y = factor(Gene, level=rev(gene_order)), fill = padj)) +
  geom_tile() +
  scale_fill_manual(values = heatmap_colors) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm')) +
  theme(axis.text.y = element_text(size = 3))



# Heatmap of expression of operon genes

operon_order <- c("isochorismatase family protein", "pirin family protein", "lysR")

# Log2FC heatmap
hmoperon <- ggplot(operon, aes(x = factor(Strain, level=St_order), y = factor(Gene, level=operon_order), fill = log2FC)) + 
  geom_tile() +
  scale_fill_gradient2(low = "#0000D5",
                       mid = "#FFFFFF",
                       high = "#D50000") +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm'))

# Significance heatmap
padjoperon <- ggplot(operon, aes(x = factor(Strain, level=St_order), y = factor(Gene, level=operon_order), fill = padj)) +
  geom_tile() +
  scale_fill_manual(values = heatmap_colors) +
  theme(axis.text.x = element_text(angle=90, vjust=0.5, hjust=1), legend.key.size = unit(0.3, 'cm'))



p <- plot_grid(hmoperon, padjoperon, hm, padj, align="v", nrow=2, ncol=2, rel_widths=c(2, 1.4), rel_heights=c(0.3, 2))



## Session Info

sessionInfo()

