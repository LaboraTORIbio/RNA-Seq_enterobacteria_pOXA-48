# Load tables
KPN15 <- read.csv("../RNAseq_KPN15/first_genes_reg_by_lysRs.txt", header=F, sep="\t")
CF13 <- read.csv("../RNAseq_CF13/first_genes_reg_by_lysRs.txt", header=F, sep="\t")

library(ggplot2)
library(cowplot)
library(dplyr)


# Plot KPN15
KPN15$V23 <- reorder(KPN15$V20, -KPN15$V3)
plot_KPN15 <- ggplot(KPN15, aes(x=V23, y=V3, color = factor(sign(V3)))) +
  geom_point(size=2, alpha=0.5) +
  coord_flip() +
  theme_classic() +
  ylim(-2, 6) +
  scale_color_manual(values = c("royalblue1", "tomato2")) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
plot_KPN15


# Plot CF13
CF13$V22 <- reorder(CF13$V20, -CF13$V3)
plot_CF13 <- ggplot(CF13, aes(x=V22, y=V3, color = factor(sign(V3)))) +
  geom_point(size=2, alpha=0.5) +
  coord_flip() +
  theme_classic() +
  ylim(-2, 6) +
  scale_color_manual(values = c("royalblue1", "tomato2")) +
  theme(axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7))
plot_CF13


p <- plot_grid(plot_CF13, plot_KPN15, align="v", nrow=2, rel_heights=c(1,1.51))
p
