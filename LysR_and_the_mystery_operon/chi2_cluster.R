strains <- rbind(c(145, 28), c(4607, 3369))
rownames(strains) <- c("pOXA48", "no_pOXA48")
colnames(strains) <- c("with_cluster", "no_cluster")
strains
chisq.test(strains)
