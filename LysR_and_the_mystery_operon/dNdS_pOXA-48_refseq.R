
setwd("/media/jorge/acce4b86-3c48-4146-9cd1-ada046e8124b/rnaseq_synteny/dnds_analysis/variant_calling_snippy/summary")

library(xlsx)
library(tidyr)
library(dplyr)

# Load variants identified by Snippy in CDSs

mutations_poxa <- read.xlsx("annotation_vcall_summary.xlsx", sheetIndex = 2)

summary_types_of_muts <- mutations_poxa %>%
  select(gene, type_mut, change_nuc) %>%
  group_by(gene, type_mut, change_nuc) %>%
  summarise(Number = n())

summary_types_of_muts <- summary_types_of_muts %>%
  select(gene, type_mut) %>%
  group_by(gene, type_mut) %>%
  summarise(Number = n())

summary_types_of_muts <- as.data.frame(summary_types_of_muts)

###write.xlsx(summary_types_of_muts, "summary_types_of_muts.xlsx")

### Load table with summary of Non-synonymous and synonymous SNPs identified per gene

summary_dNdS <- read.xlsx("summary_types_of_muts.xlsx", sheetIndex = 3) 

### Perform permutation test

# With all the strains together

pvals <- c()
p_valors <- c()

for (g in unique(summary_dNdS$gene)) {
  #print(g)
  
  # Get diference between n muts per nucl in gene and the mean of the rest of the genes as the observed statistic
  
  gene_score <- subset(summary_dNdS, gene == g)$mut_per_nucl
  rest_scores <- subset(summary_dNdS, gene != g)$mut_per_nucl
  
  obs_stat <- abs(mean(gene_score) - mean(rest_scores)) # The observed stat is the difference in mutations falling in the gene vs the mean of mutations falling per gene in the plasmid
  
  # Set the conditions for the number of observations (in this case 1 per gene) and permutations (10000 in this case)
  
  num_permutations <- 10000
  num_observations <- num_observations <- length(gene_score)
  perm_stat <- matrix(nrow=num_permutations, ncol=1) # matrix 10000 obs x 1 measurement per observation 
  
  #print(obs_stat)
  
  # Perform permutations
  
  for (i in 1:num_permutations) {
    
    shuffled_group <- sample(c(gene_score, rest_scores)) # Random shuffle
    n <- floor(runif(1, min = 1, max = length(shuffled_group))) # Random index
    perm_gene_score <- shuffled_group[n] # Pick random observation
    perm_rest_score <- shuffled_group[-n] # Generate outgroup
    perm_stat[i, ] <- abs(mean(perm_gene_score) - mean(shuffled_group)) # the observed random stat in each permutation 
  }
  
  # Count the number of times that we get a more extreme value than the one estimated for the gene
  # i.e. probability of getting the n muts per base of the gene given the distribution of all the genes (p value)
  
  p_values <- sum(perm_stat >= obs_stat)/num_permutations
  #print(p_value_gene)
  
  p_valors <- c(p_valors, p_values) # Append p value to table
  
}

stat_results <- as.data.frame(cbind(unique(summary_dNdS$gene),p_valors)) # Append p values to genes

