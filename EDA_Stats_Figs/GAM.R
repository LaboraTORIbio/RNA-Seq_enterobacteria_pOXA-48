library(flux)
library(tidyverse)
library(mgcv)

curve_data <- read.table("curves.csv", header=TRUE)

# Calculate the AUC
AUC <- curve_data %>% 
  group_by(Strain, Genotype, Drug, Concentration, Replicate) %>% 
  group_modify(~ data.frame(AUC = flux::auc(.x$Time, .x$OD600))) %>%
  ungroup()
AUC_med <- AUC %>% 
  group_by(Strain, Genotype, Drug, Concentration) %>% 
  group_modify(~ data.frame(AUC = median(.x$AUC))) %>%
  ungroup()
AUC_med$Strain <- as.factor(AUC_med$Strain)
AUC_med$Genotype <- as.factor(AUC_med$Genotype)
AUC_med$Drug <- as.factor(AUC_med$Drug)

# Plot AUC values by concentration
ggplot(AUC_med, aes(x = Concentration, y = AUC, color = interaction(Strain, Genotype, Drug))) +
  geom_line() +
  labs(x = "Concentration",
       y = "AUC",
       color = "Group") +
  theme_bw()

# Fit GAM model
gam_mod <- gam(AUC ~ s(Concentration, k=5, by=Strain) + Strain * Genotype * Drug, data = AUC_med, method = "REML")
summary(gam_mod)
plot(gam_mod, residuals=TRUE, pch=1, all.terms=TRUE, pages=1, rug=TRUE, cex=1, shade=TRUE)
gam.check(gam_mod)
concurvity(gam_mod)
