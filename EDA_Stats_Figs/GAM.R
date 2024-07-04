library(flux)
library(tidyverse)
library(mgcv)

curve_data <- read.table("curves.csv", header=TRUE)

# Calculate the AUC
AUC <- curve_data %>% 
  group_by(Strain, Genotype, Drug, Concentration, Replicate) %>% 
  group_modify(~ data.frame(AUC = flux::auc(.x$Time, .x$OD600))) %>%
  ungroup()

AUC_stats <- AUC %>%
  group_by(Strain, Genotype, Drug, Concentration) %>%
  summarise(
    median_AUC = median(AUC, na.rm = TRUE),
    se_AUC = sd(AUC, na.rm = TRUE) / sqrt(n())
  ) %>%
  ungroup()

AUC$Strain <- as.factor(AUC$Strain)
AUC$Genotype <- as.factor(AUC$Genotype)
AUC$Drug <- as.factor(AUC$Drug)


# Plot AUC values by concentration
ggplot(AUC_stats, aes(x = Concentration, y = median_AUC, color = interaction(Strain, Genotype, Drug), group = interaction(Strain, Genotype, Drug))) +
  geom_line() +
  geom_ribbon(aes(ymin = median_AUC - se_AUC, ymax = median_AUC + se_AUC, fill = interaction(Strain, Genotype, Drug)), color = NA, alpha = 0.2) +
  labs(x = "Concentration",
       y = "AUC",
       color = "Group",
       fill = "Group") +
  theme_bw()



# GAM for CF13

AUC_CF13 <- AUC[AUC$Strain == "CF13",]
AUC_CF13$Genotype_Drug <- interaction(AUC_CF13$Genotype, AUC_CF13$Drug)
AUC_CF13$Genotype <- relevel(AUC_CF13$Genotype, ref = "pOXA-48DLysR")
AUC_CF13$Drug <- relevel(AUC_CF13$Drug, ref = "QC")

gam_mod_CF13 <- gam(AUC ~ s(Concentration, k=5, by=Genotype_Drug) + Genotype * Drug, data = AUC_CF13, method = "REML")
summary.gam(gam_mod_CF13)

plot(gam_mod_CF13, residuals=TRUE, pch=1, all.terms=TRUE, pages=1, rug=TRUE, cex=1, shade=TRUE)

gam.check(gam_mod_CF13)
concurvity(gam_mod_CF13)


# GAM for KPN15

AUC_KPN15 <- AUC[AUC$Strain == "KPN15",]
AUC_KPN15$Genotype_Drug <- interaction(AUC_KPN15$Genotype, AUC_KPN15$Drug)
AUC_KPN15$Genotype <- relevel(AUC_KPN15$Genotype, ref = "pOXA-48DLysR")
AUC_KPN15$Drug <- relevel(AUC_KPN15$Drug, ref = "QC")

gam_mod_KPN15 <- gam(AUC ~ s(Concentration, k=5, by=Genotype_Drug) + Genotype * Drug, data = AUC_KPN15, method = "REML")
summary(gam_mod_KPN15)

plot(gam_mod_KPN15, residuals=TRUE, pch=1, all.terms=TRUE, pages=1, rug=TRUE, cex=1, shade=TRUE)

gam.check(gam_mod_KPN15)
concurvity(gam_mod_KPN15)

