#### Regression analyses + interaction with ethnicity

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)
library(dplyr)
library(ggsci)


## Load data
# Clinical data and phyloseq
df <- readRDS("data/clinicaldata.RDS")
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
tab <- as.data.frame(t(as(phyloseq@otu_table, 'matrix')))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/taxtable_rarefied.RDS")

# Best predictors
pred_meno <- rio::import("data/feature_importance_menopause.txt")
rownames(pred_meno) <- pred_meno$FeatName
pred_sex <- rio::import("data/feature_importance_sex.txt")
rownames(pred_sex) <- pred_sex$FeatName

## Regression analyses sex ####
pred_sex <- pred_sex %>% slice_head(n = 20)

Modelsex <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * Ethnicity"))
  coef <- data.frame(var_name, summary(lm(formula, data = dfcomplete))$coef[c(9,10,11,12,13,14), c(1, 2, 4)])
  Modelsex <- bind_rows(Modelsex, coef)
}
Modelsex <- as.data.frame(Modelsex)
Modelsex$Pr...t.. <- p.adjust(Modelsex$Pr...t.., method = "fdr", n = length(Modelsex$Pr...t..))
Modelsex$LWR <- Modelsex$Estimate - 1.96*Modelsex$Std..Error
Modelsex$UPR <- Modelsex$Estimate + 1.96*Modelsex$Std..Error
Modelsex <- Modelsex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
Modelsex$var_name <- tax$Tax[match(Modelsex$var_name, tax$ASV)]
write.csv2(Modelsex, "results/Sex ethnicity.csv")

## Regression analyses menopause ####
pred_meno <- pred_meno %>% slice_head(n = 20)

Modelmenopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * Ethnicity"))
  coef <- data.frame(var_name, summary(lm(formula, data = dfcomplete))$coef[c(9,10,11,12,13,14), c(1, 2, 4)])
  Modelmenopause <- bind_rows(Modelmenopause, coef)
}
Modelmenopause <- as.data.frame(Modelmenopause)
Modelmenopause$Pr...t.. <- p.adjust(Modelmenopause$Pr...t.., method = "fdr", n = length(Modelsex$Pr...t..))
Modelmenopause$LWR <- Modelmenopause$Estimate - 1.96*Modelmenopause$Std..Error
Modelmenopause$UPR <- Modelmenopause$Estimate + 1.96*Modelmenopause$Std..Error
Modelmenopause <- Modelmenopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
Modelmenopause$var_name <- tax$Tax[match(Modelmenopause$var_name, tax$ASV)]
write.csv2(Modelmenopause, "results/Menopause ethnicity.csv")





