#### Regression analyses + interaction with ethnicity

## Libraries
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(rio)
library(ggsci)

## Load data
# Clinical data and phyloseq
df <- readRDS("data/clinicaldata.RDS")
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
tab <- as.data.frame(t(as(phyloseq@otu_table, 'matrix')))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/taxtable_rarefied.RDS")

## Regression analyses sex ####
pred_sex <- rio::import("data/feature_importance_sex.txt") %>% slice(1:20)

Modelsex <- c()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  model <- lm(formula = paste0("log(", var_name, "+1) ~ Sex * Ethnicity + Age + BMI + 
              HT + DM + CurrSmoking"), data = dfcomplete)
  tidymodel <- tidy(model, conf.int = TRUE) %>% filter(str_detect(term, "SexFemale:"))
  print(tidymodel)
  coef <- cbind(var_name, tidymodel)
  Modelsex <- bind_rows(Modelsex, coef)
}
Modelsex$q.value <- p.adjust(Modelsex$p.value, method = "fdr")
Modelsex$var_name <- tax$Tax[match(Modelsex$var_name, tax$ASV)]
Modelsex %>% filter(q.value < 0.05)
write.csv2(Modelsex, "results/lm/interactions_ethnicity_sex.csv")

## Regression analyses menopause ####
pred_meno <- rio::import("data/feature_importance_menopause.txt") %>% slice(1:20)

Modelmenopause <- c()
for (i in 1:nrow(pred_meno)) {
    var_name <- pred_meno$FeatName[i]
    model <- lm(formula = paste0("log(", var_name, "+1) ~ MenopauseYn * Ethnicity + Age + BMI + 
              HT + DM + CurrSmoking"), data = dfcomplete)
    tidymodel <- tidy(model, conf.int = TRUE) %>% filter(str_detect(term, "MenopauseYnPostmenopausal:"))
    print(tidymodel)
    coef <- cbind(var_name, tidymodel)
    Modelmenopause <- bind_rows(Modelmenopause, coef)
}
Modelmenopause$q.value <- p.adjust(Modelmenopause$p.value, method = "fdr")
Modelmenopause$var_name <- tax$Tax[match(Modelmenopause$var_name, tax$ASV)]
Modelmenopause %>% filter(q.value < 0.05)
write.csv2(Modelmenopause, "results/lm/interactions_ethnicity_menopause.csv")
