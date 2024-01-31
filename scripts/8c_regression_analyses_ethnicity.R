#### Regression analyses + interaction with ethnicity

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)
library(dplyr)
library(compositions)
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

## Regression analyses Menopause ####
pred_sex <- pred_sex %>% slice_head(n = 20)

ModelsexDutch <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'Dutch')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexDutch <- bind_rows(ModelsexDutch, coef)
}
ModelsexDutch$Pr...t.. <- p.adjust(ModelsexDutch$Pr...t.., method = "fdr", n = length(ModelsexDutch$Pr...t..))
ModelsexDutch$LWR <- ModelsexDutch$Estimate - 1.96*ModelsexDutch$Std..Error
ModelsexDutch$UPR <- ModelsexDutch$Estimate + 1.96*ModelsexDutch$Std..Error
ModelsexDutch <- ModelsexDutch %>%
  select(ASV=var_name, Dutch_Estimate = Estimate, Dutch_LWR = LWR, Dutch_UPR = UPR, Dutch_Pvalue=Pr...t..)
write.table(ModelsexDutch, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexDutch <- ModelsexDutch %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexSASurinamese <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'South-Asian Surinamese')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexSASurinamese <- bind_rows(ModelsexSASurinamese, coef)
}
ModelsexSASurinamese$Pr...t.. <- p.adjust(ModelsexSASurinamese$Pr...t.., method = "fdr", n = length(ModelsexSASurinamese$Pr...t..))
ModelsexSASurinamese$LWR <- ModelsexSASurinamese$Estimate - 1.96*ModelsexSASurinamese$Std..Error
ModelsexSASurinamese$UPR <- ModelsexSASurinamese$Estimate + 1.96*ModelsexSASurinamese$Std..Error
ModelsexSASurinamese <- ModelsexSASurinamese %>%
  select(ASV=var_name, SASur_Estimate = Estimate, SASur_LWR = LWR, SASur_UPR = UPR, SASur_Pvalue=Pr...t..)
write.table(ModelsexSASurinamese, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexSASurinamese <- ModelsexSASurinamese %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexASurinamese <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'African Surinamese')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexASurinamese <- bind_rows(ModelsexASurinamese, coef)
}
ModelsexASurinamese$Pr...t.. <- p.adjust(ModelsexASurinamese$Pr...t.., method = "fdr", n = length(ModelsexASurinamese$Pr...t..))
ModelsexASurinamese$LWR <- ModelsexASurinamese$Estimate - 1.96*ModelsexASurinamese$Std..Error
ModelsexASurinamese$UPR <- ModelsexASurinamese$Estimate + 1.96*ModelsexASurinamese$Std..Error
ModelsexASurinamese <- ModelsexASurinamese %>%
  select(ASV=var_name, ASur_Estimate = Estimate, ASur_LWR = LWR, ASur_UPR = UPR, ASur_Pvalue=Pr...t..)
write.table(ModelsexASurinamese, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexASurinamese <- ModelsexASurinamese %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexGhanaian <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'Ghanaian')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexGhanaian <- bind_rows(ModelsexGhanaian, coef)
}
ModelsexGhanaian$Pr...t.. <- p.adjust(ModelsexGhanaian$Pr...t.., method = "fdr", n = length(ModelsexGhanaian$Pr...t..))
ModelsexGhanaian$LWR <- ModelsexGhanaian$Estimate - 1.96*ModelsexGhanaian$Std..Error
ModelsexGhanaian$UPR <- ModelsexGhanaian$Estimate + 1.96*ModelsexGhanaian$Std..Error
ModelsexGhanaian <- ModelsexGhanaian %>%
  select(ASV=var_name, Ghanaian_Estimate = Estimate, Ghanaian_LWR = LWR, Ghanaian_UPR = UPR, Ghanaian_Pvalue=Pr...t..)
write.table(ModelsexGhanaian, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexGhanaian <- ModelsexGhanaian %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexTurkish <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'Turkish')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexTurkish <- bind_rows(ModelsexTurkish, coef)
}
ModelsexTurkish$Pr...t.. <- p.adjust(ModelsexTurkish$Pr...t.., method = "fdr", n = length(ModelsexTurkish$Pr...t..))
ModelsexTurkish$LWR <- ModelsexTurkish$Estimate - 1.96*ModelsexTurkish$Std..Error
ModelsexTurkish$UPR <- ModelsexTurkish$Estimate + 1.96*ModelsexTurkish$Std..Error
ModelsexTurkish <- ModelsexTurkish %>%
  select(ASV=var_name, Turkish_Estimate = Estimate, Turkish_LWR = LWR, Turkish_UPR = UPR, Turkish_Pvalue=Pr...t..)
write.table(ModelsexTurkish, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexTurkish <- ModelsexTurkish %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexMoroccan <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'Morroccan')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexMoroccan <- bind_rows(ModelsexMoroccan, coef)
}
ModelsexMoroccan$Pr...t.. <- p.adjust(ModelsexMoroccan$Pr...t.., method = "fdr", n = length(ModelsexMoroccan$Pr...t..))
ModelsexMoroccan$LWR <- ModelsexMoroccan$Estimate - 1.96*ModelsexMoroccan$Std..Error
ModelsexMoroccan$UPR <- ModelsexMoroccan$Estimate + 1.96*ModelsexMoroccan$Std..Error
ModelsexMoroccan <- ModelsexMoroccan %>%
  select(ASV=var_name, Moroccan_Estimate = Estimate, Moroccan_LWR = LWR, Moroccan_UPR = UPR, Moroccan_Pvalue=Pr...t..)
write.table(ModelsexMoroccan, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexMoroccan <- ModelsexMoroccan %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelsexOther <- data.frame()
for (i in 1:nrow(pred_sex)) {
  var_name <- pred_sex$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ Sex * relevel(Ethnicity, ref = 'Other')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelsexOther <- bind_rows(ModelsexOther, coef)
}
ModelsexOther$Pr...t.. <- p.adjust(ModelsexOther$Pr...t.., method = "fdr", n = length(ModelsexOther$Pr...t..))
ModelsexOther$LWR <- ModelsexOther$Estimate - 1.96*ModelsexOther$Std..Error
ModelsexOther$UPR <- ModelsexOther$Estimate + 1.96*ModelsexOther$Std..Error
ModelsexOther <- ModelsexOther %>%
  select(ASV=var_name, Other_Estimate = Estimate, Other_LWR = LWR, Other_UPR = UPR, Other_Pvalue=Pr...t..)
write.table(ModelsexOther, "clipboard", sep="\t", dec=",", col.names=NA)
ModelsexOther <- ModelsexOther %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

## Combine datasets
Modelsex <- ModelsexDutch %>%
  left_join(ModelsexSASurinamese, by = "ASV") %>%
  left_join(ModelsexASurinamese, by = "ASV") %>%
  left_join(ModelsexGhanaian, by = "ASV") %>%
  left_join(ModelsexTurkish, by = "ASV") %>%
  left_join(ModelsexMoroccan, by = "ASV") %>%
  left_join(ModelsexOther, by = "ASV")

## Change ASV names to taxonomy
Modelsex$ASV <- make.unique(tax$Tax[match(Modelsex$ASV, tax$ASV)])

## Write table
write.csv2(Modelsex, "results/Sex ethnicity.csv")


## Regression analyses Menopause ####
pred_meno <- pred_meno %>% slice_head(n = 20)

ModelmenopauseDutch <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'Dutch')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseDutch <- bind_rows(ModelmenopauseDutch, coef)
}
ModelmenopauseDutch$Pr...t.. <- p.adjust(ModelmenopauseDutch$Pr...t.., method = "fdr", n = length(ModelmenopauseDutch$Pr...t..))
ModelmenopauseDutch$LWR <- ModelmenopauseDutch$Estimate - 1.96*ModelmenopauseDutch$Std..Error
ModelmenopauseDutch$UPR <- ModelmenopauseDutch$Estimate + 1.96*ModelmenopauseDutch$Std..Error
ModelmenopauseDutch <- ModelmenopauseDutch %>%
  select(ASV=var_name, Dutch_Estimate = Estimate, Dutch_LWR = LWR, Dutch_UPR = UPR, Dutch_Pvalue=Pr...t..)
write.table(ModelmenopauseDutch, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseDutch <- ModelmenopauseDutch %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseSASurinamese <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'South-Asian Surinamese')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseSASurinamese <- bind_rows(ModelmenopauseSASurinamese, coef)
}
ModelmenopauseSASurinamese$Pr...t.. <- p.adjust(ModelmenopauseSASurinamese$Pr...t.., method = "fdr", n = length(ModelmenopauseSASurinamese$Pr...t..))
ModelmenopauseSASurinamese$LWR <- ModelmenopauseSASurinamese$Estimate - 1.96*ModelmenopauseSASurinamese$Std..Error
ModelmenopauseSASurinamese$UPR <- ModelmenopauseSASurinamese$Estimate + 1.96*ModelmenopauseSASurinamese$Std..Error
ModelmenopauseSASurinamese <- ModelmenopauseSASurinamese %>%
  select(ASV=var_name, SASur_Estimate = Estimate, SASur_LWR = LWR, SASur_UPR = UPR, SASur_Pvalue=Pr...t..)
write.table(ModelmenopauseSASurinamese, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseSASurinamese <- ModelmenopauseSASurinamese %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseASurinamese <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'African Surinamese')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseASurinamese <- bind_rows(ModelmenopauseASurinamese, coef)
}
ModelmenopauseASurinamese$Pr...t.. <- p.adjust(ModelmenopauseASurinamese$Pr...t.., method = "fdr", n = length(ModelmenopauseASurinamese$Pr...t..))
ModelmenopauseASurinamese$LWR <- ModelmenopauseASurinamese$Estimate - 1.96*ModelmenopauseASurinamese$Std..Error
ModelmenopauseASurinamese$UPR <- ModelmenopauseASurinamese$Estimate + 1.96*ModelmenopauseASurinamese$Std..Error
ModelmenopauseASurinamese <- ModelmenopauseASurinamese %>%
  select(ASV=var_name, ASur_Estimate = Estimate, ASur_LWR = LWR, ASur_UPR = UPR, ASur_Pvalue=Pr...t..)
write.table(ModelmenopauseASurinamese, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseASurinamese <- ModelmenopauseASurinamese %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseGhanaian <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'Ghanaian')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseGhanaian <- bind_rows(ModelmenopauseGhanaian, coef)
}
ModelmenopauseGhanaian$Pr...t.. <- p.adjust(ModelmenopauseGhanaian$Pr...t.., method = "fdr", n = length(ModelmenopauseGhanaian$Pr...t..))
ModelmenopauseGhanaian$LWR <- ModelmenopauseGhanaian$Estimate - 1.96*ModelmenopauseGhanaian$Std..Error
ModelmenopauseGhanaian$UPR <- ModelmenopauseGhanaian$Estimate + 1.96*ModelmenopauseGhanaian$Std..Error
ModelmenopauseGhanaian <- ModelmenopauseGhanaian %>%
  select(ASV=var_name, Ghanaian_Estimate = Estimate, Ghanaian_LWR = LWR, Ghanaian_UPR = UPR, Ghanaian_Pvalue=Pr...t..)
write.table(ModelmenopauseGhanaian, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseGhanaian <- ModelmenopauseGhanaian %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseTurkish <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'Turkish')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseTurkish <- bind_rows(ModelmenopauseTurkish, coef)
}
ModelmenopauseTurkish$Pr...t.. <- p.adjust(ModelmenopauseTurkish$Pr...t.., method = "fdr", n = length(ModelmenopauseTurkish$Pr...t..))
ModelmenopauseTurkish$LWR <- ModelmenopauseTurkish$Estimate - 1.96*ModelmenopauseTurkish$Std..Error
ModelmenopauseTurkish$UPR <- ModelmenopauseTurkish$Estimate + 1.96*ModelmenopauseTurkish$Std..Error
ModelmenopauseTurkish <- ModelmenopauseTurkish %>%
  select(ASV=var_name, Turkish_Estimate = Estimate, Turkish_LWR = LWR, Turkish_UPR = UPR, Turkish_Pvalue=Pr...t..)
write.table(ModelmenopauseTurkish, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseTurkish <- ModelmenopauseTurkish %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseMoroccan <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'Morroccan')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseMoroccan <- bind_rows(ModelmenopauseMoroccan, coef)
}
ModelmenopauseMoroccan$Pr...t.. <- p.adjust(ModelmenopauseMoroccan$Pr...t.., method = "fdr", n = length(ModelmenopauseMoroccan$Pr...t..))
ModelmenopauseMoroccan$LWR <- ModelmenopauseMoroccan$Estimate - 1.96*ModelmenopauseMoroccan$Std..Error
ModelmenopauseMoroccan$UPR <- ModelmenopauseMoroccan$Estimate + 1.96*ModelmenopauseMoroccan$Std..Error
ModelmenopauseMoroccan <- ModelmenopauseMoroccan %>%
  select(ASV=var_name, Moroccan_Estimate = Estimate, Moroccan_LWR = LWR, Moroccan_UPR = UPR, Moroccan_Pvalue=Pr...t..)
write.table(ModelmenopauseMoroccan, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseMoroccan <- ModelmenopauseMoroccan %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

ModelmenopauseOther <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn * relevel(Ethnicity, ref = 'Other')"))
  coef <- data.frame(var_name, t(summary(lm(formula, data = dfcomplete))$coef[2, c(1, 2, 4)]))
  ModelmenopauseOther <- bind_rows(ModelmenopauseOther, coef)
}
ModelmenopauseOther$Pr...t.. <- p.adjust(ModelmenopauseOther$Pr...t.., method = "fdr", n = length(ModelmenopauseOther$Pr...t..))
ModelmenopauseOther$LWR <- ModelmenopauseOther$Estimate - 1.96*ModelmenopauseOther$Std..Error
ModelmenopauseOther$UPR <- ModelmenopauseOther$Estimate + 1.96*ModelmenopauseOther$Std..Error
ModelmenopauseOther <- ModelmenopauseOther %>%
  select(ASV=var_name, Other_Estimate = Estimate, Other_LWR = LWR, Other_UPR = UPR, Other_Pvalue=Pr...t..)
write.table(ModelmenopauseOther, "clipboard", sep="\t", dec=",", col.names=NA)
ModelmenopauseOther <- ModelmenopauseOther %>% mutate(across(everything(.), ~trimws(.x, which = "both")))

## Combine datasets
Modelmenopause <- ModelmenopauseDutch %>%
  left_join(ModelmenopauseSASurinamese, by = "ASV") %>%
  left_join(ModelmenopauseASurinamese, by = "ASV") %>%
  left_join(ModelmenopauseGhanaian, by = "ASV") %>%
  left_join(ModelmenopauseTurkish, by = "ASV") %>%
  left_join(ModelmenopauseMoroccan, by = "ASV") %>%
  left_join(ModelmenopauseOther, by = "ASV")

## Change ASV names to taxonomy
Modelmenopause$ASV <- make.unique(tax$Tax[match(Modelmenopause$ASV, tax$ASV)])

## Write table
write.csv2(Modelmenopause, "results/Menopause ethnicity.csv")




