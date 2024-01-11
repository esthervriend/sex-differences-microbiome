#### Regression analyses Menopause

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)

## Load data
# Clinical data and phyloseq
df <- readRDS("data/clinicaldata.RDS")
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
tab <- as.data.frame(t(as(phyloseq@otu_table, 'matrix')))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")

# Best predictors
pred_meno <- rio::import("data/feature_importance_menopause.txt")
rownames(pred_meno) <- pred_meno$FeatName


## Regression analyses Menopause
# Model 1 (Crude)
pred_meno <- pred_meno %>% slice_head(n = 20)

Model1Menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model1Menopause<- bind_rows(Model1Menopause,coef)
}
Model1Menopause$LWR <- Model1Menopause$Estimate - 1.96*Model1Menopause$Std..Error
Model1Menopause$UPR <- Model1Menopause$Estimate + 1.96*Model1Menopause$Std..Error
Model1Menopause <- Model1Menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model1Menopause[c("Estimate", "LWR", "UPR")] <- exp(Model1Menopause[c("Estimate", "LWR", "UPR")])
write.table(Model1Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model1Menopause <- Model1Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model1Menopause, "results/Model1Menopause.csv")

# Model 2 (+ Age)
Model2Menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + Age"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model2Menopause <- bind_rows(Model2Menopause,coef)
}
Model2Menopause$LWR <- Model2Menopause$Estimate - 1.96*Model2Menopause$Std..Error
Model2Menopause$UPR <- Model2Menopause$Estimate + 1.96*Model2Menopause$Std..Error
Model2Menopause <- Model2Menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model2Menopause[c("Estimate", "LWR", "UPR")] <- exp(Model2Menopause[c("Estimate", "LWR", "UPR")])
write.table(Model2Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model2Menopause <- Model2Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model2Menopause, "results/Model2Menopause.csv")


# Model 3 (+ Med)
Model3Menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + AntiHT + DMMed + LipidLowering + SystSteroids"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model3Menopause <- bind_rows(Model3Menopause,coef)
}
Model3Menopause$LWR <- Model3Menopause$Estimate - 1.96*Model3Menopause$Std..Error
Model3Menopause$UPR <- Model3Menopause$Estimate + 1.96*Model3Menopause$Std..Error
Model3Menopause <- Model3Menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model3Menopause[c("Estimate", "LWR", "UPR")] <- exp(Model3Menopause[c("Estimate", "LWR", "UPR")])
write.table(Model3Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model3Menopause <- Model3Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model3Menopause, "results/Model3Menopause.csv")


# Model 4 (+ CVD)
Model4Menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + BMI + HT + DM + CurrSmoking"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model4Menopause <- bind_rows(Model4Menopause,coef)
}
Model4Menopause$LWR <- Model4Menopause$Estimate - 1.96*Model4Menopause$Std..Error
Model4Menopause$UPR <- Model4Menopause$Estimate + 1.96*Model4Menopause$Std..Error
Model4Menopause <- Model4Menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model4Menopause[c("Estimate", "LWR", "UPR")] <- exp(Model4Menopause[c("Estimate", "LWR", "UPR")])
write.table(Model4Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model4Menopause <- Model4Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model4Menopause, "results/Model4Menopause.csv")



# Model 5 (+ dieet (later met zout))
Model5Menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + Alcohol + TotalCalories + Fibre"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model5Menopause <- bind_rows(Model5Menopause,coef)
}
Model5Menopause$LWR <- Model5Menopause$Estimate - 1.96*Model5Menopause$Std..Error
Model5Menopause$UPR <- Model5Menopause$Estimate + 1.96*Model5Menopause$Std..Error
Model5Menopause <- Model5Menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model5Menopause[c("Estimate", "LWR", "UPR")] <- exp(Model5Menopause[c("Estimate", "LWR", "UPR")])
write.table(Model5Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model5Menopause <- Model5Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model5Menopause, "results/Model5Menopause.csv")


