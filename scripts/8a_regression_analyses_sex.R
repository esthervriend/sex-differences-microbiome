#### Regression analyses Sex

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
pred_sex <- rio::import("data/feature_importance_sex.txt")
rownames(pred_sex) <- pred_sex$FeatName


## Regression analyses sex
# Model 1 (Crude)
pred_sex <- pred_sex %>% slice_head(n = 20)

Model1sex <- data.frame()
for (i in 1:nrow(pred_sex)) {
 var_name <- pred_sex$FeatName[i]
 formula <- as.formula(paste("log(", var_name, "+1) ~ Sex"))
 coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
 Model1sex <- bind_rows(Model1sex,coef)
}
Model1sex$LWR <- Model1sex$Estimate - 1.96*Model1sex$Std..Error
Model1sex$UPR <- Model1sex$Estimate + 1.96*Model1sex$Std..Error
Model1sex <- Model1sex %>%
    select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model1sex[c("Estimate", "LWR", "UPR")] <- exp(Model1sex[c("Estimate", "LWR", "UPR")])
write.table(Model1sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model1sex <- Model1sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model1sex, "results/Model1sex.csv")

# Model 2 (+ Age)
Model2sex <- data.frame()
for (i in 1:nrow(pred_sex)) {
    var_name <- pred_sex$FeatName[i]
    formula <- as.formula(paste("log(", var_name, "+1) ~ Sex + Age"))
    coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
    Model2sex <- bind_rows(Model2sex,coef)
}
Model2sex$LWR <- Model2sex$Estimate - 1.96*Model2sex$Std..Error
Model2sex$UPR <- Model2sex$Estimate + 1.96*Model2sex$Std..Error
Model2sex <- Model2sex %>%
    select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model2sex[c("Estimate", "LWR", "UPR")] <- exp(Model2sex[c("Estimate", "LWR", "UPR")])
write.table(Model2sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model2sex <- Model2sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model2sex, "results/Model2sex.csv")


# Model 3 (+ Med)
Model3sex <- data.frame()
for (i in 1:nrow(pred_sex)) {
    var_name <- pred_sex$FeatName[i]
    formula <- as.formula(paste("log(", var_name, "+1) ~ Sex + AntiHT + DMMed + LipidLowering + SystSteroids"))
    coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
    Model3sex <- bind_rows(Model3sex,coef)
}
Model3sex$LWR <- Model3sex$Estimate - 1.96*Model3sex$Std..Error
Model3sex$UPR <- Model3sex$Estimate + 1.96*Model3sex$Std..Error
Model3sex <- Model3sex %>%
    select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model3sex[c("Estimate", "LWR", "UPR")] <- exp(Model3sex[c("Estimate", "LWR", "UPR")])
write.table(Model3sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model3sex <- Model3sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model3sex, "results/Model3sex.csv")


# Model 4 (+ CVD)
Model4sex <- data.frame()
for (i in 1:nrow(pred_sex)) {
    var_name <- pred_sex$FeatName[i]
    formula <- as.formula(paste("log(", var_name, "+1) ~ Sex + BMI + HT + DM + CurrSmoking"))
    coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
    Model4sex <- bind_rows(Model4sex,coef)
}
Model4sex$LWR <- Model4sex$Estimate - 1.96*Model4sex$Std..Error
Model4sex$UPR <- Model4sex$Estimate + 1.96*Model4sex$Std..Error
Model4sex <- Model4sex %>%
    select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model4sex[c("Estimate", "LWR", "UPR")] <- exp(Model4sex[c("Estimate", "LWR", "UPR")])
write.table(Model4sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model4sex <- Model4sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model4sex, "results/Model4sex.csv")



# Model 5 (+ dieet (later met zout))
Model5sex <- data.frame()
for (i in 1:nrow(pred_sex)) {
    var_name <- pred_sex$FeatName[i]
    formula <- as.formula(paste("log(", var_name, "+1) ~ Sex + Alcohol + TotalCalories + Fibre"))
    coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
    Model5sex <- bind_rows(Model5sex,coef)
}
Model5sex$LWR <- Model5sex$Estimate - 1.96*Model5sex$Std..Error
Model5sex$UPR <- Model5sex$Estimate + 1.96*Model5sex$Std..Error
Model5sex <- Model5sex %>%
    select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model5sex[c("Estimate", "LWR", "UPR")] <- exp(Model5sex[c("Estimate", "LWR", "UPR")])
write.table(Model5sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model5sex <- Model5sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model5sex, "results/Model5sex.csv")


