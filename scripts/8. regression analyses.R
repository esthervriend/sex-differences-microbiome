#### Regression analyses

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
regsex <- rio::import("data/feature_importance_sex.txt")
rownames(regsex) <- regsex$FeatName
regmenopause <- rio::import("data/feature_importance_menopause.txt")


## Regression analyses sex
# Model 1
regsex <- filter(regsex, !FeatName %in% c("random_variable1", "random_variable2"))
Model1sex <- data.frame()
for (i in 1:nrow(regsex)) {
 var_name <- regsex$FeatName[i]
 formula <- as.formula(paste("Sex ~", var_name))
 coef <- data.frame(var_name,t(summary(glm(formula, family = "binomial",data = dfcomplete))$coef[2,c(1,2,4)]))
 Model1sex <- bind_rows(Model1sex,coef)
}
Model1sex$LWR <- Model1sex$Estimate - 1.96*Model1sex$Std..Error
Model1sex$UPR <- Model1sex$Estimate + 1.96*Model1sex$Std..Error

# Model 2
Model2sex <- data.frame()
for (i in 1:nrow(regsex)) {
    var_name <- regsex$FeatName[i]
    formula <- as.formula(paste("Sex ~ Age + HT + DM + Smoking + ", var_name))
    coef <- data.frame(var_name,t(summary(glm(formula, family = "binomial",data = dfcomplete))$coef[2,c(1,2,4)]))
    Model2sex <- bind_rows(Model2sex,coef)
}
Model2sex$LWR <- Model2sex$Estimate - 1.96*Model2sex$Std..Error
Model2sex$UPR <- Model2sex$Estimate + 1.96*Model2sex$Std..Error


## Regression analyses menopause
regmenopause <- filter(regmenopause, !FeatName %in% c("random_variable1", "random_variable2"))
Model1menopause <- data.frame()
for (i in 1:nrow(regmenopause)) {
    var_name <- regmenopause$FeatName[i]
    formula <- as.formula(paste("MenopauseYN ~", var_name))
    coef <- data.frame(var_name,t(summary(glm(formula, family = "binomial",data = dfcomplete))$coef[2,c(1,2,4)]))
    Model1menopause <- bind_rows(Model1menopause,coef)
}
Model1menopause$LWR <- Model1menopause$Estimate - 1.96*Model1menopause$Std..Error
Model1menopause$UPR <- Model1menopause$Estimate + 1.96*Model1menopause$Std..Error

Model2menopause <- data.frame()
for (i in 1:nrow(regmenopause)) {
    var_name <- regmenopause$FeatName[i]
    formula <- as.formula(paste("MenopauseYN ~ Age + HT + DM + Smoking + ", var_name))
    coef <- data.frame(var_name,t(summary(glm(formula, family = "binomial",data = dfcomplete))$coef[5,c(1,2,4)]))
    Model2menopause <- bind_rows(Model2menopause,coef)
}
Model2menopause$LWR <- Model2menopause$Estimate - 1.96*Model2menopause$Std..Error
Model2menopause$UPR <- Model2menopause$Estimate + 1.96*Model2menopause$Std..Error

