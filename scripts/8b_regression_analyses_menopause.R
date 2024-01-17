#### Regression analyses Menopause

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


## ggplot theme
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  library(stringr)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "bottom",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing  = unit(0, "cm"),
            # legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))}

## Load data
# Clinical data and phyloseq
df <- readRDS("data/clinicaldata.RDS")
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
tab <- as.data.frame(t(as(phyloseq@otu_table, 'matrix')))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/tax_table.RDS")

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
Model1Menopause$Model <- "Model 1"
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
Model2Menopause$Model <- "Model 2"
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
Model3Menopause$Model <- "Model 3"
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
Model4Menopause$Model <- "Model 4"
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
Model5Menopause$Model <- "Model 5"
write.table(Model5Menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model5Menopause <- Model5Menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model5Menopause, "results/Model5Menopause.csv")


## Combine datasets
ModelMenopause <- rbind(Model1Menopause, Model2Menopause, Model3Menopause, Model4Menopause, Model5Menopause)
ModelMenopause[, c(2,3,4,5)] <- lapply(ModelMenopause[, c(2,3,4,5)], as.numeric)


## Change zotu names to taxonomy
ModelMenopause$ASV <- tax$Tax[match(ModelMenopause$ASV, tax$ASV)]


## Make ggplot
forest_plot_menopause <- ggplot(ModelMenopause, aes(x = Estimate, y = ASV, color = Model)) +
  geom_point() +
  geom_vline(aes(xintercept = 1), size = .50, linetype = "dashed") + 
  geom_errorbarh(aes(xmin = LWR, xmax = UPR), height = 0.5) +
  facet_wrap(~ Model, ncol = 5) +
  theme_Publication() +
  theme(axis.ticks.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'none') +
  scale_color_nejm() +
  scale_x_continuous(breaks = c(0,1,2), labels = c(0,1,2))

forest_plot_menopause



