#### Regression analyses Sex

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
pred_sex <- rio::import("data/feature_importance_sex.txt")
rownames(pred_sex) <- pred_sex$FeatName

## Determine transformation
hist(dfcomplete$Zotu7202, breaks = 60)
hist(log(dfcomplete$Zotu7202 + 1), breaks = 60)
dfcomplete$CLRZotu7207 <- clr(dfcomplete$Zotu7207 + 1)
hist(dfcomplete$CLRZotu7207, breaks = 60)



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
Model1sex$Model <- "Model 1"
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
Model2sex$Model <- "Model 2"
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
Model3sex$Model <- "Model 3"
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
Model4sex$Model <- "Model 4"
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
Model5sex$Model <- "Model 5"
write.table(Model5sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model5sex <- Model5sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model5sex, "results/Model5sex.csv")

## Combine datasets
Modelsex <- rbind(Model1sex, Model2sex, Model3sex, Model4sex, Model5sex)
Modelsex[, c(2,3,4,5)] <- lapply(Modelsex[, c(2,3,4,5)], as.numeric)


## Change zotu names to taxonomy
Modelsex$ASV <- tax$Tax[match(Modelsex$ASV, tax$ASV)]


## Make ggplot
forest_plot_sex <- ggplot(Modelsex, aes(x = Estimate, y = ASV, color = Model)) +
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

forest_plot_sex



