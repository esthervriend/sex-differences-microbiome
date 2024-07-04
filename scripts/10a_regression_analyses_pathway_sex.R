#### Regression analyses pathway Sex

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)
library(dplyr)
library(ggsci)
library(forcats)
library(cowplot)


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
tab <- readRDS("data/pathways_filtered.RDS")
tab <- as.data.frame(t(tab))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/pathway_keys.RDS")

dfcompletediet <- dfcomplete[!is.na(dfcomplete$Sodium),]

# Best predictors
pred_sex_pathway <- rio::import("data/sex_pathways_feature_importance.txt")
rownames(pred_sex_pathway) <- pred_sex_pathway$FeatName

## Regression analyses sex
# Model 1 (Crude)
pred_sex_pathway <- pred_sex_pathway %>% slice_head(n = 20)

Model1sex <- data.frame()
for (i in 1:nrow(pred_sex_pathway)) {
  var_name <- pred_sex_pathway$FeatName[i] 
  formula <- paste0("log(`", var_name, "` + 1) ~ Sex")
  coef <- summary(lm(formula, data = dfcompletediet))$coef[2, c(1, 2, 4)]
  coef <- data.frame(var_name, t(coef))
  Model1sex <- bind_rows(Model1sex,coef)
}
Model1sex$Pr...t.. <- p.adjust(Model1sex$Pr...t.., method = "fdr", n = length(Model1sex$Pr...t..))
Model1sex$LWR <- Model1sex$Estimate - 1.96*Model1sex$Std..Error
Model1sex$UPR <- Model1sex$Estimate + 1.96*Model1sex$Std..Error
Model1sex <- Model1sex %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model1sex$Model <- "Unadjusted"
write.table(Model1sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model1sex <- Model1sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model1sex, "results/Model1sexpathway.csv")

# Model 2 (+ Age & CVD)
Model2sex <- data.frame()
for (i in 1:nrow(pred_sex_pathway)) {
  var_name <- pred_sex_pathway$FeatName[i]
  formula <- paste0("log(`", var_name, "` + 1) ~ Sex + Age + BMI + HT + DM + CurrSmoking")
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcompletediet))$coef[2,c(1,2,4)]))
  Model2sex <- bind_rows(Model2sex,coef)
}
Model2sex$Pr...t.. <- p.adjust(Model2sex$Pr...t.., method = "fdr", n = length(Model2sex$Pr...t..))
Model2sex$LWR <- Model2sex$Estimate - 1.96*Model2sex$Std..Error
Model2sex$UPR <- Model2sex$Estimate + 1.96*Model2sex$Std..Error
Model2sex <- Model2sex %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model2sex$Model <- "+Age, BMI, HT, DM, smoking"
write.table(Model2sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model2sex <- Model2sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model2sex, "results/Model2sexpathway.csv")


# Model 3 (+ Age & CVD & diet)
Model3sex <- data.frame()
for (i in 1:nrow(pred_sex_pathway)) {
  var_name <- pred_sex_pathway$FeatName[i]
  formula <- paste0("log(`", var_name, "` + 1) ~ Sex + Age + BMI + HT + DM + CurrSmoking + Alcohol + TotalCalories + Fibre")
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcompletediet))$coef[2,c(1,2,4)]))
  Model3sex <- bind_rows(Model3sex,coef)
}
Model3sex$Pr...t.. <- p.adjust(Model3sex$Pr...t.., method = "fdr", n = length(Model3sex$Pr...t..))
Model3sex$LWR <- Model3sex$Estimate - 1.96*Model3sex$Std..Error
Model3sex$UPR <- Model3sex$Estimate + 1.96*Model3sex$Std..Error
Model3sex <- Model3sex %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model3sex$Model <- "+Diet"
write.table(Model3sex, "clipboard", sep="\t", dec=",", col.names=NA)
Model3sex <- Model3sex %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model3sex, "results/Model3sexpathway.csv")

## Combine datasets
Modelsex <- rbind(Model1sex, Model2sex, Model3sex)
Modelsex[, c(2,3,4,5)] <- lapply(Modelsex[, c(2,3,4,5)], as.numeric)

## Change ASV names to taxonomy
Modelsex$ASV <- tax$label[match(Modelsex$ASV, tax$key)]


## Make ggplot
forest_plot_sex <- ggplot(Modelsex, aes(x = Estimate, y = fct_rev(fct_inorder(ASV)), color = fct_rev(fct_inorder(Model)), shape = Pvalue > 0.05)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.0) +
  geom_vline(aes(xintercept = 0), size = .50, linetype = "dashed") + 
  geom_errorbarh(aes(xmin = LWR, xmax = UPR), height = 0.5, position = position_dodge(width = 0.6)) +
  theme_Publication() +
  theme(axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right') +
  scale_color_nejm(breaks=c('Unadjusted', '+Age, BMI, HT, DM, smoking', '+Diet')) +
  labs(x = "Log-transformed estimate and 95% CI for women") + 
  scale_shape_manual(values = c(16, 1)) +
  guides(color = guide_legend(title = NULL), shape = "none") +
  scale_x_continuous(breaks = seq(-2,2, by = 0.5))
forest_plot_sex

forest_plot_sex <- forest_plot_sex + ggtitle("Best predicting pathways for sex") + theme(plot.title = element_text(size = 15))

forest_plot_sex <- plot_grid(
  forest_plot_sex,
  label_size = 15)
forest_plot_sex
ggsave("results/forest_plot_sex_pathways.tiff", plot=forest_plot_sex, units="in", width=13, height=7, dpi=600, compression = 'lzw')




