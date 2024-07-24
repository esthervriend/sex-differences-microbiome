#### Regression analyses Menopause

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)
library(dplyr)
library(ggsci)
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
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
tab <- as.data.frame(t(as(phyloseq@otu_table, 'matrix')))
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/taxtable_rarefied.RDS")

# Best predictors
pred_meno <- rio::import("data/feature_importance_menopause.txt")
rownames(pred_meno) <- pred_meno$FeatName


## Regression analyses Menopause
# Model 1 (Unadjusted)
pred_meno <- pred_meno %>% slice_head(n = 20)

Model1menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model1menopause <- bind_rows(Model1menopause,coef)
}
Model1menopause$Pr...t.. <- p.adjust(Model1menopause$Pr...t.., method = "fdr", n = length(Model1menopause$Pr...t..))
Model1menopause$LWR <- Model1menopause$Estimate - 1.96*Model1menopause$Std..Error
Model1menopause$UPR <- Model1menopause$Estimate + 1.96*Model1menopause$Std..Error
Model1menopause <- Model1menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model1menopause$Model <- "Unadjusted"
Model1menopause$ASV <- make.unique(tax$Tax[match(Model1menopause$ASV, tax$ASV)])
write.table(Model1menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model1menopause <- Model1menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model1menopause, "results/Model1menopause.csv")

# Model 2 (+ Age & CVD)
Model2menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + HT + DM + CurrSmoking"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model2menopause <- bind_rows(Model2menopause,coef)
}
Model2menopause$Pr...t.. <- p.adjust(Model2menopause$Pr...t.., method = "fdr", n = length(Model2menopause$Pr...t..))
Model2menopause$LWR <- Model2menopause$Estimate - 1.96*Model2menopause$Std..Error
Model2menopause$UPR <- Model2menopause$Estimate + 1.96*Model2menopause$Std..Error
Model2menopause <- Model2menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model2menopause$Model <- "+Age, BMI, HT, DM, smoking"
Model2menopause$ASV <- make.unique(tax$Tax[match(Model2menopause$ASV, tax$ASV)])
write.table(Model2menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model2menopause <- Model2menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model2menopause, "results/Model2menopause.csv")


# Model 3 (+ Age & CVD & diet)
Model3menopause <- data.frame()
for (i in 1:nrow(pred_meno)) {
  var_name <- pred_meno$FeatName[i]
  formula <- as.formula(paste("log(", var_name, "+1) ~ MenopauseYn + Age + BMI + HT + DM + CurrSmoking + Alcohol + TotalCalories + Fibre"))
  coef <- data.frame(var_name,t(summary(lm(formula, data = dfcomplete))$coef[2,c(1,2,4)]))
  Model3menopause <- bind_rows(Model3menopause,coef)
}
Model3menopause$Pr...t.. <- p.adjust(Model3menopause$Pr...t.., method = "fdr", n = length(Model3menopause$Pr...t..))
Model3menopause$LWR <- Model3menopause$Estimate - 1.96*Model3menopause$Std..Error
Model3menopause$UPR <- Model3menopause$Estimate + 1.96*Model3menopause$Std..Error
Model3menopause <- Model3menopause %>%
  select(ASV=var_name, Estimate, LWR, UPR, Pvalue=Pr...t..)
Model3menopause$Model <- "+Diet"
Model3menopause$ASV <- make.unique(tax$Tax[match(Model3menopause$ASV, tax$ASV)])
write.table(Model3menopause, "clipboard", sep="\t", dec=",", col.names=NA)
Model3menopause <- Model3menopause %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Model3menopause, "results/Model3menopause.csv")

## Combine datasets
Modelmenopause <- rbind(Model1menopause, Model2menopause, Model3menopause)
Modelmenopause[, c(2,3,4,5)] <- lapply(Modelmenopause[, c(2,3,4,5)], as.numeric)

## Make ggplot
forest_plot_menopause <- ggplot(Modelmenopause, aes(x = Estimate, y = fct_rev(fct_inorder(ASV)), color = fct_rev(fct_inorder(Model)), shape = Pvalue > 0.05)) +
  geom_point(position = position_dodge(width = 0.6), size = 2.0) +
  geom_vline(aes(xintercept = 0), size = .50, linetype = "dashed") + 
  geom_errorbarh(aes(xmin = LWR, xmax = UPR), height = 0.5, position = position_dodge(width = 0.6)) +
  theme_Publication() +
  theme(axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = 'right') +
  scale_color_nejm(breaks=c('Unadjusted', '+Age, BMI, HT, DM, smoking', '+Diet')) +
  labs(x = "Estimate and 95% CI for postmenopausal women") +
  scale_shape_manual(values = c(16, 1)) +
  guides(color = guide_legend(title = NULL), shape = "none") +
  scale_x_continuous(breaks = seq(-2,2, by = 0.5)) 
forest_plot_menopause

forest_plot_menopause <- forest_plot_menopause + ggtitle("Best predicting microbes for menopause") + theme(plot.title = element_text(size = 15))
forest_plot_menopause

ggsave("results/forest_plot_menopause.tiff", plot=forest_plot_menopause, units="in", width=10, height=7, dpi=600, compression = 'lzw')

## Combine with plot for sex
forest_plot_sexmenopause <- ggarrange(
  forest_plot_sex,
  forest_plot_menopause,
  labels = c("A", "B"),
  common.legend = TRUE, 
  legend = "right",
  ncol = 1,
   align = "v"
)
forest_plot_sexmenopause
ggsave("results/forest_plot_sexmenopause.tiff", plot=forest_plot_sexmenopause, units="in", width=12, height=12, dpi=600, compression = 'lzw')

