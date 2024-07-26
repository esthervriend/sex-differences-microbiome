## Process results machine learning

library(tidyverse)
library(ggsci)
library(stringr)

options(scipen=999)

source("scripts/functions.R")

## Plot assembled figure composition sex
path_true <- 'male_female/output_XGB_class_male_female_2024_01_26__23-56-15'
path_permuted <- 'male_female/output_XGB_class_male_female_2023_12_22__20-25-22_PERMUTED'
data_path <- 'male_female/input_data'
labels <- c("Men", "Women")

compared_to_permuted_class(path_true, path_permuted)

pl2 <- plot_feature_importance_class(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/16s_comp_sex.pdf",
       width = 14, height = 18)
ggsave(plarr1, filename = "results/ml_figures/16s_comp_sex.png",
       width = 14, height = 18)

## Plot assembled figure composition menopause
path_true <- 'menopause/output_XGB_class_menopause_2024_01_24__10-09-15'
path_permuted <- 'menopause/output_XGB_class_menopause_2024_01_26__23-56-32_PERMUTED'
data_path <- 'menopause/input_data'
labels <- c("Postmenopausal", "Premenopausal")

compared_to_permuted_class(path_true, path_permuted)

pl2 <- plot_feature_importance_class(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_tests_top(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/16s_comp_meno.pdf", width = 14, height = 18)
ggsave(plarr1, filename = "results/ml_figures/16s_comp_meno.png", width = 14, height = 18)
