## Process results machine learning

library(tidyverse)
library(ggsci)
library(stringr)
library(ggpubr)

source("scripts/functions_shotgun.R")

#### Composition ####
## Plot assembled figure composition sex
path_true <- 'sex_metagen/output_XGB_class_sex_metagen_2024_07_25__21-26-07'
data_path <- 'sex_metagen/input_data'
labels <- c("Men", "Women")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/comp_sex.pdf",
       width = 14, height = 18)
ggsave(plarr1, filename = "results/ml_figures/comp_sex.png",
       width = 14, height = 18)


## Plot assembled figure composition menopause
path_true <- 'menopause_metagen/output_XGB_class_menopause_metagen_2024_07_25__22-30-05'
data_path <- 'menopause_metagen/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/comp_menopause.pdf",
       width = 14, height = 15)
ggsave(plarr1, filename = "results/ml_figures/comp_menopause.png",
       width = 14, height = 15)


#### Pathways ####
## Plot assembled figure pathways sex
path_true <- 'sex_pathways/output_XGB_class_sex_pathways_2024_03_03__20-18-59'
data_path <- 'sex_pathways/input_data'
labels <- c("Men", "Women")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr1 <- ggarrange(ggarrange(svg_grob, NULL, nrow = 2, heights = c(2.0,0.3)), pl2, 
                    nrow = 1, widths = c(0.8, 1.3), labels = c("A", "B"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_sex.pdf",
       width = 15, height = 20)
ggsave(plarr1, filename = "results/ml_figures/pathways_sex.png",
       width = 15, height = 20)

## Plot assembled figure pathways menopause
path_true <- 'menopause_pathways/output_XGB_class_menopause_pathways_2024_03_03__20-28-16'
data_path <- 'menopause_pathways/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_pathways(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/ml_figures/pathways_menopause.pdf",
       width = 14, height = 15)
ggsave(plarr1, filename = "results/ml_figures/pathways_menopause.png",
       width = 14, height = 15)
