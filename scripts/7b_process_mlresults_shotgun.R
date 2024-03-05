## Process results machine learning

library(tidyverse)
library(ggplot2)
library(ggsci)
library(stringr)

source("scripts/functions_shotgun.R")

path_true <- 'sex_metagen/output_XGB_class_sex_metagen_2024_03_03__01-40-28'
data_path <- 'sex_metagen/input_data'
labels <- c("Male", "Female")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'menopause_metagen/output_XGB_class_menopause_metagen_2024_03_03__01-59-48'
data_path <- 'menopause_metagen/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'sex_pathways/output_XGB_class_sex_pathways_2024_03_03__20-18-59'
data_path <- 'sex_pathways/input_data'
labels <- c("Male", "Female")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'menopause_pathways/output_XGB_class_menopause_pathways_2024_03_03__20-28-16'
data_path <- 'menopause_pathways/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)
