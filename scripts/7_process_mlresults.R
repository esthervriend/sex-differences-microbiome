## Process results machine learning

library(dplyr)
library(ggplot2)
library(ggsci)
library(stringr)

source("scripts/functions.R")

path_true <- 'male_female/output_XGB_class_male_female_2023_12_20__22-56-43'
path_permuted <- 'male_female/output_XGB_class_male_female_2023_12_22__20-25-22_PERMUTED'
data_path <- 'male_female/input_data'
labels <- c("Male", "Female")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)

path_true <- 'menopause/output_XGB_class_menopause_2023_12_20__22-57-00'
path_permuted <- 'menopause/output_XGB_class_menopause_2023_12_21__20-51-36_PERMUTED'
data_path <- 'menopause/input_data'
labels <- c("Postmenopausal", "Premenopausal")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
