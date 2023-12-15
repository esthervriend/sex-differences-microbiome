#### Diversity indices

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)
library(rio)

## Load data
phyloseq <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio:: import("data/clinicaldata.RDS")
tax <- readRDS("data/tax_table.RDS")

tab <- as(phyloseq@otu_table, 'matrix')
counts <- sample_sums(phyloseq@otu_table)
tab <- as.data.frame(tab/sample_sums(phyloseq)*100)
rowSums(tab) # samples should all sum up to 100%

## Diversity metrics
# Shannon plots
shannon <- vegan::diversity(tab, index = 'shannon')
df_shan <- data.frame(ID = as.integer(names(shannon)), shannon = shannon)
shanpl <- ggplot(df_shan, aes(x = shannon)) +
    geom_histogram(color = "black", 
                   fill = "royalblue", alpha = 0.8) + 
    theme_classic() + 
    theme(axis.text = element_text(color = "black", size = 16),
          axis.title = element_text(color = "black", size = 16)) + 
    xlab("Shannon index") 
shanpl
ggsave("results/shannon_index_hist.tiff")

df_shan <- left_join(df_shan, df_new, by = "ID")
shanviolin <- ggplot(data = df_shan, aes(x = Sex, y = shannon, fill = Sex)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_classic() + 
    labs( y = "Shannon index", x="") + 
    theme(axis.text = element_text(color = "black", size = 16),
      axis.title = element_text(color = "black", size = 16),
      legend.text = element_text(color = "black", size = 16),
      legend.title = element_text(color = "black", size = 16)) +
    stat_compare_means(method = "wilcox.test", size = 5, label.x.npc = 0.31)
shanviolin  
ggsave("results/shannon.tiff", width = 6, height = 5)


## Ordination plots
# Bray-Curtis PCoA
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)
dbray <- pcoord$vectors[, c(x_comp, y_comp)]
dbray <- as.data.frame(dbray)

# add metadata / covariates
dbray$ID <- as.integer(rownames(dbray))
dbray <- left_join(dbray, df_new, by = 'ID')
head(dbray)

Braycurtis <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
Braycurtis
ggsave("results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)