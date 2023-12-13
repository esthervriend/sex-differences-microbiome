## Diversity indices

## Libraries
library(phyloseq)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(vegan)

## Load data
p <- readRDS("data/phyloseq_sampledata.RDS")
df <- readRDS("")
tax <- readRDS("data/tax_table.RDS")

tab <- as(p@otu_table, 'matrix')
counts <- sample_sums(p@otu_table)
tab <- as.data.frame(tab/sample_sums(p)*100)
rowSums(tab) # samples should all sum up to 100%

## Diversity metrics
# Shannon
shannon <- vegan::diversity(tab, index = 'shannon')
df_shan <- data.frame(ID = as.integer(names(shannon)), shannon = shannon)
shanpl <- ggplot(df_shan, aes(x = shannon)) +
    geom_histogram(color = "black", 
                   fill = "royalblue", alpha = 0.8) + 
    theme_Publication() +
    xlab("Shannon index") +
    ggtitle("Shannon index")
ggsave("results/shannon_index_hist.pdf")

df_shan <- left_join(df_shan, ma, by = "sampleID")
ggplot(data = df_shan, aes(x = Sex, y = shannon, fill = Sex)) +
    geom_violin() +
    geom_boxplot(width = 0.1, fill = "white", outlier.shape = NA) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Alpha diversity (Shannon)", y = "Shannon index", x="") +
    stat_compare_means(method = "wilcox.test")
ggsave("results/shannon.pdf", device = "pdf", width = 6, height = 5)

## Species richness
specrich <- specnumber(tab)
dfspec <- data.frame(ID = as.integer(names(richness)), richness = richness)
specpl <- ggplot(dfspec, aes(x = specrich)) +
    geom_histogram(color = "black", 
                   fill = "darkgreen", alpha = 0.8) + 
    theme_Publication() +
    xlab("Number of species") +
    ggtitle("Species richness")
ggsave("results/species_richness_hist.pdf")

dfspec <- left_join(dfspec, ma, by = "ID")
ggplot(data = dfveg, aes(x = Sex, y = richness, fill = Sex)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_lancet(guide = FALSE) + 
    labs(title = "Species richness", y = "Number of species", x = "") +
    stat_compare_means(method = "wilcox.test")
ggsave("results/richness.pdf", device = "pdf", width = 6, height = 5)

## Faith's PD
#faith <- picante::pd(tab@otu_table, tree = tab@phy_tree)
dffai <- as.data.frame(faith)
dffai$ID <- as.integer(rownames(faith))
dffai <- left_join(dffai, ma, by = "ID")
ggplot(data = dffai, aes(x = Sex, y = PD, fill = Sex)) +
    geom_violin()+
    geom_boxplot(outlier.shape = NA, fill = "white", width = 0.1) +
    theme_Publication() + 
    scale_fill_lancet(guide = "none") + 
    labs(title = "Alpha diversity (Faith's PD)", y = "Faith's phylogenetic diversity") +
    stat_compare_means(method = "wilcox.test")
ggsave("results/faiths.pdf", device = "pdf", width = 4, height = 5)

#### Ordination plots ####
## Bray-Curtis PCoA
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)
dbray <- pcoord$vectors[, c(x_comp, y_comp)]
dbray <- as.data.frame(dbray)

# add metadata / covariates
dbray$sampleID <- as.integer(rownames(dbray))
dbray <- left_join(dbray, ma, by = 'sampleID')
head(dbray)

pl2 <- dbray %>% 
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
pl2
ggsave("results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)

## PERMANOVA / adonis
set.seed(1234)
# distance matrix and metadata (df with the outcome / covariates) must have the same sample order as bray distance matrix / distance object
all(ma$sampleID == sample_names(tc)) # FALSE
mb <- ma %>%
    slice(match(sample_names(tc), sampleID))
all(mb$sampleID == sample_names(tc)) # TRUE
dim(mb)
names(mb)
res <- adonis(bray ~ group, data = mb) # PERMANOVA

## Bray curtis plot with PERMANOVA annotation
pl <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    ggtitle("Bray curtis distance PCoA") +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
pl <- pl + annotate("text", x = 0.2, y = 0.4, label = paste0("PERMANOVA p = ", res$aov.tab[3,6]))

ggsave("results/PCA_BrayCurtis_permanova.pdf", device = "pdf", width = 6, height = 5)

## From another script: 
## CLR-transformed PCA
pseudocount <- min(tab[tab != 0]) / 2 # could also just be 1
pc <- mixOmics::pca(mat + pseudocount, center = T, scale = F, logratio = 'CLR') # scale should be F when CLR is used
expl_variance <- pc$explained_variance * 100
x_comp <- paste0('PC',1)
y_comp <- paste0('PC',2)
dfclr <- pc$x[, c(x_comp, y_comp)]
dfclr <- as.data.frame(dfclr)
pc$cum.var
# add metadata
rownames(dfclr)
dfclr$Sex <- ma$Sex[match(rownames(dfclr), ma$sampleID)]

# plot PCA CLR-transformed
pl <- dfclr %>% 
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(color = Sex), size = 2) +
    xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PComp, not PCoord
    ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
    #scale_shape_manual(values = c(21, 22, 23, 24)) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed') + 
    stat_ellipse(aes(color = group), type = "norm")
pl
ggsave("results/PCA_CLR.pdf", device = "pdf", width = 6, height = 5)

## Weighted UniFrac
wunifrac <- UniFrac(tab, normalized = T, weighted = T)
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)

# get PCoA coordinates
dfpc <- pcoord$vectors[, c(x_comp, y_comp)]
dfpc <- as.data.frame(dfpc)

# add metadata / covariates
dfpc$ID <- as.integer(rownames(dfpc))
dfpc <- left_join(dfpc, ma, by = 'ID')

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = group), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_lancet() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = group), type = "norm")
pl
ggsave("results/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 6, height = 5)
