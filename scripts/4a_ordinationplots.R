## Calculate distances and save

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)

theme_Publication <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    library(stringr)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(0.8)),
                axis.title.y = element_text(angle=90, vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(size = rel(0.7)),
                axis.text.x = element_text(angle = 0), 
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
                strip.text = element_text(face="bold"),
                plot.caption = element_text(size = rel(0.5), face = "italic")
        ))
    
} 

## Load data
phydata <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio:: import("data/clinicaldata.RDS")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phyloseq@otu_table, 'matrix'))

# Bray-Curtis PCoA
bray <- vegan::vegdist(tab, method = 'bray')
saveRDS(bray, 'data/braycurtis.RDS')
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

braycurt <- dbray %>% 
        ggplot(aes(Axis.1, Axis.2)) +
        geom_point(aes(color = Sex), size = 2) +
        ggtitle("PCoA Bray-Curtis distance") +
        xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
        ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
        scale_fill_manual(values = rev(pal_nejm()(2))) +
        theme_Publication() +
        guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
        stat_ellipse(aes(color = Sex), type = "norm")

ggsave(braycurt, "results/PCA_BrayCurtis.pdf", device = "pdf", width = 6, height = 5)
ggsave(braycurt, "results/PCA_BrayCurtis.svg", device = "svg", width = 6, height = 5)

## PERMANOVA / adonis
set.seed(1234)
# distance matrix and metadata (df with the outcome / covariates) must have the same sample order as bray distance matrix / distance object
all(df_new$sampleID == sample_names(tab)) # FALSE
dfanova <- df_new %>%
    slice(match(sample_names(tab), ID))
all(dfanova$sampleID == sample_names(tab)) # TRUE
dim(dfanova)
names(dfanova)
res <- adonis(bray ~ Sex, data = dfanova) # PERMANOVA
print(res)

## Bray curtis plot with PERMANOVA annotation
planova <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    scale_fill_nejm() +
    geom_point(aes(color = Sex), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    ggtitle("Bray curtis distance PCoA") +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = Sex), type = "norm")
planova <- planova + annotate("text", x = 0.2, y = 0.4, label = paste0("PERMANOVA p = ", res$aov.tab[3,6]))

ggsave(planova, "results/PCA_BrayCurtis_permanova.pdf", device = "pdf", width = 6, height = 5)

## Weighted UniFrac
wunifrac <- UniFrac(phydata, normalized = T, weighted = T)
saveRDS(wunifrac, "data/weighted_unifrac.RDS")
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
x_comp <- paste0('Axis.',1)
y_comp <- paste0('Axis.',2)

# get PCoA coordinates
dfpc <- pcoord$vectors[, c(x_comp, y_comp)]
dfpc <- as.data.frame(dfpc)

# add metadata / covariates
dfpc$ID <- as.integer(rownames(dfpc))
dfpc <- left_join(dfpc, df_new, by = 'ID')

pl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 2) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(pl, "results/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 6, height = 5)

## CLR-transformed PCA
pseudocount <- 1
pc <- mixOmics::pca(tab_matrix + pseudocount, center = T, scale = F, logratio = 'CLR') # scale should be F when CLR is used
expl_variance <- pc$explained_variance * 100
x_comp <- paste0('PC',1)
y_comp <- paste0('PC',2)
dfclr <- pc$x[, c(x_comp, y_comp)]
dfclr <- as.data.frame(dfclr)
print(pc$cum.var)
# add metadata
dfclr$Sex <- df_new$Sex[match(rownames(dfclr), df_new$ID)]

# Plot PCA CLR-transformed
plclr <- dfclr %>% 
    ggplot(aes(PC1, PC2)) +
    geom_point(aes(color = Sex), size = 2) +
    xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PComp, not PCoord
    ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    guides(fill = guide_legend(override.aes = list(shape = 21))) +
    ggtitle('PCA CLR-transformed') + 
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(plclr, "results/PCA_CLR.pdf", device = "pdf", width = 6, height = 5)
