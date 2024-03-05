## Calculate distances, plot PCoA and PCA

## Libraries
library(phyloseq)
library(vegan)
library(tidyverse)
library(ggplot2)
library(ggpubr)
library(ggsci)
library(doParallel)
registerDoParallel(24)

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

#### Load data ####
phydata <- readRDS("data/phyloseq_sampledata.RDS")
df_new <- rio::import("data/clinicaldata.RDS")
df_meno <- df_new %>% filter(Sex == "Female")
tab <- as.data.frame(t(as(phydata@otu_table, 'matrix')))
tab_matrix <- t(as(phydata@otu_table, 'matrix'))
phydata_meno <- prune_samples(df_meno$ID, phydata)
tab_meno <- as.data.frame(t(as(phydata_meno@otu_table, 'matrix')))

#### Output folder ####
resultsfolder <- "results/ordination"
dir.create(resultsfolder, showWarnings = FALSE)

#### Bray-Curtis distance: male-female ####
print('Bray-Curtis distance total dataset')
bray <- vegan::vegdist(tab, method = 'bray')
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance <- pcoord$values$Rel_corr_eig * 100
head(expl_variance)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, df_new, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
set.seed(1234)
dfanova <- df_new %>%
    slice(match(sample_names(phydata), ID)) # distance matrix and metadata must have the same sample order
all(dfanova$ID == sample_names(phydata)) # TRUE
dim(dfanova)
res1 <- adonis2(bray ~ Sex, data = dfanova) # PERMANOVA
print(res1)

print('plotting..')
braycurt <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    theme_Publication() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm") + 
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
        label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", format(round(res1$R2[1],3), nsmall = 3)))
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis.pdf", device = "pdf", width = 8, height = 8)

#### Bray-Curtis menopausal status ####
print('Bray-Curtis distance female subjects')
bray_meno <- vegan::vegdist(tab_meno, method = 'bray')
pcoord <- ape::pcoa(bray_meno, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
dbray_meno <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray_meno <- as.data.frame(dbray_meno)
dbray_meno$ID <- rownames(dbray_meno)
dbray_meno <- left_join(dbray_meno, df_meno, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
set.seed(1234)
dfanova_meno <- df_meno %>%
    slice(match(sample_names(phydata_meno), ID)) # distance matrix and metadata must have the same sample order
all(dfanova_meno$ID == sample_names(phydata_meno)) # TRUE
res2 <- adonis2(bray_meno ~ MenopauseYn, data = dfanova_meno) # PERMANOVA
print(res2)

print('plotting..')
braycurt_meno <- dbray_meno %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_nejm()(4)[3:4]) +
    theme_Publication() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = MenopauseYn), type = "norm") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
        label = str_c("PERMANOVA: p = ", res2$`Pr(>F)`, ", r2 = ", format(round(res2$R2[1],3), nsmall = 3)))
ggsave(braycurt_meno, filename = "results/ordination/PCoA_BrayCurtis_meno.pdf", device = "pdf", width = 8, height = 8)

#### Weighted UniFrac ####
print('Weighted UniFrac')
wunifrac <- UniFrac(phydata, normalized = T, weighted = T, parallel = T)
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
head(expl_variance)
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, df_new, by = 'ID') # add metadata / covariates

print('plotting..')
unifracpl <- dfpc %>% 
                ggplot(aes(Axis.1, Axis.2)) +
                geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
                xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
                ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
                theme_Publication() +
                scale_color_manual(values = rev(pal_nejm()(2))) +
                ggtitle('PCoA Weighted UniFrac') +
                guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
                labs(color = "") +
                stat_ellipse(aes(color = Sex), type = "norm")
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", device = "pdf", width = 8, height = 8)

print('plotting..')
unifrac_meno <- dfpc %>% filter(Sex == "Female") %>% 
                    ggplot(aes(Axis.1, Axis.2)) +
                    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
                    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
                    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
                    theme_Publication() +
                    scale_color_manual(values = pal_nejm()(4)[3:4]) +
                    labs(color = "") +
                    ggtitle('PCoA Weighted UniFrac') +
                    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
                    stat_ellipse(aes(color = MenopauseYn), type = "norm")
ggsave(unifrac_meno, filename = "results/ordination/PCoA_WeightedUnifrac_meno.pdf", device = "pdf", width = 8, height = 8)

#### CLR-transformed PCA ####
print('CLR-transformation and PCA')
pseudocount <- 1
pc <- mixOmics::pca(tab_matrix + pseudocount, center = T, scale = F, logratio = 'CLR') # scale should be F when CLR is used
expl_variance <- pc$explained_variance * 100
dfclr <- pc$x[, c('PC1', 'PC2')]
dfclr <- as.data.frame(dfclr)
print(pc$cum.var)
dfclr$Sex <- df_new$Sex[match(rownames(dfclr), df_new$ID)] 
dfclr$MenopauseYn <- df_new$MenopauseYn[match(rownames(dfclr), df_new$ID)]

print('plotting..')
plclr <- dfclr %>% 
            ggplot(aes(PC1, PC2)) +
            geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
            xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PComp, not PCoord
            ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
            theme_Publication() +
            scale_color_manual(values = rev(pal_nejm()(2))) +
            guides(fill = guide_legend(override.aes = list(shape = 21))) +
            ggtitle('PCA CLR-transformed') + 
            labs(color = "") +
            stat_ellipse(aes(color = Sex), type = "norm")
ggsave(plclr, filename = "results/ordination/PCA_CLR.pdf", device = "pdf", width = 8, height = 8)

plclr_meno <- dfclr %>% 
                filter(Sex == "Female") %>% 
                    ggplot(aes(PC1, PC2)) +
                    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
                    xlab(paste0('PC1 (', round(expl_variance[1], digits = 1),'%)')) + # PComp, not PCoord
                    ylab(paste0('PC2 (', round(expl_variance[2], digits = 1),'%)')) +
                    theme_Publication() +
                    scale_color_manual(values = pal_nejm()(4)[3:4]) +
                    guides(fill = guide_legend(override.aes = list(shape = 21))) +
                    ggtitle('PCA CLR-transformed') + 
                    labs(color = "") +
                    stat_ellipse(aes(color = MenopauseYn), type = "norm")
ggsave(plclr_meno, filename = "results/ordination/PCA_CLR_meno.pdf", device = "pdf", width = 8, height = 8)


#### Shotgun beta diversity ####
dfshot <- rio::import("data/clinicaldata_shotgun.RDS")
wunifrac <- rio::import("data/metaphlan/diversity/combined_table_weighted-unifrac.tsv") 
rownames(wunifrac) <- wunifrac$V1
wunifrac$V1 <- NULL
all(rownames(wunifrac) == colnames(wunifrac))
IDs <- str_remove(rownames(wunifrac), "_T1")
colnames(wunifrac) <- rownames(wunifrac) <- IDs

#### Bray-Curtis distance: male-female ####
print('Bray-Curtis distance total dataset')
bray <- rio::import("data/metaphlan/diversity/combined_table_bray-curtis.tsv")
rownames(bray) <- bray$V1
bray$V1 <- NULL
all(rownames(bray) == colnames(bray))
IDs <- str_remove(rownames(bray), "_T1")
colnames(bray) <- rownames(bray) <- IDs
pcoord <- ape::pcoa(bray, correction = "cailliez")
str(pcoord$values)
expl_variance <- pcoord$values$Rel_corr_eig * 100
head(expl_variance)
dbray <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray <- as.data.frame(dbray)
dbray$ID <- rownames(dbray)
dbray <- left_join(dbray, dfshot, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
set.seed(1234)
dfanova <- dfshot %>%
    slice(match(rownames(bray), ID)) # distance matrix and metadata must have the same sample order
all(dfanova$ID == rownames(bray)) # TRUE
dim(dfanova)
res1 <- adonis2(bray ~ Sex, data = dfanova) # PERMANOVA
print(res1)

print('plotting..')
braycurt <- dbray %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    theme_Publication() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm") + 
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res1$`Pr(>F)`, ", r2 = ", 
                           format(round(res1$R2[1],3), nsmall = 3)))
ggsave(braycurt, filename = "results/ordination/PCoA_BrayCurtis_shotgun.pdf", 
       device = "pdf", width = 5, height = 5)

#### Bray-Curtis menopausal status ####
print('Bray-Curtis distance female subjects')
menopause <- dfshot %>% filter(Sex == "Female")
IDs_meno <- menopause$ID
bray_meno <- bray[which(rownames(bray) %in% IDs_meno), which(colnames(bray) %in% IDs_meno)]
pcoord <- ape::pcoa(bray_meno, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
dbray_meno <- pcoord$vectors[, c('Axis.1', 'Axis.2')]
dbray_meno <- as.data.frame(dbray_meno)
dbray_meno$ID <- rownames(dbray_meno)
dbray_meno <- left_join(dbray_meno, menopause, by = 'ID') # add metadata / covariates

print('PERMANOVA..')
set.seed(1234)
dfanova_meno <- menopause %>%
    slice(match(rownames(bray_meno), ID)) # distance matrix and metadata must have the same sample order
all(dfanova_meno$ID == rownames(bray_meno)) # TRUE
res2 <- adonis2(bray_meno ~ MenopauseYn, data = dfanova_meno) # PERMANOVA
print(res2)

print('plotting..')
braycurt_meno <- dbray_meno %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
    ggtitle("PCoA Bray-Curtis distance") +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    scale_color_manual(values = pal_nejm()(4)[3:4]) +
    theme_Publication() +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = MenopauseYn), type = "norm") +
    annotate("text", x= Inf, y = Inf, hjust = 1, vjust = 1,
             label = str_c("PERMANOVA: p = ", res2$`Pr(>F)`, ", r2 = ", 
                           format(round(res2$R2[1],3), nsmall = 3)))
ggsave(braycurt_meno, filename = "results/ordination/PCoA_BrayCurtis_meno_shotgun.pdf", 
       device = "pdf", width = 5, height = 5)

#### Weighted UniFrac ####
print('Weighted UniFrac')
wunifrac <- UniFrac(phydata, normalized = T, weighted = T, parallel = T)
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
head(expl_variance)
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, df_new, by = 'ID') # add metadata / covariates

print('plotting..')
unifracpl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    ggtitle('PCoA Weighted UniFrac') +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac.pdf", 
       device = "pdf", width = 8, height = 8)

print('plotting..')
unifrac_meno <- dfpc %>% filter(Sex == "Female") %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = pal_nejm()(4)[3:4]) +
    labs(color = "") +
    ggtitle('PCoA Weighted UniFrac') +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = MenopauseYn), type = "norm")
ggsave(unifrac_meno, filename = "results/ordination/PCoA_WeightedUnifrac_meno.pdf", 
       device = "pdf", width = 5, height = 5)

print('Weighted UniFrac')
pcoord <- ape::pcoa(wunifrac, correction = "cailliez")
expl_variance <- pcoord$values$Rel_corr_eig * 100
head(expl_variance)
dfpc <- pcoord$vectors[, c('Axis.1', 'Axis.2')] # get PCoA coordinates
dfpc <- as.data.frame(dfpc)
dfpc$ID <- rownames(dfpc)
dfpc <- left_join(dfpc, dfshot, by = 'ID') # add metadata / covariates

print('plotting..')
unifracpl <- dfpc %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = Sex), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = rev(pal_nejm()(2))) +
    ggtitle('PCoA Weighted UniFrac') +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    labs(color = "") +
    stat_ellipse(aes(color = Sex), type = "norm")
ggsave(unifracpl, filename = "results/ordination/PCoA_WeightedUnifrac_shotgun.pdf", 
       device = "pdf", width = 5, height = 5)

print('plotting..')
unifrac_meno <- dfpc %>% filter(Sex == "Female") %>% 
    ggplot(aes(Axis.1, Axis.2)) +
    geom_point(aes(color = MenopauseYn), size = 1, alpha = 0.7) +
    xlab(paste0('PCo1 (', round(expl_variance[1], digits = 1),'%)')) +
    ylab(paste0('PCo2 (', round(expl_variance[2], digits = 1),'%)')) +
    theme_Publication() +
    scale_color_manual(values = pal_nejm()(4)[3:4]) +
    labs(color = "") +
    ggtitle('PCoA Weighted UniFrac') +
    guides(fill = guide_legend(override.aes = list(shape = 21, size = 2))) +
    stat_ellipse(aes(color = MenopauseYn), type = "norm")
ggsave(unifrac_meno, filename = "results/ordination/PCoA_WeightedUnifrac_meno_shotgun.pdf", 
       device = "pdf", width = 5, height = 5)
