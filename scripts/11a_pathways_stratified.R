# Pathway plots: which species are driving pathway differences
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)

## Output folder
resultsfolder <- "results/pathways"
dir.create(resultsfolder, showWarnings = FALSE)

## Clinical data
df <- readRDS('data/clinicaldata_shotgun.RDS') %>% 
    mutate(Sex = fct_recode(Sex, "Men" = "Male", "Women" = "Female"),
           Sex = fct_rev(Sex)) # women left and men right in the plot (to align with other plots)

pathways <- rio::import("data/humann/merged_pathway_renorm_stratified.tsv")
names(pathways)[2:ncol(pathways)] <- str_remove(names(pathways)[2:ncol(pathways)], "_T1_Abundance-CPM")
rownames(pathways) <- pathways$`# Pathway`
pathways$`# Pathway` <- NULL
pathways <- as.data.frame(t(as.matrix(pathways)))
dim(pathways)
fi <- rio::import("data/sex_pathways_feature_importance.txt") %>% arrange(-RelFeatImp) %>% slice(1:5)
fi

#### Pyridoxal pathway ####
pyridox <- pathways %>% dplyr::select(contains(fi$FeatName[1]))
tk <- apply(pyridox, 2, function(x) sum(x > 10) > (0.2*length(x)))
pyridox <- pyridox[,tk] # apply
## tidy rownames
colnames(pyridox) <- str_remove_all(colnames(pyridox), "PYRIDOXSYN-PWY: pyridoxal 5'-phosphate biosynthesis I\\|")
colnames(pyridox) <- str_remove_all(colnames(pyridox), "g__[A-z]*.")

## Merge with clinical data
pyridox$ID <- rownames(pyridox)
dfpaths <- inner_join(df, pyridox, by = "ID")
dim(dfpaths)

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4, method = "wilcox.test", tip.length = 0)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 3, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/pyridoxalphos.pdf", width = 12, height = 15)
ggsave("results/pathways/pyridoxalphos.png", width = 12, height = 15)


#### ppGpp pathway ####
ppgpp <- pathways %>% dplyr::select(contains(fi$FeatName[2]))
tk <- apply(ppgpp, 2, function(x) sum(x > 5) > (0.1*length(x)))
ppgpp <- ppgpp[,tk] # apply
## tidy rownames
colnames(ppgpp) <- str_remove_all(colnames(ppgpp), "PPGPPMET-PWY: ppGpp metabolism\\|")
colnames(ppgpp) <- str_remove_all(colnames(ppgpp), "g__[A-z]*.")

## Merge with clinical data
ppgpp$ID <- rownames(ppgpp)
dfpaths <- inner_join(df, ppgpp, by = "ID")
dim(dfpaths)

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4, method = "wilcox.test", tip.length = 0)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 2, nrow = 1, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/ppgpp.pdf", width = 8, height = 5)
ggsave("results/pathways/ppgpp.png", width = 8, height = 5)


#### Stachyose pathway ####
stach <- pathways %>% dplyr::select(contains(fi$FeatName[3]))
tk <- apply(stach, 2, function(x) sum(x > 20) > (0.2*length(x)))
stach <- stach[,tk] # apply
## tidy rownames
colnames(stach) <- str_remove_all(colnames(stach), "PWY-6527: stachyose degradation\\|")
colnames(stach) <- str_remove_all(colnames(stach), "g__[A-z]*.")

## Merge with clinical data
stach$ID <- rownames(stach)
dfpaths <- inner_join(df, stach, by = "ID")
dim(dfpaths)

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4, method = "wilcox.test", tip.length = 0)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 4, nrow = 5, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/stachyose.pdf", width = 15, height = 22)
ggsave("results/pathways/stachyose.png", width = 15, height = 22)


#### Glucose-xylose pathway ####
gluc <- pathways %>% dplyr::select(contains(fi$FeatName[4]))
tk <- apply(gluc, 2, function(x) sum(x > 5) > (0.1*length(x)))
gluc <- gluc[,tk] # apply
## tidy rownames
colnames(gluc) <- str_remove_all(colnames(gluc), "PWY-6901: superpathway of glucose and xylose degradation\\|")
colnames(gluc) <- str_remove_all(colnames(gluc), "g__[A-z]*.")

## Merge with clinical data
gluc$ID <- rownames(gluc)
dfpaths <- inner_join(df, gluc, by = "ID")
dim(dfpaths)

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4, method = "wilcox.test", tip.length = 0)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 1, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/glucxyl.pdf", width = 8, height = 5)
ggsave("results/pathways/glucxyl.png", width = 8, height = 5)

#### Phospholipid pathway ####
phos <- pathways %>% dplyr::select(contains(fi$FeatName[5]))
tk <- apply(phos, 2, function(x) sum(x > 5) > (0.1*length(x)))
phos <- phos[,tk] # apply
## tidy rownames
colnames(phos) <- str_remove_all(colnames(phos), 
                                 "PHOSLIPSYN-PWY: superpathway of phospholipid biosynthesis I \\(bacteria\\)\\|")
colnames(phos) <- str_remove_all(colnames(phos), "g__[A-z]*.")

## Merge with clinical data
phos$ID <- rownames(phos)
dfpaths <- inner_join(df, phos, by = "ID")
dim(dfpaths)

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_nejm(guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance (cpm)', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4, method = "wilcox.test", tip.length = 0)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 2, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/phos.pdf", width = 10, height = 10)
ggsave("results/pathways/phos.png", width = 10, height = 10)
