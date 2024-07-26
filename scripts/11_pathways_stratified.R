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

#### Pyridoxal pathway ####
## Open pathway data and tidy
pathways <- rio::import("data/pyridoxal_pathways.txt")
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))] # select cpm
pathways <- pathways[which(str_detect(rownames(pathways), "PYRIDOXSYN")),]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x))) # filter >0.5 in 10%
pathways <- pathways[tk,] # apply
## tidy rownames
rownames(pathways) <- str_remove_all(rownames(pathways), "pyridoxal 5'-phosphate biosynthesis I\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

## Merge with clinical data
pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
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
ggarrange(plotlist = pltli, ncol = 3, nrow = 5, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)
ggsave("results/pathways/pyridoxalphos.pdf", width = 12, height = 22)
ggsave("results/pathways/pyridoxalphos.png", width = 12, height = 22)


#### Histidine pathway ####
## Open pathway data and tidy
pathways <- rio::import("data/histidine_pathways.txt")
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))] # select cpm
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x))) # filter >0.5 in 10%
pathways <- pathways[tk,] # apply filter
rownames(pathways) <- str_remove_all(rownames(pathways), "L-histidine degradation III\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

## Merge with clinical data
pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")

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
ggsave("results/pathways/histidine_pathways.pdf", width = 13, height = 25)
ggsave("results/pathways/histidine_pathways.png", width = 13, height = 25)


#### Arginine pathway ####
## Open pathway data and tidy
pathways <- rio::import("data/arginine_pathways.txt")
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))] # select cpm
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x))) # filter >0.5 in 10%
pathways <- pathways[tk,] # apply filter
rownames(pathways) <- str_remove_all(rownames(pathways), "L-arginine biosynthesis III \\(via N-acetyl-L-citrulline\\)\\|") # because they all start with this
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

## Merge with clinical data
pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")

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
ggsave("results/pathways/arginine_pathways.pdf", width = 11, height = 14)

#### pGpp pathway ####
## Open pathway data and tidy
pathways <- rio::import("data/pgpp_pathways.txt")
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))] # select cpm
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x))) # filter
pathways <- pathways[tk,] # apply filter
rownames(pathways) <- str_remove_all(rownames(pathways), "ppGpp metabolism\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

## Merge with clinical data
pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")

## Make violin plots
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
ggsave("results/pathways/pgpp_pathways.pdf", width = 7, height = 5)


#### Lactate/ethnanol pathway ####
## Open pathway data and tidy
pathways <- rio::import("data/lacteth_pathways.txt")
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))] # select cpm
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x))) # filter >0.5 in 10%
pathways <- pathways[tk,] # apply
rownames(pathways) <- str_remove_all(rownames(pathways), "hexitol fermentation to lactate, formate, ethanol and acetate\\|") # because they all start with this
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

## Merge with clinical data
pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")

## Make list of violin plots
pltli <- list()
for(path in names(dfpaths)[58:ncol(dfpaths)]){
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
ggarrange(plotlist = pltli, ncol = 4, nrow = 2, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)  
ggsave("results/pathways/lacteth_pathways.pdf", width = 9, height = 5)
    
