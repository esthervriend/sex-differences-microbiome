#### Pathway plots

## Libraries
library(tidyverse)
library(ggplot2)
library(ggsci)
library(ggpubr)

## Data
df <- readRDS('data/clinicaldata_shotgun.RDS') %>% mutate(Sex = fct_recode(Sex, "Men" = "Male", "Women" = "Female"))
pathways <- rio::import("data/pyridoxal_pathways.txt")
head(pathways)
names(pathways)
dim(pathways)
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
names(pathways)
rownames(pathways)
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))]
pathways <- pathways[which(str_detect(rownames(pathways), "PYRIDOXSYN")),]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x)))
pathways <- pathways[tk,]
dim(pathways)
rownames(pathways) <- str_remove_all(rownames(pathways), "pyridoxal 5'-phosphate biosynthesis I\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
dim(dfpaths)

colorguide <- c(pal_nejm()(2)[c(2,1)])

pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4)+
        scale_y_log10() + 
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 5, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)    
ggsave("results/pyridoxalphos.pdf", width = 12, height = 22)
ggsave("results/pyridoxalphos.png", width = 12, height = 22)

abundance <- rio::import("data/metaphlan/merged/combined_table.tsv")
abundance2 <- abundance %>% filter(clade_name == "k__Bacteria|p__Verrucomicrobia|c__Verrucomicrobiae|o__Verrucomicrobiales|f__Akkermansiaceae|g__Akkermansia|s__Akkermansia_muciniphila")
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:298], "_T1", ""))
rownames(abundance2) <- abundance2$clade_name
abundance2$clade_name <- NULL
akk <- as.data.frame(t(as.matrix(abundance2)))
colnames(akk)[1] <- "Akkermansia_muciniphila"
akk$ID <- rownames(akk)

pathwaysfilt <- readRDS("data/pathways_filtered.RDS")
pathwaysofinterest <- c("PPGPPMET-PWY", "PHOSLIPSYN-PWY", "PYRIDOXSYN-PWY", "PWY-5030")
pathwaysfilt <- pathwaysfilt %>% filter(rownames(.) %in% pathwaysofinterest)
pyri <- as.data.frame(t(as.matrix(pathwaysfilt)))
pyri$ID <- rownames(pyri)
pyri <- pyri %>% select(everything(.), "ppGp metabolism" = `PPGPPMET-PWY`,
                                        "pyridoxal 5'−phosphate biosynthesis I" = `PYRIDOXSYN-PWY`,
                                        "superpathway of phospholipid biosynthesis I (bacteria)" = `PHOSLIPSYN-PWY`,
                                        "L−histidine degradation III" = `PWY-5030`)

dfpaths <- inner_join(dfpaths, akk, by = "ID") %>% inner_join(., pyri, by = "ID")
dfpaths <- dfpaths %>% mutate(Menopause_Sex = case_when(
        Sex == "Male" & Age < 50 ~ "Younger men",
        Sex == "Male" & Age >= 50 ~ "Older men",
        MenopauseYn == "Postmenopausal" ~ "Postmenopausal",
        MenopauseYn == "Premenopausal" ~ "Premenopausal"
        ),
        Menopause_Sex = as.factor(Menopause_Sex),
        Menopause_Sex = fct_relevel(Menopause_Sex, "Postmenopausal", after = 2L)
    )
names(dfpaths)
ps1 <- min(dfpaths$`pyridoxal 5'−phosphate biosynthesis I`[which(dfpaths$`pyridoxal 5'−phosphate biosynthesis I` != 0)])/2
ps2 <- min(dfpaths$Akkermansia_muciniphila[which(dfpaths$Akkermansia_muciniphila != 0)])/2
dfpaths$`pyridoxal 5'−phosphate biosynthesis I` <- dfpaths$`pyridoxal 5'−phosphate biosynthesis I` + ps1
dfpaths$Akkermansia_muciniphila <- dfpaths$Akkermansia_muciniphila + ps2
options(scipen=999)
ggplot(data = dfpaths %>% filter(Akkermansia_muciniphila  > ps2), 
       aes(x = `pyridoxal 5'−phosphate biosynthesis I`, 
                           y = Akkermansia_muciniphila)) +
    scale_x_log10()+
    scale_y_log10()+
    geom_jitter(aes(color = Sex)) +
    geom_smooth(method = "lm", aes(color = Sex)) +
    scale_color_manual(values = colorguide, guide = "none") +
    stat_cor() +
    facet_wrap(~Sex) +
    theme_Publication()
ggsave("results/akkermansia_pyridoxal.pdf", width = 7, height = 5)
ggsave("results/akkermansia_pyridoxal.png", width = 7, height = 5)

ggplot(data = dfpaths %>% filter(Akkermansia_muciniphila  > ps2), 
       aes(x = `pyridoxal 5'−phosphate biosynthesis I`, 
                           y = Akkermansia_muciniphila)) +
    scale_x_log10()+
    scale_y_log10()+
    geom_jitter(aes(color = MenopauseYn)) +
    geom_smooth(method = "lm", aes(color = MenopauseYn)) +
    scale_color_manual(values = pal_nejm()(4)[c(4,3)]) +
    stat_cor() +
    facet_wrap(~Sex) +
    labs(color = "") +
    theme_Publication()
ggsave("results/akkermansia_pyridoxal_meno.pdf", width = 7, height = 5)
ggsave("results/akkermansia_pyridoxal_meno.png", width = 7, height = 5)

#### Other pathways ####
pathways <- rio::import("data/histidine_pathways.txt")
head(pathways)
names(pathways)
dim(pathways)
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
names(pathways)
rownames(pathways)
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x)))
pathways <- pathways[tk,]
dim(pathways)
rownames(pathways) <- str_remove_all(rownames(pathways), "L-histidine degradation III\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
dim(dfpaths)

colorguide <- c(pal_nejm()(2)[c(2,1)])

pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4)+
        scale_y_log10() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 4, nrow = 5, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)  
ggsave("results/histidine_pathways.pdf", width = 13, height = 25)
ggsave("results/histidine_pathways.png", width = 13, height = 25)


pathways <- rio::import("data/arginine_pathways.txt")
head(pathways)
names(pathways)
dim(pathways)
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
names(pathways)
rownames(pathways)
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x)))
pathways <- pathways[tk,]
dim(pathways)
rownames(pathways) <- str_remove_all(rownames(pathways), "L-arginine biosynthesis III \\(via N-acetyl-L-citrulline\\)\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
dim(dfpaths)

colorguide <- c(pal_nejm()(2)[c(2,1)])

pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4)+
        scale_y_log10() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 3, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)  
ggsave("results/arginine_pathways.pdf", width = 11, height = 14)


pathways <- rio::import("data/pgpp_pathways.txt")
head(pathways)
names(pathways)
dim(pathways)
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
names(pathways)
rownames(pathways)
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x)))
pathways <- pathways[tk,]
dim(pathways)
rownames(pathways) <- str_remove_all(rownames(pathways), "ppGpp metabolism\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
dim(dfpaths)

colorguide <- c(pal_nejm()(2)[c(2,1)])

pltli <- list()
for(path in names(dfpaths)[59:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4)+
        scale_y_log10() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 2, nrow = 1, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)  
ggsave("results/pgpp_pathways.pdf", width = 7, height = 5)


pathways <- rio::import("data/lacteth_pathways.txt")
head(pathways)
names(pathways)
dim(pathways)
colnames(pathways) <- str_remove_all(names(pathways), "_T1_Abundance")
rownames(pathways) <- pathways$`# Gene Family`
pathways$`# Gene Family` <- NULL
names(pathways)
rownames(pathways)
pathways <- pathways[,which(!str_detect(names(pathways), "-RPK"))]
tk <- apply(pathways, 1, function(x) sum(x > 0.5) > (0.1*length(x)))
pathways <- pathways[tk,]
dim(pathways)
rownames(pathways) <- str_remove_all(rownames(pathways), "hexitol fermentation to lactate, formate, ethanol and acetate\\|")
rownames(pathways) <- str_remove_all(rownames(pathways), "g__[A-z]*.")

pathways2 <- as.data.frame(t(as.matrix(pathways)))
pathways2$ID <- rownames(pathways2)
dfpaths <- inner_join(df, pathways2, by = "ID")
dim(dfpaths)

colorguide <- c(pal_nejm()(2)[c(2,1)])

pltli <- list()
for(path in names(dfpaths)[58:ncol(dfpaths)]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        geom_violin(aes(fill=Sex), trim = TRUE)+
        scale_fill_manual(values = colorguide, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white")+
        theme_Publication()+
        theme(legend.position = 'none')+
        labs(x='Group', y = 'Relative abundance', title = path)+
        ggpubr::stat_compare_means(comparisons = list(c("Men", "Women")), label = "p.signif",
                                   paired = F, size = 4)+
        scale_y_log10() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pltli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pltli, ncol = 3, nrow = 1, labels = LETTERS[1:length(pltli)],
          common.legend = TRUE)  
ggsave("results/lacteth_pathways.pdf", width = 9, height = 5)



    