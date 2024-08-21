# Pathway plots: which species are driving pathway differences
# Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

# Libraries
library(tidyverse)
library(ggsci)
library(ggpubr)

## Clinical data
df <- readRDS('data/clinicaldata_shotgun.RDS') %>% 
    mutate(Sex = fct_recode(Sex, "Men" = "Male", "Women" = "Female"))
pathwaysfilt <- readRDS("data/pathways_filtered.RDS")
fi <- rio::import("data/sex_pathways_feature_importance.txt") %>% arrange(-RelFeatImp) %>% slice(1:5)
pathwaysofinterest <- fi$FeatName
pathwaysfilt <- pathwaysfilt %>% filter(rownames(.) %in% pathwaysofinterest)
pyri <- as.data.frame(t(as.matrix(pathwaysfilt)))
pyri$ID <- rownames(pyri)
pyri <- pyri %>% dplyr::select(everything(.),"pyridoxal 5'−phosphate biosynthesis I" = `PYRIDOXSYN-PWY`, 
                               "ppGp metabolism" = `PPGPPMET-PWY`,
                               "stachyose degradation" = `PWY-6527`,
                               "superpathway of glucose and xylose degradation" = `PWY-6901`,
                               "superpathway of phospholipid biosynthesis I (bacteria)" = `PHOSLIPSYN-PWY`)
dfpaths <- inner_join(df, pyri, by = "ID")
dfpaths <- dfpaths %>% mutate(Menopause_Sex = case_when(
    Sex == "Men" & Age < 50 ~ "Men <50",
    Sex == "Men" & Age >= 50 ~ "Men >=50",
    MenopauseYn == "Postmenopausal" ~ "Postmenopausal",
    MenopauseYn == "Premenopausal" ~ "Premenopausal"
),
Menopause_Sex = as.factor(Menopause_Sex),
Menopause_Sex = fct_relevel(Menopause_Sex, "Postmenopausal", after = 3L),
Menopause_Sex = fct_relevel(Menopause_Sex, "Men <50", after = 0L)
)

pathli <- list()
for(path in names(dfpaths)[59:62]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    lab_y <- log10(max(dfpaths$pathid))*1.05
    comp <- list(c("Men <50", "Men >=50"), 
                 c("Men <50", "Premenopausal"), 
                 c("Premenopausal", "Postmenopausal"),
                 c("Men >=50", "Postmenopausal"))
    pl <- ggplot(data = dfpaths, aes(x = Menopause_Sex, 
                                     y = pathid)) +
        scale_y_log10()+
        geom_violin(aes(fill = Menopause_Sex)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal_nejm()(4)[c(4,3,2,1)], guide = "none") +
        stat_compare_means(label.y = lab_y, comparisons = comp, method = "wilcox.test",
                           # label = "p.signif",
                           step.increase = 0.07, tip.length = 0) +
        labs(fill = "", x = "", y = "Relative abundance (cpm)", title = path) +
        theme_Publication() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1), 
              plot.title = element_text(size = rel(0.8)))
    print(pl)
    pathli[[path]] <- pl
    dfpaths$pathid <- NULL
}
ggarrange(plotlist = pathli, ncol = 2, nrow = 2, labels = LETTERS[1:4])
ggsave("results/pathways/diffpathways_menopause.pdf", width = 10, height = 12)
ggsave("results/pathways/diffpathways_menopause.png", width = 10, height = 12)

dfpaths$Sex <- fct_rev(dfpaths$Sex)
pathli2 <- list()
for(path in names(dfpaths)[59:62]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    lab_y <- log10(max(dfpaths$pathid))*1.2
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        scale_y_log10()+
        geom_violin(aes(fill = Sex)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal_nejm()(2), guide = "none") +
        stat_compare_means(method = "t.test") +
        labs(fill = "", x = "", y = "Relative abundance (cpm)", title = path) +
        theme_Publication() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pathli2[[path]] <- pl
    dfpaths$pathid <- NULL
}