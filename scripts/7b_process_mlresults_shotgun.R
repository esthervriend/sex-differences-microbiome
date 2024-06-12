## Process results machine learning

library(tidyverse)
library(ggplot2)
library(ggsci)
library(stringr)

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

source("scripts/functions_shotgun.R")

#### Composition ####
## Plot assembled figure composition sex
path_true <- 'sex_metagen/output_XGB_class_sex_metagen_2024_03_11__12-49-23'
data_path <- 'sex_metagen/input_data'
labels <- c("Men", "Women")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), 
                    heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/comp_sex.pdf",
       width = 14, height = 18)
ggsave(plarr1, filename = "results/comp_sex.png",
       width = 14, height = 18)


## Plot assembled figure composition menopause
path_true <- 'menopause_metagen/output_XGB_class_menopause_metagen_2024_03_11__13-00-41'
data_path <- 'menopause_metagen/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_shotgun(path_true, 20)
plot_feature_importance_color_shotgun(path_true, 20)
plot_features_tests_shotgun(data_path, path_true, top_n=20, labels)
plot_features_top_shotgun(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_shotgun(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
pl3 <- plot_features_top_shotgun(data_path, path_true, top_n=5, nrow = 1, labels)
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, pl3,
                    nrow = 3, labels = c("A", "B", "C"), heights = c(1.2,1.3,0.8))
ggsave(plarr1, filename = "results/comp_menopause.pdf",
       width = 14, height = 15)
ggsave(plarr1, filename = "results/comp_menopause.png",
       width = 14, height = 15)


#### PATHWAYS ####
## Plot differences between sex pathways of interest
path_true <- 'sex_pathways/output_XGB_class_sex_pathways_2024_03_03__20-18-59'
data_path <- 'sex_pathways/input_data'
labels <- c("Men", "Women")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)

df <- readRDS('data/clinicaldata_shotgun.RDS') %>% mutate(Sex = fct_recode(Sex, "Men" = "Male", "Women" = "Female"))
pathwaysfilt <- readRDS("data/pathways_filtered.RDS")
pathwaysofinterest <- c("PPGPPMET-PWY", "PHOSLIPSYN-PWY", "PYRIDOXSYN-PWY", "PWY-5030")
pathwaysfilt <- pathwaysfilt %>% filter(rownames(.) %in% pathwaysofinterest)
pyri <- as.data.frame(t(as.matrix(pathwaysfilt)))
pyri$ID <- rownames(pyri)
pyri <- pyri %>% dplyr::select(everything(.), "ppGp metabolism" = `PPGPPMET-PWY`,
                        "pyridoxal 5'−phosphate biosynthesis I" = `PYRIDOXSYN-PWY`,
                        "superpathway of phospholipid biosynthesis I (bacteria)" = `PHOSLIPSYN-PWY`,
                        "L−histidine degradation III" = `PWY-5030`)

dfpaths <- inner_join(df, pyri, by = "ID")
dfpaths <- dfpaths %>% mutate(Menopause_Sex = case_when(
    Sex == "Men" & Age < 50 ~ "Men <50",
    Sex == "Men" & Age >= 50 ~ "Men >50",
    MenopauseYn == "Postmenopausal" ~ "Postmenopausal",
    MenopauseYn == "Premenopausal" ~ "Premenopausal"
),
Menopause_Sex = as.factor(Menopause_Sex),
Menopause_Sex = fct_relevel(Menopause_Sex, "Postmenopausal", after = 3L),
Menopause_Sex = fct_relevel(Menopause_Sex, "Younger men", after = 0L)
)

pathli <- list()
for(path in names(dfpaths)[59:62]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    lab_y <- log10(max(dfpaths$pathid))*1.05
    comp <- list(c("Men <50", "Men >50"), c("Men <50", "Premenopausal"), c("Premenopausal", "Postmenopausal"),
                 c("Men >50", "Postmenopausal"))
    pl <- ggplot(data = dfpaths, aes(x = Menopause_Sex, 
                                     y = pathid)) +
        scale_y_log10()+
        geom_violin(aes(fill = Menopause_Sex)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal_nejm()(4)[c(4,3,2,1)], guide = "none") +
        stat_compare_means(label.y = lab_y, comparisons = comp, label = "p.signif",
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
ggsave("results/diffpathways_menopause.pdf", width = 10, height = 12)
ggsave("results/diffpathways_menopause.png", width = 10, height = 12)

pathli2 <- list()
for(path in names(dfpaths)[59:62]){
    dfpaths$pathid <- dfpaths[[path]]
    print(path)
    lab_y <- log10(max(dfpaths$pathid))*1.2
    pl <- ggplot(data = dfpaths, aes(x = Sex, y = pathid)) +
        scale_y_log10()+
        geom_violin(aes(fill = Sex)) +
        geom_boxplot(fill = "white", width = 0.2, outlier.shape = NA) +
        scale_fill_manual(values = pal_nejm()(2)[c(2,1)], guide = "none") +
        stat_compare_means(method = "t.test") +
        labs(fill = "", x = "", y = "Relative abundance (cpm)", title = path) +
        theme_Publication() +
        theme(plot.title = element_text(size = rel(0.8)))
    print(pl)
    pathli2[[path]] <- pl
    dfpaths$pathid <- NULL
}

path_true <- 'sex_pathways/output_XGB_class_sex_pathways_2024_03_03__20-18-59'
pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr1 <- ggarrange(ggarrange(svg_grob, NULL, nrow = 2, heights = c(2.0,0.3)), pl2, 
                    nrow = 1, widths = c(0.8, 1.3), labels = c("A", "B"))
plarr2 <- ggarrange(plarr1, 
                    ggarrange(plotlist = pathli[c(2,1,4,3)], nrow = 1, ncol = 4, labels = LETTERS[3:7]),
                    nrow = 2, heights = c(1.0, 0.6))
plarr2b <- ggarrange(svg_grob, pl2, 
                    ggarrange(plotlist = pathli2[c(2,1,4,3)], nrow = 1, ncol = 4, labels = LETTERS[3:7]),
                    nrow = 3, heights = c(1.0, 1.2, 1.0), labels = c(LETTERS[1:2], ""))
ggsave(plarr2b, filename = "results/pathways_sexdiff.pdf",
       width = 15, height = 20)
ggsave(plarr2b, filename = "results/pathways_sexdiff.png",
       width = 15, height = 20)

## Plot assembled figure pathways menopause
path_true <- 'menopause_pathways/output_XGB_class_menopause_pathways_2024_03_03__20-28-16'
data_path <- 'menopause_pathways/input_data'
labels <- c("Postmenopausal", "Premenopausal")

plot_feature_importance_pathways(path_true, 20)
plot_features_tests_pathways(data_path, path_true, top_n=20, labels)
plot_features_top_pathways(data_path, path_true, top_n=20, nrow=4, labels)

pl2 <- plot_feature_importance_pathways(path_true, 20)
grConvert::convertPicture(file.path(path_true,"Plot_AUC.pdf"), file.path(path_true,"auc.svg"))
svg_grob <- svgparser::read_svg(file.path(path_true,"auc.svg"))
plarr1 <- ggarrange(ggarrange(svg_grob), pl2, 
                    nrow = 2, labels = c("A", "B"), heights = c(1.0,1.0))
ggsave(plarr1, filename = "results/pathways_menopause.pdf",
       width = 12, height = 15)
ggsave(plarr1, filename = "results/pathways_menopause.png",
       width = 12, height = 15)

