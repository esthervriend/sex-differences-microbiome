#### Dietary data analysis

## Libraries
library(tableone)
library(rio)
library(dplyr)
library(mixOmics)

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

## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")


#### PCA sex differences in diet ####
df_diet <- df_new %>% dplyr::select(ID, Sex, MenopauseYn, TotalCalories, Fibre, Proteins, 
                                    Carbohydrates, FattyAcids, SatFat, AlcoholIntake, Sodium) %>% 
    mutate(Sex = fct_recode(Sex, "Men" = "Male", "Women" = "Female")) %>% 
    filter(!is.na(TotalCalories)) %>% 
    filter(TotalCalories < 7500)
df_diet2 <- df_diet %>% dplyr::select(-ID, -Sex, -MenopauseYn, -TotalCalories) %>% mutate(across(everything(.), scale))
matdiet <- as.matrix(df_diet2)
tunediet <- tune.pca(matdiet, ncomp = 5, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matdiet, ncomp = 4)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_diet$ID, Sex = df_diet$Sex, MenopauseYn = df_diet$MenopauseYn)
expvar_diet <- pc$prop_expl_var$X[1:4]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)

(pcadiet <- pcs %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = Sex), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = Sex, fill = Sex), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = rev(pal_nejm()(2))) +
        scale_fill_manual(values = rev(pal_nejm()(2)), guide = "none") +
        labs(color = "", title = "Sex differences diet")
        # geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
        #              arrow = arrow(length = unit(1/2, "picas")),
        #              color = "black", linewidth = 0.9) +
        # annotate("text", x = (loadings$PC1*9), y = (loadings$PC2*9),
        #          label = loadings$Variables)
)
ggsave(pcadiet, filename = "results/diet/PCA_diet.pdf",
       width = 5, height = 6)
ggsave(pcadiet, filename = "results/diet/PCA_diet.png",
       width = 5, height = 6)

colpal <- rev(pal_nejm()(2))
vars <- c("TotalCalories", "Proteins", "Fibre", "Carbohydrates", "FattyAcids", "SatFat", "Sodium",
          "AlcoholIntake")
plist <- c()
for(a in vars){
    df_diet$var <- df_diet[[a]]
    unit <- case_when(a == "TotalCalories" ~ "Kcal", .default = "gram")
    pl <- ggplot(df_diet, aes(x=Sex, y=var)) +
        geom_violin(aes(fill=Sex)) +
        scale_fill_manual(values = colpal, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
        theme_Publication() +
        theme(legend.position = 'none') +
        labs(x='', y = unit, title = a) +
        ggpubr::stat_compare_means(method = "wilcox.test", tip.length = 0, label="p.signif", hide.ns = TRUE,
                                   comparisons = list(c(levels(df_diet$Sex))))
    plist[[a]] <- pl
}
pldiff <- ggarrange(plotlist=plist, labels = LETTERS[1:length(plist)], ncol = 4, nrow = 2)
ggarrange(pldiff, ggarrange(pcadiet, NULL), labels = c("", LETTERS[length(plist)+1]), nrow = 2, heights = c(2,1))
ggsave("results/diet/sex_diff_pca.pdf", height = 14, width = 12)
ggsave("results/diet/sex_diff_pca.png", height = 14, width = 12)


#### PCA menopausal differences diet ####
df_diet1 <- df_diet %>% filter(!is.na(MenopauseYn))
df_diet2 <- df_diet1 %>% 
    dplyr::select(-ID, -Sex, -MenopauseYn, -TotalCalories) %>% mutate(across(everything(.), scale))
matdiet <- as.matrix(df_diet2)
tunediet <- tune.pca(matdiet, ncomp = 5, scale = TRUE)
plot(tunediet)
pc <- mixOmics::pca(matdiet, ncomp = 2)
pcs <- as.data.frame(pc$variates$X)
pcs <- pcs %>% mutate(ID = df_diet1$ID, Sex = df_diet1$Sex, MenopauseYn = df_diet1$MenopauseYn)
expvar_diet <- pc$prop_expl_var$X[1:4]
loadings <- as.data.frame(pc$loadings$X)
loadings$Variables <- rownames(loadings)
(pcadiet <- pcs %>% filter(!is.na(MenopauseYn)) %>% 
        ggplot(aes(PC1, PC2)) +
        geom_point(aes(color = MenopauseYn), size = 1, alpha = 1.0) +
        xlab(paste0('PC1 (', round(expvar_diet[1]*100, digits = 1),'%)')) +
        ylab(paste0('PC2 (', round(expvar_diet[2]*100, digits = 1),'%)')) +
        theme_Publication() +
        stat_ellipse(geom = "polygon", aes(color = MenopauseYn, fill = MenopauseYn), linewidth = 1.0,
                     alpha = 0.1, type = "norm")+
        scale_color_manual(values = pal_nejm()(4)[3:4]) +
        scale_fill_manual(values = pal_nejm()(4)[3:4], guide = "none") +
        labs(color = "", title = "Menopausal differences diet")
        # geom_segment(data = loadings, aes(x = 0, y = 0, xend = (PC1*8), yend = (PC2*8)), 
        #              arrow = arrow(length = unit(1/2, "picas")),
        #              color = "black", linewidth = 0.9) +
        # annotate("text", x = (loadings$PC1*9), y = (loadings$PC2*9),
        #          label = loadings$Variables)
)
ggsave(pcadiet, filename = "results/diet/PCA_diet_menopause.pdf",
       width = 5, height = 6)
ggsave(pcadiet, filename = "results/diet/PCA_diet_menopause.png",
       width = 5, height = 6)

colpal <- pal_nejm()(4)[3:4]
vars <- c("TotalCalories", "Proteins", "Fibre", "Carbohydrates", "FattyAcids", "SatFat", "Sodium",
          "AlcoholIntake")
plist <- c()
for(a in vars){
    df_diet1$var <- df_diet1[[a]]
    unit <- case_when(a == "TotalCalories" ~ "Kcal", .default = "gram")
    pl <- ggplot(df_diet1, aes(x=MenopauseYn, y=var)) +
        geom_violin(aes(fill=MenopauseYn)) +
        scale_fill_manual(values = colpal, guide = FALSE) +
        geom_boxplot(width=0.1, fill="white", outlier.shape = NA) +
        theme_Publication() +
        theme(legend.position = 'none') +
        labs(x='', y = unit, title = a) +
        ggpubr::stat_compare_means(method = "wilcox.test", tip.length = 0, label="p.signif", hide.ns = TRUE,
                                   comparisons = list(c(levels(df_diet$MenopauseYn))))
    plist[[a]] <- pl
}
pldiff <- ggarrange(plotlist=plist, labels = LETTERS[1:length(plist)], ncol = 4, nrow = 2)
ggarrange(pldiff, ggarrange(pcadiet, NULL), labels = c("", LETTERS[length(plist)+1]), nrow = 2, heights = c(2,1))
ggsave("results/diet/menopause_diff_pca.pdf", height = 14, width = 13)
ggsave("results/diet/menopause_diff_pca.png", height = 14, width = 13)
