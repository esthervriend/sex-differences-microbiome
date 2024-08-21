# Regression analyses pathways - Sex

## Libraries
library(tidyverse)
library(broom)
library(ggpubr)
library(ggsci)

## Functions
theme_Publication <- function(base_size=14, base_family="sans") {
  library(grid)
  library(ggthemes)
  library(stringr)
  (theme_foundation(base_size=base_size, base_family=base_family)
    + theme(plot.title = element_text(face = "bold",
                                      size = rel(1.2), hjust = 0.5),
            text = element_text(),
            panel.background = element_rect(colour = NA),
            plot.background = element_rect(colour = NA),
            panel.border = element_rect(colour = NA),
            axis.title = element_text(face = "bold",size = rel(1)),
            axis.title.y = element_text(angle=90,vjust =2),
            axis.title.x = element_text(vjust = -0.2),
            axis.text = element_text(), 
            axis.line = element_line(colour="black"),
            axis.ticks = element_line(),
            panel.grid.major = element_line(colour="#f0f0f0"),
            panel.grid.minor = element_blank(),
            legend.key = element_rect(colour = NA),
            legend.position = "right",
            # legend.direction = "horizontal",
            legend.key.size= unit(0.2, "cm"),
            legend.spacing  = unit(0, "cm"),
            # legend.title = element_text(face="italic"),
            plot.margin=unit(c(10,5,5,5),"mm"),
            strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
            strip.text = element_text(face="bold")
    ))}
afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
afronden3 <- function(x) return(as.numeric(format(round(x, 3),3)))

# Load data
df <- readRDS("data/clinicaldata.RDS")
tab <- readRDS("data/pathways_filtered.RDS")
tab$ID <- rownames(tab)
dfcomplete <- merge(df, tab, by = "ID")
tax <- readRDS("data/pathway_keys.RDS")

# Best predictors
pred_sex_pathways <- rio::import("data/sex_pathways_feature_importance.txt") %>% slice_head(n = 20)

# Regression analyses sex
res <- c()
for (i in 1:nrow(pred_sex_pathways)) {
    var_name <- pred_sex_pathways$FeatName[i]
    print(var_name) # print name
    dfcomplete$dep <- scale(dfcomplete[[var_name]])
    dfcomplete$logdep <- log(dfcomplete$dep + 1)
    # inspect distribution
    print(dfcomplete %>% gghistogram("dep", title = str_c(var_name, " - no log"), fill = "royalblue4"))
    print(dfcomplete %>% gghistogram("logdep", title = str_c(var_name, " - log+1"), fill = "firebrick"))
    
    # run models
    m0 <- lm(dep ~ Sex, data = dfcomplete)
    m1 <- lm(dep ~ Sex + Age + BMI + HT + DM + CurrSmoking, data = dfcomplete)
    m2 <- lm(dep ~ Sex + Age + BMI + HT + DM + CurrSmoking + Alcohol + TotalCalories + 
                 Fibre + Proteins, data = dfcomplete)
    
    # extract estimates for variable sex
    m0 <- tidy(m0, conf.int=T)[2,]
    m1 <- tidy(m1, conf.int=T)[2,]
    m2 <- tidy(m2, conf.int = T)[2,]
    
    # bind results together
    resRow <- cbind(var_name, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                    m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                    m2$estimate, m2$conf.low, m2$conf.high, m2$p.value)
    colnames(resRow) <- c("pathway", 
                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                          "m2-est", "m2-l95", "m2-u95", "m2-p")
    
    # add results row to final table
    res <- rbind(res, resRow)
    dfcomplete$dep <- NULL
}

res2 <- as.data.frame(res) %>% 
    mutate(across(c(2:13), ~as.numeric(as.character(.x)))) %>% # first to char, then num
    mutate(`m0-q` = p.adjust(`m0-p`, 'fdr'), # p adjust 
           `m1-q` = p.adjust(`m1-p`, 'fdr'),
           `m2-q` = p.adjust(`m2-p`, 'fdr')) %>% 
    mutate(across(c(5,9,13), afronden3), # round p values at 3 digits
           across(c(2:4, 6:8, 10:12), afronden2)) # round estimates at 2 digits
openxlsx::write.xlsx(res2, "results/lm/lm_sex_pathways.xlsx")

res3 <- res2 %>% pivot_longer(c(2:16), names_to=c("model", "cat"), names_prefix="m", 
                              names_sep='-', values_to="value") %>% # to have models in long format
    pivot_wider(names_from = cat, values_from = value) %>% # but est+l95/u95 to wide
    mutate(model = factor(model, levels = c("0", "1", "2"), # label models
                          labels = c("Unadjusted", "+Age, BMI, HT, DM, Smoking", "+Diet")),
           model = fct_inorder(model),
           pathway = tax$label[match(pathway, tax$key)], # change name to tidy tax
           pathway = fct_rev(fct_inorder(pathway)),
           sigq = case_when(q < 0.05 ~ paste0("q<0.05"), q >= 0.05 ~ paste0("not sig")),
           sigq = as.factor(sigq)
    )

# Make forest plot
(forest_path_sex <- ggplot(res3, aes(x = est, y = pathway, color = model, shape = sigq)) +
        geom_point(position = position_dodge(-0.5), size = 2.0) +
        geom_vline(aes(xintercept = 0), linewidth = .50, linetype = "dashed") + 
        geom_errorbarh(aes(xmin = l95, xmax = u95), height = 0.5, 
                       position = position_dodge(-0.5)) +
        scale_x_continuous(breaks = seq(-60, 60, by = 20)) +
        scale_color_nejm() +
        scale_shape_manual(values = c(1, 16)) +
        labs(title = "Best predicting pathways for sex",
             x = "Estimate and 95% CI for women", y = "", color = "", shape = "") + 
        theme_Publication()
)
ggsave("results/lm/lm_pathways_sex.pdf", width = 12, height = 8)
