## Compositional plots and diversity metrics

## Libraries
library(phyloseq)
library(tidyverse)
library(ggpubr)
library(vegan)
library(circlize)

## Functions
cols <- colorRampPalette(c(pal_npg()(10)))

theme_composition <- function(base_size=14, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold",
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text.x =  element_text(angle = 45, hjust = 1),
                axis.text = element_text(), 
                axis.line = element_line(colour="black"),
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "right",
                legend.key.size= unit(0.4, "cm"),
                legend.spacing  = unit(0, "cm"),
                legend.text = element_text(size = rel(0.7)),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
}

## Output folder
resultsfolder <- "results/descriptives"
dir.create(resultsfolder, showWarnings = FALSE)

## Open data
phyloseqrare <- readRDS("data/phyloseq_sampledata.RDS")
tax <- readRDS("data/taxtable_rarefied.RDS")

tab <- as(phyloseqrare@otu_table, 'matrix')
counts <- sample_sums(phyloseqrare@otu_table)
tab <- as.data.frame(t(tab/sample_sums(phyloseqrare))*100)
rowSums(tab) # samples should all sum up to 100%
clindata <- readRDS("data/clinicaldata.RDS")

#### Genus-level ####
N <- 20
d <- tab %>% 
    rownames_to_column(var = 'Sample') %>% 
    pivot_longer(-Sample, names_to = 'ASV', values_to = 'Abundance') %>% 
    mutate(sampleID = Sample)
d$Genus <- tax$Genus[match(d$ASV, tax$ASV)]

top_gen <- d %>%
    group_by(Genus, Sample) %>%
    summarise(Abundance = sum(Abundance)) %>%
    group_by(Genus) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    arrange(-Abundance) %>%
    dplyr::select(Genus) %>%
    filter(Genus != 'ambiguous') %>%
    head(N) %>%
    unlist()
top_gen

dx_genus <- d %>% mutate(
        Genus2 = case_when(
            Genus %in% top_gen ~ paste(Genus),
            is.na(Genus) ~ paste("Unknown"),
            !(Genus %in% top_gen) ~ paste("Other genera")
        ),
        Genus2 = as.factor(Genus2)
    ) %>% 
    group_by(Genus2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(Sex = clindata$Sex[match(Sample, clindata$ID)]) %>% 
    group_by(Genus2, Sex) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    mutate(group = "all samples",
           Genus2 = fct_reorder(Genus2, Abundance),
           Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
           Genus2 = fct_relevel(Genus2, "Unknown", after = 0L)
           )

lev <- levels(dx_genus$Genus2)

set.seed(1234)
comp_genus <- dx_genus %>% 
    ggplot(aes(x = Sex, y = Abundance, fill = Genus2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(cols(20)), "grey70", "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (%)", x = "", title = "Sex", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(comp_genus, filename = "results/descriptives/composition_sex.pdf", width = 8, height = 8)
ggsave(comp_genus, filename = "results/descriptives/composition_sex.svg", width = 8, height = 8)

dx_genus <- d %>% mutate(
    Genus2 = case_when(
        Genus %in% top_gen ~ paste(Genus),
        is.na(Genus) ~ paste("Unknown"),
        !(Genus %in% top_gen) ~ paste("Other genera")
    ),
    Genus2 = as.factor(Genus2)
) %>% 
    group_by(Genus2, Sample) %>% 
    summarise(Abundance = sum(Abundance)) %>% 
    mutate(MenopauseYn = clindata$MenopauseYn[match(Sample, clindata$ID)],
           MenopauseYn = as.factor(MenopauseYn)) %>% 
    filter(!is.na(MenopauseYn)) %>% 
    group_by(Genus2, MenopauseYn) %>% 
    summarise(Abundance = mean(Abundance)) %>% 
    mutate(group = "all samples",
           Genus2 = fct_reorder(Genus2, Abundance),
           Genus2 = fct_relevel(Genus2, "Other genera", after = 0L),
           Genus2 = fct_relevel(Genus2, "Unknown", after = 0L)
    ) 

lev <- levels(dx_genus$Genus2)

set.seed(1234)
comp_genus <- dx_genus %>% 
    ggplot(aes(x = MenopauseYn, y = Abundance, fill = Genus2)) +
    geom_bar(stat = "identity", color = "black") +
    scale_fill_manual(values = rev(c(sample(cols(20)), "grey70", "grey90")), labels = lev) +
    guides(fill = guide_legend(ncol = 1)) +
    labs(y="Composition (%)", x = "Menopausal status", title = "Menopause", fill = "") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_composition()
ggsave(comp_genus, filename = "results/descriptives/composition_genus_meno.pdf", width = 8, height = 8)
ggsave(comp_genus, filename = "results/descriptives/composition_genus_meno.svg", width = 8, height = 8)
