#### Data cleaning

## Libraries
library(tidyverse)
library(rio) # kun je alles mee openen zonder te definieren wat voor bestand het is
library(ggsci) # voor kleurensets van de lancet, nejm en jama
library(mixOmics)
library(ggpubr) # voor statistiek bij plots
library(vegan) # voor diversiteitsindices
library(dplyr)
library(haven)
library(phyloseq) # voor handling van phyloseq object

## Open phyloseq object
otu <- rio::import("data/HELUIS.Verhaar.2023.otu_table.csv", header = TRUE)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax <- rio::import("data/HELUIS.Verhaar.2023.tax_table.csv", header = TRUE)
rownames(tax) <- tax$V1
tax$V1 <- NULL
tax <- as.matrix(tax)
tree <- ape::read.tree("data/HELUIS.Verhaar.2023.phy_tree.tree")

## Make phyloseq
heliusmb <- phyloseq(otu_table(otu, taxa_are_rows = T), tax_table(tax), phy_tree(tree))

## Inspect phyloseq object
ntaxa(heliusmb)
nsamples(heliusmb)
sample_names(heliusmb)
depth <- colSums(heliusmb@otu_table)
all(depth == 14932)
as.data.frame(as(heliusmb@otu_table, "matrix"))

### Fix ASV taxonomy
tax <- as.data.frame(as(heliusmb@tax_table, 'matrix'))
head(tax)
sum(!is.na(tax$Species)) / nrow(tax) * 100 # 2.0 % of all ASVs have species level
sum(!is.na(tax$Genus)) / nrow(tax) * 100 # 81.4 % of all ASVs have genus level
sum(!is.na(tax$Family)) / nrow(tax) * 100 # 96.9 % of all ASVs have family level
sum(!is.na(tax$Phylum)) / nrow(tax) * 100 # 99.9 % of all ASVs have phylum level
nrow(tax)

### Get top 300 ASV by abundance 
ss <- taxa_sums(heliusmb)
ss <- ss[order(ss, decreasing = T)]
ss <- ss[1:300]
top300 <- names(ss)
tax300 <- tax[rownames(tax) %in% top300, ]

sum(!is.na(tax300$Species)) / nrow(tax300) * 100 # 28.7 % of top 300 ASVs have species level
sum(!is.na(tax300$Genus)) / nrow(tax300) * 100 # 84.7 % of top 300 ASVs  have genus level
sum(!is.na(tax300$Family)) / nrow(tax300) * 100 # 97.3 % of top 300 ASVs have family level
sum(!is.na(tax300$Phylum)) / nrow(tax300) * 100 # 100 % of top 300 ASVs have phylum level

### Get taxonomy for ASVs
tax <- tax %>% 
  mutate(Tax = case_when(
    !is.na(Genus) & !is.na(Species) ~ paste(Genus, Species),
    !is.na(Genus) & is.na(Species) ~ paste(Genus, 'spp.'),
    !is.na(Family) & is.na(Genus) ~ paste(Family, 'spp.'),
    !is.na(Order) & is.na(Family) ~ paste(Order, 'spp.'),
    !is.na(Class) & is.na(Order) ~ paste(Class, 'spp.'),
    !is.na(Phylum) & is.na(Class) ~ paste(Phylum, 'spp.'),
    !is.na(Kingdom) & is.na(Phylum) ~ paste(Kingdom, 'spp.'),
    is.na(Kingdom) ~ 'unclassified'),
    ASV = rownames(.)
  )
unique(tax$Tax)    

saveRDS(heliusmb, file = "data/phyloseq.RDS")
saveRDS(tax, file = "data/tax_table.RDS")

## Open HELIUS clinical data
df <- haven::read_sav("data/231109_HELIUS data Barbara Verhaar.sav")

# Change type of variable
df <- df %>% mutate_at(c("H1_lft", "H1_PackYears", "H1_LO_BMI", "H1_LO_WHR", "H1_LO_GemBPSysZit", "H1_LO_GemBPDiaZit", "H1_LO_GemBPHRZit", "H1_Lab_UitslagGLUC", "H1_Lab_UitslagIH1C", "H1_Lab_UitslagCHOL", "H1_Lab_UitslagTRIG","H1_Lab_UitslagHDLS", "H1_Lab_uitslagRLDL"), as.numeric)
df <- df %>% mutate_at(c("H1_LichamelijkOnderzoekJN", "H1_VragenlijstJN", "H1_geslacht", "H1_EtnTotaal", "H1_Roken","H1_AlcoholJN", "H1_AlcoholConsumption", "H1_HT_BPMed","H1_HT_BP", "H1_Antihypertensiva", "H1_Diabetes_GlucMed", "H1_Diabetesmiddelen", "H1_Antilipaemica", "H1_Antibiotica"), as_factor)

# Sex
df$H1_geslacht <- factor(df$H1_geslacht, levels = c("man","vrouw"), labels = c('Male', 'Female'))

# Ethnicity
df$H1_EtnTotaal <- droplevels(df$H1_EtnTotaal)
df$df <- factor(df$H1_EtnTotaal, levels = c("NL","Hind","Creools","Javaans", "Sur anders/onbekend", "Ghanees", "Turks", "Marokkaans", "Anders/onbekend"), labels = c('Dutch','South-Asian Surinamese', 'African Surinamese', 'Other', 'Other', 'Ghanaian', 'Turkish', 'Moroccan', 'Other'))
df$H1_EtnTotaal <- factor(df$H1_EtnTotaal, levels = c("Dutch", "African Surinamese", "South-Asian Surinamese", "Ghanaian", "Turkish", "Moroccan", "Other"))

# Hypertension and BP-lowering medication
df$H1_HT_BPMed <- factor(df$H1_HT_BPMed, levels = c("Nee","Ja"), labels = c('No', 'Yes'))
df$H1_Antihypertensiva <- factor(df$H1_Antihypertensiva, levels = c("Nee","Ja"), labels = c('No', 'Yes'))

# Antibiotics use
df$H1_Antibiotica <- factor(df$H1_Antibiotica, levels = c("Nee","Ja"), labels = c('No', 'Yes'))


# Diabetes mellitus
df$H1_Diabetes_GlucMed <- factor(df$H1_Diabetes_GlucMed, levels = c("Nee","Ja"), labels = c('No', 'Yes'))
df$H1_Diabetesmiddelen <- factor(df$H1_Diabetesmiddelen, levels = c("Nee","Ja"), labels = c('No', 'Yes'))

# Lipid profile and lipid lowering medication
df$H1_Antilipaemica <- factor(df$H1_Antilipaemica, levels = c("Nee","Ja"), labels = c('No', 'Yes'))

# Smoking
df$H1_Roken <- factor(df$H1_Roken, levels = c("Ja", "Nee, ik heb nooit gerookt", "Nee, maar vroeger wel"), labels = c('Yes', 'Never', 'Former'))
df$H1_Roken <- relevel(df$H1_Roken, ref = "Never")  

# Alcohol use
df$H1_AlcoholJN <- factor(df$H1_AlcoholJN, levels = c("ja", "nee"), labels = c('Yes', 'No'))
df$H1_AlcoholConsumption <- factor(df$H1_AlcoholConsumption, levels = c("low (men 0-4 gl/w, women 0-2 gl/w)", "moderate (men 5-14 gl/w, women 3-7 gl/w)", "high (men >14 gl/w, women >7 gl/wk)"), labels = c('Low', 'Moderate', 'High'))

## Select ID's present in both dataframes
sample_names_to_keep <- df$Heliusnr
sample_names_to_keep <- as.character(sample_names_to_keep)
heliusmb <- prune_samples(sample_names_to_keep, heliusmb)

sample_ids <- colnames(otu_table(heliusmb))
df <- df[df$Heliusnr %in% sample_ids, ]

rownames(df) <- df$Heliusnr
df <- df[,-1]

heliusmb@otu_table <- t(heliusmb@otu_table)
p <- merge_phyloseq(heliusmb, df)


## Save files
saveRDS(df, file = "data/clinicaldata.RDS")
saveRDS(p, file = "data/phyloseq_sampledata.RDS")

####################### check missing heliusnrs
heliusdata <- df[,1]
heliusdata$Heliusnr <- as.character(heliusdata$Heliusnr)
otudata <- as.data.frame(colnames(otu))
colnames(otudata)[1] <- "Heliusnr"

join <- anti_join(otudata, heliusdata)
