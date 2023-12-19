#### Data cleaning

## Libraries
library(tidyverse)
library(rio) # kun je alles mee openen zonder te definieren wat voor bestand het is
library(haven)
library(phyloseq) # voor handling van phyloseq object

## Open phyloseq object
otu <- rio::import("data/HELUIS.Verhaar.2023.otu_table.csv", header = TRUE)
rownames(otu) <- otu[,1]
otu <- otu[,-1]
tax <- rio::import("Data/HELUIS.Verhaar.2023.tax_table.csv", header = TRUE)
rownames(tax) <- tax$V1
tax$V1 <- NULL
tax <- as.matrix(tax)
tree <- ape::read.tree("data/HELUIS.Verhaar.2023.phy_tree.tree")
# sampledata <- rio::import("data/HELUIS.Verhaar.2023.sam_data.csv", header = TRUE) = empty

## Make phyloseq
heliusmb <- phyloseq(otu_table(otu, taxa_are_rows = T), tax_table(tax), phy_tree(tree))
sample_names(heliusmb) <-str_c("S",sample_names(heliusmb)) # safeguard that sample names cannot be used as numeric in matrices
sample_names(heliusmb)

## Inspect phyloseq object
ntaxa(heliusmb)
nsamples(heliusmb)
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
yesnosmall <- function(x) xnew = fct_recode(as_factor(x, levels=c("labels")), "No"="nee", "Yes"="ja")
yesnocaps <- function(x) xnew = fct_recode(as_factor(x, levels=c("labels")), "No"="Nee", "Yes"="Ja")

## Clean HELIUS dataframe
df_new <- df %>% 
    dplyr::select(ID=Heliusnr, Age=H1_lft, Sex=H1_geslacht, Ethnicity=H1_EtnTotaal, 
           Smoking=H1_Roken, PackYears=H1_PackYears, Alcohol=H1_AlcoholJN, 
           AlcCons=H1_AlcoholConsumption, Questionnaire=H1_VragenlijstJN, 
           MenopauseAge=H1_LftMenstrWeg,
           # Physical exam
           Physical=H1_LichamelijkOnderzoekJN,
           BMI=H1_LO_BMI, WHR=H1_LO_WHR, SBP=H1_LO_GemBPSysZit, 
           DBP=H1_LO_GemBPDiaZit, HR=H1_LO_GemBPHRZit,
           # Cardiometabolic
           HT=H1_HT_BPMed, HT_SelfBP=H1_HT_SelfBP, 
           MetSyn=H1_MetSyn_MetabolicSyndrome, DM=H1_Diabetes_GlucMed,
           FramRisk=H1_Fram_CVD, SCORENL=H1_SCORE_CVDmort_NL,
           # Lab
           TC = H1_Lab_UitslagCHOL, LDL=H1_Lab_uitslagRLDL, HDL=H1_Lab_UitslagHDLS, 
           Trig=H1_Lab_UitslagTRIG, HbA1c=H1_Lab_UitslagIH1C, Glucose = H1_Lab_UitslagGLUC,
           # Diet
           TotalCalories=ENKcal_Sum, Proteins=Prot_Sum, FattyAcids=FattyAcidsTot_Sum,
           SatFat=SFA_Sum, MonoUnsatFat=MUFA_cis_Sum, PolyUnsatFat=PUFA_Sum,
           Fibre=Fibre_Sum, Carbohydrates=Carbo_Sum, AlcoholIntake=Alcohol_Sum,
           # Medication
           DMMed=H1_Diabetesmiddelen, AntiHT=H1_Antihypertensiva, 
           AB=H1_Antibiotica, LipidLowering=H1_Antilipaemica, Antihistaminics=H1_Antihistaminica,
           Corticosteroids=H1_Corticosteroiden,SystSteroids=H1_SystSteroiden,
           Antithromb=H1_Antithrombotica,
           # Antihypertensive medication
           AllAntiHT=H1_Antihypertensiva, Diuretics=H1_AntihypertensivaC03,
           CalciumAnt=H1_AntihypertensivaC08, BetaBlocker=H1_AntihypertensivaC07,
           RAASi=H1_AntihypertensivaC09,
           # Fecal sample 
           FecalSample_AB=H1_Feces_q2, FecalSample_Diarrhoea=H1_Feces_q3
           ) %>% 
            mutate(across(where(is.character), ~na_if(., c("Missing", "niet ingevuld","nvt"))),
                    across(c("Physical", "Questionnaire", "Sex", "Ethnicity", "Smoking",
                             "AlcCons","Antihistaminics","SystSteroids"), as_factor),
                    across(c("FecalSample_AB","FecalSample_Diarrhoea"),
                           ~as_factor(.x, levels = c("labels"))),
                   across(c("HT_SelfBP", "HT", "DM", "AntiHT", "AllAntiHT", "Diuretics",
                            "CalciumAnt", "BetaBlocker", "RAASi", "Antithromb", 
                            "Corticosteroids",  "LipidLowering","DMMed","AB"), yesnocaps)
                   ) %>% 
            droplevels(.) %>% 
            mutate(Sex = factor(Sex, levels = c("man","vrouw"), 
                                 labels = c('Male', 'Female')),
                   Ethnicity = forcats::fct_recode(Ethnicity, "Dutch"="NL",
                                                   "South-Asian Surinamese"="Hind",
                                                   "African Surinamese"="Creools",
                                                   "Ghanaian"="Ghanees",
                                                   "Turkish"="Turks",
                                                   "Morroccan"="Marokkaans",
                                                   "Other"="Anders/onbekend",
                                                   "Other"="Sur anders/onbekend",
                                                   "Other"="Javaans"),
                   Alcohol = yesnosmall(Alcohol), # coded without caps
                   CurrSmoking = fct_recode(Smoking, "Yes" = "Ja", 
                                            "No" = "Nee, ik heb nooit gerookt", 
                                            "No" = "Nee, maar vroeger wel"),
                   MenopauseAge = case_when(MenopauseAge > Age ~ NA, .default = MenopauseAge),
                   CurrSmoking = fct_infreq(CurrSmoking),
                   AlcCons = factor(AlcCons, levels = c("low (men 0-4 gl/w, women 0-2 gl/w)",
                                                        "moderate (men 5-14 gl/w, women 3-7 gl/w)",
                                                        "high (men >14 gl/w, women >7 gl/wk)"), 
                                        labels = c('Low', 'Moderate', 'High')),
                   across(where(is.numeric), as.numeric), # all other vars to numeric, do this last,
                   MenopauseYn = case_when(
                       is.na(MenopauseAge) & Sex == "Female" ~ "No",
                       MenopauseAge == "Yes" ~ "Yes",
                       .default = NA
                   ),
                   MenopauseDuration = case_when(
                       MenopauseYn == "Yes" ~ Age - MenopauseAge,
                       .default = NA
                   ),
                   ID = str_c("S", ID)
            ) %>% 
            # remove unused levels
            droplevels(.) %>% 
            # filter out all participants that used antibiotics at the moment of sample collection
            filter(FecalSample_AB == "No") 
dim(df_new)

## Select IDs present in both dataframes
sample_names_to_keep <- df_new$ID
sample_names_to_keep <- as.character(sample_names_to_keep)
heliusmb <- prune_samples(sample_names_to_keep, heliusmb)

sample_ids <- colnames(otu_table(heliusmb))
df_new <- df_new[df_new$ID %in% sample_ids, ]
dim(df_new)
rownames(df_new) <- df_new$ID
df_new <- df_new[,-1]

heliusmb@otu_table <- t(heliusmb@otu_table)
heliusmbcomplete <- merge_phyloseq(heliusmb, df_new)

## Save files
saveRDS(df_new, file = "data/clinicaldata.RDS")
saveRDS(heliusmbcomplete, file = "data/phyloseq_sampledata.RDS")

# ####################### check missing heliusnrs
# heliusdata <- df_new[,1]
# heliusdata$ID <- as.character(heliusdata$ID)
# otudata <- as.data.frame(colnames(otu))
# colnames(otudata)[1] <- "ID"
# 
# join <- anti_join(otudata, heliusdata)
