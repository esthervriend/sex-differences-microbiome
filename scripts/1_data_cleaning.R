#### Data cleaning

## Libraries
library(tidyverse)
library(dplyr)
library(rio) 
library(haven)
library(phyloseq)
library(forcats)
library(stringr)

## Open HELIUS clinical data
df <- haven::read_sav("data/231109_HELIUS data Barbara Verhaar_240702.sav")
df$H1_LftMenstrWeg <- as.numeric(df$H1_LftMenstrWeg)

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
                Sodium=Natrium_intake_totaal_gram,
                # Medication
                DMMed=H1_Diabetesmiddelen, AntiHT=H1_Antihypertensiva, 
                AB=H1_Antibiotica, LipidLowering=H1_Antilipaemica, Antihistaminics=H1_Antihistaminica,
                Corticosteroids=H1_Corticosteroiden,SystSteroids=H1_SystSteroiden,
                Antithromb=H1_Antithrombotica, Anticon=H1_Anticon2, HRT=H1_HRT_type_Barbara, PPI=H1_ProtPumpInh_Barbara,
                # Antihypertensive medication
                AllAntiHT=H1_Antihypertensiva, Diuretics=H1_AntihypertensivaC03,
                CalciumAnt=H1_AntihypertensivaC08, BetaBlocker=H1_AntihypertensivaC07,
                RAASi=H1_AntihypertensivaC09,
                # Fecal sample 
                FecalSample_AB=H1_Feces_q2, FecalSample_Diarrhoea=H1_Feces_q3
  ) %>% 
  mutate(across(where(is.character), ~na_if(., c("Missing", "niet ingevuld","nvt"))),
         across(c("Physical", "Questionnaire", "Sex", "Ethnicity", "Smoking",
                  "AlcCons","Antihistaminics","SystSteroids", "HRT", "Anticon"), as_factor),
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
         HRT = fct_recode(HRT, "No"= "Did not have hormone replacement therapy",
                          "Yes" = "Oestrogen only (ATC-codes G03C)",
                          "Yes" = "Oestrogen and progestogen combination (ATC-codes G03F or G03C&G03D or G03HB01)",
                          "Yes" = "Tibolone (ATC-codes G03CX01 or G03DC05)",
                          "Yes" = "Other types of HRT (incl. progestogen only, vaginal and other local treatments, and combinations of other subtypes)"),
         Anticon = fct_recode(Anticon, "No"="Niet geselecteerd",
                              "Yes"="Ja"),
         Alcohol = yesnosmall(Alcohol), # coded without caps
         CurrSmoking = fct_recode(Smoking, "Yes" = "Ja", 
                                  "No" = "Nee, ik heb nooit gerookt", 
                                  "No" = "Nee, maar vroeger wel"),
         MenopauseAge = case_when(
           MenopauseAge > Age ~ NA_real_,
           TRUE ~ MenopauseAge),
         CurrSmoking = fct_infreq(CurrSmoking),
         AlcCons = factor(AlcCons, levels = c("low (men 0-4 gl/w, women 0-2 gl/w)",
                                              "moderate (men 5-14 gl/w, women 3-7 gl/w)",
                                              "high (men >14 gl/w, women >7 gl/wk)"), 
                          labels = c('Low', 'Moderate', 'High')),
         across(where(is.numeric), as.numeric), # all other vars to numeric, do this last,
         MenopauseYn = case_when(
           is.na(MenopauseAge) & Sex == "Female" ~ "Premenopausal",
           !is.na(MenopauseAge) ~ "Postmenopausal",
           TRUE ~ NA_character_),
         MenopauseYn = fct_rev(as.factor(MenopauseYn)),
         MenopauseDuration = case_when(
           MenopauseYn == "Postmenopausal" ~ Age - MenopauseAge,
           TRUE ~ NA_real_
         ),
         ID = str_c("S", ID)
  ) %>% 
  # remove unused levels
  droplevels(.)

df_new <- df_new %>% 
  # filter out all participants that used antibiotics at the moment of sample collection
  filter(FecalSample_AB == "No")
dim(df_new)

# Open phyloseq object
heliusmb <- readRDS("data/phyloseq_rarefied.RDS")
sample_names(heliusmb) <- str_c("S", sample_names(heliusmb))

# Select samples that are in dataset
heliusmb <- prune_samples(sample_names(heliusmb) %in% df_new$ID, heliusmb)
heliusmb
df_new <- df_new %>% filter(ID %in% sample_names(heliusmb))
all(df_new$ID %in% sample_names(heliusmb)) # TRUE

## Save files
saveRDS(df_new, file = "data/clinicaldata.RDS")
saveRDS(heliusmb, file = "data/phyloseq_sampledata.RDS")
