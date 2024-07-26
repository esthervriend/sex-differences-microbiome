#### Baseline characteristics

## Libraries
library(tableone)
library(rio)
library(dplyr)

## Output folder
resultsfolder <- "results/tables"
dir.create(resultsfolder, showWarnings = FALSE)

## 16S data
## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")

## Table 1 sex
nonnormal <- c("Age", "Sodium", "AlcoholIntake", "TotalCalories")
Table1 <- CreateTableOne(data = df_new, 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "SBP", "DBP", "HT", "AntiHT", "TC", "LDL", "LipidLowering",
                                  "Glucose", "DM", "DMMed",
                                  "Sodium", "AlcoholIntake", "TotalCalories", "Fibre",
                                  "MenopauseYn", "HRT", "Anticon"), 
                         strata = c("Sex"),
                         addOverall = T)
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal, missing = TRUE)
Table1 <- as.data.frame(Table1)
write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
Table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table1, "results/tables/Table1.csv")

## Supplementary table 1 menopausal status
nonnormal <- c("Age", "MenopauseDuration", "Sodium", "AlcoholIntake", "TotalCalories")
Table2 <- CreateTableOne(data = df_new %>% filter(Sex == "Female"), 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "SBP", "DBP", "HT", "AntiHT", "TC", "LDL", "LipidLowering",
                                  "Glucose", "DM", "DMMed",
                                  "Sodium", "AlcoholIntake", "TotalCalories", "Fibre",
                                  "MenopauseDuration", "HRT", "Anticon"), 
                         strata = c("MenopauseYn"))
Table2 <- print(Table2, contDigits = 1, nonnormal = nonnormal, missing = TRUE)
Table2 <- as.data.frame(Table2)
write.table(Table2, "clipboard", sep="\t", dec=",", col.names=NA)
Table2 <- Table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table2, "results/tables/Table2.csv")


## Shotgun data ####
## Load dataset
df <- readRDS("data/clinicaldata.RDS")
tab <- readRDS("data/shotgun_abundance.RDS")
tab <- as.data.frame(tab)
tab$ID <- rownames(tab)
df_shotgun <- merge(df, tab, by = "ID")


## Table 1 sex
df_shotgun$Ethnicity <- droplevels(df_shotgun$Ethnicity )
nonnormal <- c("Age")
Table1 <- CreateTableOne(data = df_shotgun, 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking", "SBP", "DBP", "HT", 
                                  "AntiHT", "TC", "LDL", "LipidLowering", "Glucose", "DM", "DMMed",
                                  "Sodium", "AlcoholIntake", "TotalCalories",
                                  "Fibre"),
                         strata = c("Sex"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal, missing = TRUE)
Table1 <- as.data.frame(Table1)
#write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
Table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table1, "results/tables/Table1_shotgun.csv")

## Table 2 menopausal status
nonnormal <- c("Age", "MenopauseDuration")
Table2 <- CreateTableOne(data = df_shotgun %>% filter(Sex == "Female"), 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking", "SBP", "DBP", "HT", 
                                  "AntiHT", "TC", "LDL", "LipidLowering", "Glucose", "DM", "DMMed",
                                  "Sodium", "AlcoholIntake", "TotalCalories",
                                  "Fibre",
                                  "MenopauseDuration", "HRT", "Anticon"), 
                         strata = c("MenopauseYn"))
Table2 <- print(Table2, contDigits = 1, nonnormal = nonnormal, missing = TRUE)
Table2 <- as.data.frame(Table2)
#write.table(Table2, "clipboard", sep="\t", dec=",", col.names=NA)
Table2 <- Table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table2, "results/tables/Table2_shotgun.csv")
