#### Baseline characteristics

## Libraries
library(tableone)
library(rio)
library(dplyr)

## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")

## Table 1 sex
nonnormal <- c("Age")
Table1 <- CreateTableOne(data = df_new, 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseYn", "HRT", "Anticon"), 
                         strata = c("Sex"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal)
Table1 <- as.data.frame(Table1)
write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
Table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table1, "results/Table1.csv")

## Table 2 menopausal status
nonnormal <- c("Age", "MenopauseDuration")
Table2 <- CreateTableOne(data = df_new %>% filter(Sex == "Female"), 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseDuration", "HRT", "Anticon"), 
                         strata = c("MenopauseYn"))
Table2 <- print(Table2, contDigits = 1, nonnormal = nonnormal)
Table2 <- as.data.frame(Table2)
write.table(Table2, "clipboard", sep="\t", dec=",", col.names=NA)
Table2 <- Table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table2, "results/Table2.csv")



## Load dataset
df_new <- rio:: import("data/clinicaldata_shotgun.RDS")

## Table 1 sex
nonnormal <- c("Age")
Table1 <- CreateTableOne(data = df_new, 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseYn", "HRT", "Anticon"), 
                         strata = c("Sex"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal)
Table1 <- as.data.frame(Table1)
#write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
Table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table1, "results/Table1_shotgun.csv")

## Table 2 menopausal status
nonnormal <- c("Age", "MenopauseDuration")
Table2 <- CreateTableOne(data = df_new %>% filter(Sex == "Female"), 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseDuration", "HRT", "Anticon"), 
                         strata = c("MenopauseYn"))
Table2 <- print(Table2, contDigits = 1, nonnormal = nonnormal)
Table2 <- as.data.frame(Table2)
#write.table(Table2, "clipboard", sep="\t", dec=",", col.names=NA)
Table2 <- Table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(Table2, "results/Table2_shotgun.csv")
