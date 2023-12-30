#### Baseline characteristics

## Libraries
library(tableone)
library(rio)

## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")

## Table 1 sex
nonnormal <- c("Age")
Table1 <- CreateTableOne(data = df_new, 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseYn"), 
                         strata = c("Sex"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal)
Table1 <- as.data.frame(Table1)
write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
table1 <- Table1 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(table1, "results/table1.csv")

## Table 2 menopausal status
nonnormal <- c("Age", "MenopauseDuration")
table2 <- CreateTableOne(data = df_new %>% filter(Sex == "Female"), 
                         vars = c("Age", "Ethnicity", "BMI", "CurrSmoking",
                                  "AlcoholIntake", "TotalCalories", "Fibre",
                                  "SBP", "DBP", "TC", "LDL", "Glucose",
                                  "HT", "AntiHT", "DM", "DMMed", "LipidLowering", 
                                  "MenopauseDuration"), 
                         strata = c("MenopauseYn"))
table2 <- print(table2, contDigits = 1, nonnormal = nonnormal)
table2 <- as.data.frame(table2)
write.table(table2, "clipboard", sep="\t", dec=",", col.names=NA)
table2 <- table2 %>% mutate(across(everything(.), ~trimws(.x, which = "both")))
write.csv2(table2, "results/table2.csv")
