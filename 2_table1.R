#### Baseline characteristics

## Libraries
library(tableone)
library(rio)


## Load dataset
df_new <- rio:: import("data/clinicaldata.RDS")


## Create tableone
nonnormal <- c("Lft", "MenopauseDuration")
Table1 <- CreateTableOne(data = df_new, vars = c("Age", "Ethnicity", "BMI", "HT", "AntiHT", "Smoking", "DM", "DMMed", "LipidLowering", "AB", "MenopauseYn", "MenopauseDuration"), strata = c("Sex"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal)
Table1 <- as.data.frame(Table1)
write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
