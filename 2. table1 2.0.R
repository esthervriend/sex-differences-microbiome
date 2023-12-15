#### Baseline characteristics

## Libraries
library(tableone)
library(rio)

## Load dataset
df <- rio:: import("data/clinicaldata.RDS")


## Create tableone
nonnormal <- c("H1_lft")
Table1 <- CreateTableOne(data = df, vars = c("H1_lft", "H1_EtnTotaal", "H1_LO_BMI", "H1_HT_BPMed", "H1_Antihypertensiva", "H1_Roken", "H1_Diabetes_GlucMed", "H1_Diabetesmiddelen", "H1_Lab_UitslagCHOL","H1_Antilipaemica", "H1_Antibiotica"), strata = c("H1_geslacht"))
Table1 <- print(Table1, contDigits = 1, nonnormal = nonnormal)
Table1 <- as.data.frame(Table1)
write.table(Table1, "clipboard", sep="\t", dec=",", col.names=NA)
