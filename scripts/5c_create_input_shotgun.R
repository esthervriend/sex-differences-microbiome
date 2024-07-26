# create XGB input files for male/female model

library(dplyr)
library(stringr)
library(ggplot2)
rm(list=ls())

# make data for machine learning XGB classification models

# writes input data files for XGB models as tab-delimited 
# subject ids and feature ids are written as separate tab-delimited files
# write X data / predictors
write_data <- function(x, data_path){
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write y / predicted outcome
write_y <- function(x, name_y, data_path){
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}


## Open dataframe
df <- readRDS('data/clinicaldata.RDS') %>% 
    mutate(Sex = case_when(
        Sex == "Female" ~ 0,
        Sex == "Male" ~ 1
    ))
abundance <- rio::import("data/metaphlan/merged/combined_table.tsv")
abundance2 <- abundance %>% filter(str_detect(clade_name, "s__") & !str_detect(clade_name, "t__"))
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:298], "_T1", ""))
rownames(abundance2) <- abundance2$clade_name

# Short labeling bugs
clade <- abundance2$clade_name
cladesplit <- str_split(clade, "\\|", n = 8, simplify = TRUE)
cladesplit <- as.data.frame(cladesplit[,-8])
colnames(cladesplit) <- c("Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
cladesplit <- cladesplit %>% mutate(across(everything(.), ~str_remove_all(.x, "[a-z]__")))
cladesplit$rowname <- clade
saveRDS(cladesplit, "data/shotgun_taxtable.RDS")
abundance2$clade_name <- NULL

# inspect abundance
tk <- apply(abundance2, 1, function(x) sum(x > 0.1) > (0.25*length(x)))
summary(tk)

# select pathways
abundance2 <- abundance2[tk,]
abundance2 <- t(as.matrix(abundance2))
saveRDS(abundance2, "data/shotgun_abundance.RDS")

# Put clinical data and microbiome data in same sequence of IDs
abundance2 <- abundance2[which(rownames(abundance2) %in% df$ID),]
df <- df[match(rownames(abundance2), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(abundance2)) # TRUE
df$ID
rownames(abundance2)

# make input data
path <- 'sex_metagen'
dir.create(path)
dir.create("sex_metagen/input_data")
write_data(abundance2, file.path(path, 'input_data'))
y <- as.data.frame(df$Sex)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))



## Open dataframe
df <- readRDS('data/clinicaldata.RDS') %>% filter(Sex == "Female")
df$MenopauseYn <- case_when(
                df$MenopauseYn=="Premenopausal" ~ 0,
                df$MenopauseYn=="Postmenopausal" ~ 1)
summary(as.factor(df$MenopauseYn))
df <- df %>% dplyr::select(ID, MenopauseYn)
dim(df)
abundance <- rio::import("data/metaphlan/merged/combined_table.tsv")
abundance2 <- abundance %>% filter(str_detect(clade_name, "s__") & !(str_detect(clade_name, "t__")))
colnames(abundance2) <- c(colnames(abundance2)[1], str_replace(colnames(abundance2)[2:298], "_T1", ""))
rownames(abundance2) <- abundance2$clade_name
abundance2$clade_name <- NULL

# inspect abundance
tk <- apply(abundance2, 1, function(x) sum(x > 0.1) > (0.25*length(x)))
summary(tk)

# select pathways
abundance2 <- abundance2[tk,]
abundance2 <- t(as.matrix(abundance2))

# Put clinical data and microbiome data in same sequence of IDs
abundance2 <- abundance2[which(rownames(abundance2) %in% df$ID),]
df <- df[match(rownames(abundance2), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df$ID == rownames(abundance2)) # TRUE
df$ID
rownames(abundance2)

# make input data
path <- 'menopause_metagen'
dir.create(path)
dir.create("menopause_metagen/input_data")
write_data(abundance2, file.path(path, 'input_data'))
y <- as.data.frame(df$MenopauseYn)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
