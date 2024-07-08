# create XGB input files for male/female model

library(dplyr)
library(phyloseq)
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
df <- readRDS('data/clinicaldata.RDS')

## Male vs premenopausal women model
df2 <- df %>% filter(Sex == "Male" | MenopauseYn == "Premenopausal")
dim(df2)
any(is.na(df2$Sex)) # FALSE
df2$Sex <- case_when(
                df2$Sex=="Female" ~ 0,
                df2$Sex=="Male" ~ 1)
summary(as.factor(df2$Sex))
df2 <- df2 %>% dplyr::select(ID, Sex)

## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
mb <- prune_samples(sample_names(mb) %in% df2$ID, mb)
otu <- t(as(mb@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.3*length(x)))
mbdf <- otu[,tk]
dim(mbdf)
mbdf <- as.data.frame(mbdf)

# Put clinical data and microbiome data in same sequence of IDs
df2 <- df2[match(rownames(mbdf), df2$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df2$ID == rownames(mbdf)) # TRUE
df2$ID
rownames(mbdf)

# make input data
path <- 'male_premen'
dir.create(path)
dir.create("male_premen/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df2$Sex)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))


## Male vs postmenopausal women model
df2 <- df %>% filter(Sex == "Male" | MenopauseYn == "Postmenopausal")
dim(df2)
any(is.na(df2$Sex)) # FALSE
df2$Sex <- case_when(
    df2$Sex=="Female" ~ 0,
    df2$Sex=="Male" ~ 1)
summary(as.factor(df2$Sex))
df2 <- df2 %>% dplyr::select(ID, Sex)

## Open RDS file with OTU table
mb <- readRDS('data/phyloseq_sampledata.RDS')
mb <- prune_samples(sample_names(mb) %in% df2$ID, mb)
otu <- t(as(mb@otu_table, "matrix"))
tk <- apply(otu, 2, function(x) sum(x > 5) > (0.3*length(x)))
mbdf <- otu[,tk]
dim(mbdf)
mbdf <- as.data.frame(mbdf)

# Put clinical data and microbiome data in same sequence of IDs
df2 <- df2[match(rownames(mbdf), df2$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(df2$ID == rownames(mbdf)) # TRUE
df2$ID
rownames(mbdf)

# make input data
path <- 'male_postmen'
dir.create(path)
dir.create("male_postmen/input_data")
write_data(mbdf, file.path(path, 'input_data'))
y <- as.data.frame(df2$Sex)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
