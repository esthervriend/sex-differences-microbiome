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
pathways <- rio::import("data/humann/merged/pathway_abundance_cpm_unstratified.txt")
pathways2 <- pathways %>% select(1, !(contains("Abundance-RPKs") | contains("Coverage")))
colnames(pathways2) <- c(colnames(pathways2)[1], str_replace(colnames(pathways2)[2:298], "_T1_Abundance", ""))

# Short labeling pathways
genes <- pathways2$`# Gene Family`
geneskey <- str_split(genes, ": ", n = 2, simplify = TRUE)
colnames(geneskey) <- c("key","label")
geneskey <- as.data.frame(geneskey)
print('Check if any keys are duplicated:')
any(duplicated(geneskey$key)) # unique keys
rownames(pathways2) <- geneskey$key

# Remove zero abundance pathways
print('Remove zero abundance pathways..')
names(pathways2)
pathways2 <- pathways2 %>% filter(!`# Gene Family` %in% c("UNMAPPED", "UNINTEGRATED"))
pathways2[,1] <- NULL

# inspect abundance
print('Select pathways with abundance >50 cpm')
sums <- as.data.frame(rowSums(pathways2))
colnames(sums) <- "prev"
sums <- sums %>% mutate(prev = prev / 297)
sumsprev <- sums %>% filter(prev > 50) %>% print()
pathwaysprev <- rownames(sumsprev)

# select pathways
print('Select pathways..')
pathways2 <- pathways2 %>% filter(rownames(.) %in% pathwaysprev)
pathways2 <- t(as.matrix(pathways2))
saveRDS(pathways2, "data/pathways_filtered.RDS")

# Put clinical data and microbiome data in same sequence of IDs
print('Put clinical data in same sequence as pathways..')
pathways2 <- pathways2[which(rownames(pathways2) %in% df$ID),]
dim(pathways2)
df <- df[match(rownames(pathways2), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
print('Check IDs..')
all(df$ID == rownames(pathways2)) # TRUE
df$ID
rownames(pathways2)

# make input data
path <- 'sex_pathways'
dir.create(path)
dir.create("sex_pathways/input_data")
write_data(pathways2, file.path(path, 'input_data'))
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
pathways <- rio::import("data/humann/merged/pathway_abundance_cpm_unstratified.txt")
pathways2 <- pathways %>% select(1, !(contains("Abundance-RPKs") | contains("Coverage")))
colnames(pathways2) <- c(colnames(pathways2)[1], str_replace(colnames(pathways2)[2:298], "_T1_Abundance", ""))

# Short labeling pathways
genes <- pathways2$`# Gene Family`
geneskey <- str_split(genes, ": ", n = 2, simplify = TRUE)
colnames(geneskey) <- c("key","label")
geneskey <- as.data.frame(geneskey)
any(duplicated(geneskey$key)) # unique keys
rownames(pathways2) <- geneskey$key
pathways2[,1] <- NULL

# Remove zero abundance pathways
print('Remove zero abundance pathways..')
pathways2 <- pathways2 %>% filter(!rownames(.) %in% c("UNMAPPED", "UNINTEGRATED")) 
# filter(rowSums(.) != 0)

# inspect abundance
print('Select pathways with abundance >50 cpm')
sums <- as.data.frame(rowSums(pathways2))
colnames(sums) <- "prev"
sums <- sums %>% mutate(prev = prev / 297)
sumsprev <- sums %>% filter(prev > 50) %>% print()
pathwaysprev <- rownames(sumsprev)

# select pathways
print('Select pathways..')
pathways2 <- pathways2 %>% filter(rownames(.) %in% pathwaysprev)
pathways2 <- t(as.matrix(pathways2))

# Put clinical data and microbiome data in same sequence of IDs
print('Put clinical data in same sequence as pathways..')
pathways2 <- pathways2[which(rownames(pathways2) %in% df$ID),]
df <- df[match(rownames(pathways2), df$ID), ]

# check that outcome subject ids match metabolite subjects ids
print('Check IDs..')
all(df$ID == rownames(pathways2)) # TRUE
df$ID
rownames(pathways2)

# make input data
path <- 'menopause_pathways'
dir.create(path)
dir.create("menopause_pathways/input_data")
write_data(pathways2, file.path(path, 'input_data'))
y <- as.data.frame(df$MenopauseYn)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
