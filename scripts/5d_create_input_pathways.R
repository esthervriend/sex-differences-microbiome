# create XGB input files for pred eth/timepoint from pathways

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
pathways <- rio::import("data/humann/merged_pathway_renorm_unstratified.tsv")
rownames(pathways) <- pathways$`# Pathway`
pathways$`# Pathway` <- NULL
pathways <- as.data.frame(t(as.matrix(pathways)))
rownames(pathways) <- str_remove(rownames(pathways), "_T1_Abundance-CPM")

# Short labeling pathways
paths <- colnames(pathways)
pathskey <- str_split(paths, ": ", n = 2, simplify = TRUE)
colnames(pathskey) <- c("key","label")
pathskey <- as.data.frame(pathskey)
print('Check if any keys are duplicated:')
any(duplicated(pathskey$key)) # FALSE = unique keys
colnames(pathways) <- pathskey$key
saveRDS(pathskey, "data/pathway_keys.RDS")

# inspect abundance
print('Select pathways with abundance >50 cpm')
sums <- as.data.frame(colSums(pathways))
colnames(sums) <- "prev"
sums <- sums %>% mutate(prev = prev / 297)
hist(sums$prev)
tk <- apply(pathways, 2, function(x) sum(x > 200) > (0.3*length(x)))
summary(tk)
# sumsprev <- sums %>% filter(prev > 50) %>% print()
# pathwaysprev <- rownames(sumsprev)
print('Select pathways..')
pathways2 <- pathways[,tk] # select pathways
rownames(pathways2) <- rownames(pathways)
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
dir.create(path, showWarnings = FALSE)
dir.create("sex_pathways/input_data", showWarnings = FALSE)
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
pathways <- rio::import("data/humann/merged_pathway_renorm_unstratified.tsv")
rownames(pathways) <- pathways$`# Pathway`
pathways$`# Pathway` <- NULL
pathways <- as.data.frame(t(as.matrix(pathways)))
rownames(pathways) <- str_remove(rownames(pathways), "_T1_Abundance-CPM")

# Short labeling pathways
paths <- colnames(pathways)
pathskey <- str_split(paths, ": ", n = 2, simplify = TRUE)
colnames(pathskey) <- c("key","label")
pathskey <- as.data.frame(pathskey)
print('Check if any keys are duplicated:')
any(duplicated(pathskey$key)) # unique keys
colnames(pathways) <- pathskey$key

# inspect abundance
print('Select pathways with abundance >50 cpm')
sums <- as.data.frame(colSums(pathways))
colnames(sums) <- "prev"
sums <- sums %>% mutate(prev = prev / 297)
hist(sums$prev)
tk <- apply(pathways, 2, function(x) sum(x > 200) > (0.3*length(x)))
summary(tk)
# sumsprev <- sums %>% filter(prev > 50) %>% print()
# pathwaysprev <- rownames(sumsprev)
print('Select pathways..')
pathways2 <- pathways[,tk] # select pathways
rownames(pathways2) <- rownames(pathways)

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
path <- 'menopause_pathways'
dir.create(path, showWarnings = FALSE)
dir.create("menopause_pathways/input_data", showWarnings = FALSE)
write_data(pathways2, file.path(path, 'input_data'))
y <- as.data.frame(df$MenopauseYn)
y
write_y(y, name_y = 'y_binary.txt', file.path(path, 'input_data'))
