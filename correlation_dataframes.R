#Generate dataframes ready for co-occurence network analysis
library(tidyverse)

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/metadata")


nice_mat <- read.table(file="nice_matrix.Rdata")
head(nice_mat)

spread_mat <- spread(nice_mat, key=assay, value=log_conc_truesoil) 

rownames(spread_mat) <- spread_mat$sample

final_mat <- spread_mat[,(2:21)] %>%
  as.matrix()

nice_mat_2 <- read.table(file="nice_matrix.Rdata")
spread_df <- spread(nice_mat, sample, log_conc_truesoil)
abund_df <- spread_df[,-1]
assays <- spread_df[,1]
rownames(abund_df) <- assays

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")

nice <- read_csv("final_nice_std_data.csv") %>%
  select(assay, sample, log_conc_truesoil, primer_F, primer_R) %>%
  mutate(primer_pair = paste(primer_F, primer_R, sep= " / ")) %>%
  select(-c(primer_F, primer_R)) %>%
  mutate(timepoint = case_when(sample %in% 201:256 ~ 1,
                               sample %in% 301:356 ~ 2, 
                               sample %in% 401:456 ~ 3, 
                               sample %in% 501:556 ~ 4,
                               sample %in% 601:656 ~ 5,
                               sample %in% 701:756 ~ 6))

#add filtering by location here?
nice_1 <- filter(nice, timepoint ==1) 
nice_2 <- filter(nice, timepoint ==2)
nice_3 <- filter(nice, timepoint ==3)
nice_4 <- filter(nice, timepoint ==4)
nice_5 <- filter(nice, timepoint ==5)
nice_6 <- filter(nice, timepoint ==6)


#get data into correct format for correlation analysis
df_to_mat <- function(data) {
  data <- data %>%
    select(-c(timepoint, assay))
  
  spread_data <- spread(data, key=primer_pair, value=log_conc_truesoil) %>%
    column_to_rownames('sample')
  
  return(spread_data)
}

nice_1_cor <- df_to_mat(nice_1) 
nice_2_cor <- df_to_mat(nice_2)
nice_3_cor <- df_to_mat(nice_3) 
nice_4_cor <- df_to_mat(nice_4)
nice_5_cor <- df_to_mat(nice_5)
nice_6_cor <- df_to_mat(nice_6) 
nice_all_cor <- df_to_mat(nice)

#log gene quantity matix > spearman correlation matrix

make_cor_mat <- function(dat) {
  cor_dat <- Hmisc::rcorr(as.matrix(dat), type="spearman")
  return(cor_dat)
}

nice_all_mat <- Hmisc::rcorr(final_mat, type="spearman") 
nice_all_cor_mat <- nice_all_mat$r

nice_1_mat <- make_cor_mat(nice_1_cor)
nice_1_cor_mat <- nice_1_mat$r

nice_2_mat <- make_cor_mat(nice_2_cor)
nice_2_cor_mat <- nice_2_mat$r

nice_3_mat <- make_cor_mat(nice_3_cor)
nice_3_cor_mat <- nice_3_mat$r

nice_4_mat <- make_cor_mat(nice_4_cor)
nice_4_cor_mat <- nice_4_mat$r

nice_5_mat <- make_cor_mat(nice_5_cor)
nice_5_cor_mat <- nice_5_mat$r

nice_6_mat <- make_cor_mat(nice_6_cor)
nice_6_cor_mat <- nice_6_mat$r


#For network analysis, we need to change back to a long df format so that each pairwise comparison has an edge value

#we want to have the spearman coefficient and p-value associated with each comparison such that:
#x1  x2  rho  p
# #  #   #   #
# #  #   #   #
# #  #   #   #
# #  #   #   #

make_cor_df <- function(mat) {
  
  cormat <- mat$r
  pmat <- mat$P
  
  ut <- upper.tri(cormat)
  
  df <- data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
  
  df <- df %>% 
    unique() %>%
    rename(a = row,
           b= column,
           r = cor)
  
  return(df)
}

nice_1_cor_df <- make_cor_df(nice_1_mat)
nice_2_cor_df <- make_cor_df(nice_2_mat)
nice_3_cor_df <- make_cor_df(nice_3_mat)
nice_4_cor_df <- make_cor_df(nice_4_mat)
nice_5_cor_df <- make_cor_df(nice_5_mat)
nice_6_cor_df <- make_cor_df(nice_6_mat)
nice_all_cor_df <- make_cor_df(nice_all_mat)


