#Generate dataframes ready for co-occurence network analysis
library(tidyverse)
library(bcdstats)

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/metadata")


nice_mat <- read.table(file="nice_matrix.Rdata")

spread_mat <- spread(nice_mat, key=assay, value=log_conc_truesoil) 

rownames(spread_mat) <- spread_mat$sample

final_mat <- spread_mat[,(2:21)] %>%
  as.matrix()

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")

assays <- read_csv("assay_list_2.csv") %>% select(primer_pair, acronym)

nice <- read_csv("final_nice_std_data.csv") %>%
  select(assay, sample, log_conc_truesoil, primer_F, primer_R, location) %>%
  mutate(primer_pair = paste(primer_F, primer_R, sep= " / ")) %>%
  select(-c(primer_F, primer_R)) %>%
  mutate(timepoint = case_when(sample %in% 201:256 ~ 1,
                               sample %in% 301:356 ~ 2, 
                               sample %in% 401:456 ~ 3, 
                               sample %in% 501:556 ~ 4,
                               sample %in% 601:656 ~ 5,
                               sample %in% 701:756 ~ 6)) %>%
  merge(assays, by="primer_pair")

#loc x time graphs
k_1 <- filter(nice, timepoint ==1 & location=="Karama")
r_1 <- filter(nice, timepoint ==1 & location=="Rubona")

k_2 <- filter(nice, timepoint ==2 & location=="Karama")
r_2 <- filter(nice, timepoint ==2 & location=="Rubona")

k_3 <- filter(nice, timepoint ==3 & location=="Karama")
r_3 <- filter(nice, timepoint ==3 & location=="Rubona")

k_4 <- filter(nice, timepoint ==4 & location=="Karama")
r_4 <- filter(nice, timepoint ==4 & location=="Rubona")

k_5 <- filter(nice, timepoint ==5 & location=="Karama")
r_5 <- filter(nice, timepoint ==5 & location=="Rubona")

k_6 <- filter(nice, timepoint ==6 & location=="Karama")
r_6 <- filter(nice, timepoint ==6 & location=="Rubona")

#get data into correct format for correlation analysis
df_to_mat <- function(data) {
  data <- data %>%
    select(-c(timepoint, assay, location, primer_pair))
  
  spread_data <- spread(data, key=acronym, value=log_conc_truesoil) %>%
    column_to_rownames('sample')
  
  return(spread_data)
}


k_1_cor <- df_to_mat(k_1)
r_1_cor <- df_to_mat(r_1)
k_2_cor <- df_to_mat(k_2)
r_2_cor <- df_to_mat(r_2)
k_3_cor <- df_to_mat(k_3)
r_3_cor <- df_to_mat(r_3)
k_4_cor <- df_to_mat(k_4)
r_4_cor <- df_to_mat(r_4)
k_5_cor <- df_to_mat(k_5)
r_5_cor <- df_to_mat(r_5)
k_6_cor <- df_to_mat(k_6)
r_6_cor <- df_to_mat(r_6)

#log gene quantity per g soil matix > spearman correlation matrix

make_cor_mat <- function(dat) {
  
  cor_dat <- adjust.corr(dat, type="spearman", adjust="BH")
  
  #cor_dat <- Hmisc::rcorr(as.matrix(dat), type="spearman")
  return(cor_dat)
}

#use adjust.corr() from bcdstats package to get BH - adjusted p-values
#necessary to adjust for error introduced from multiple comparisons

k_1_mat <- make_cor_mat(k_1_cor)
r_1_mat <- make_cor_mat(r_1_cor)

k_2_mat <- make_cor_mat(k_2_cor)
r_2_mat <- make_cor_mat(r_2_cor)

k_3_mat <- make_cor_mat(k_3_cor)
r_3_mat <- make_cor_mat(r_3_cor)

k_4_mat <- make_cor_mat(k_4_cor)
r_4_mat <- make_cor_mat(r_4_cor)

k_5_mat <- make_cor_mat(k_5_cor)
r_5_mat <- make_cor_mat(r_5_cor)

k_6_mat <- make_cor_mat(k_6_cor)
r_6_mat <- make_cor_mat(r_6_cor)

#For network analysis, we need to change back to a long df format so that each pairwise comparison has an edge value

#we want to have the spearman coefficient and p-value associated with each comparison such that:
#x1  x2  rho  p
# #  #   #   #
# #  #   #   #
# #  #   #   #
# #  #   #   #

make_cor_df <- function(mat) {
  
  cormat <- mat[["R"]][["r"]]
  pmat <- mat[["P"]]
  
  ut <- upper.tri(cormat)
  
  df <- data.frame(
    row = rownames(cormat)[row(cormat)[ut]],
    column = rownames(cormat)[col(cormat)[ut]],
    cor  =(cormat)[ut],
    p = pmat[ut]
  )
  
  df <- df %>% 
    unique() %>%
    dplyr::rename(a = row,
           b= column,
           r = cor) %>%
    mutate(p = as.numeric(substr(p, 2,6)))
  
  return(df)
}

k_1_cor_df <- make_cor_df(k_1_mat)
r_1_cor_df <- make_cor_df(r_1_mat)

k_2_cor_df <- make_cor_df(k_2_mat)
r_2_cor_df <- make_cor_df(r_2_mat)

k_3_cor_df <- make_cor_df(k_3_mat)
r_3_cor_df <- make_cor_df(r_3_mat)

k_4_cor_df <- make_cor_df(k_4_mat)
r_4_cor_df <- make_cor_df(r_4_mat)

k_5_cor_df <- make_cor_df(k_5_mat)
r_5_cor_df <- make_cor_df(r_5_mat)

k_6_cor_df <- make_cor_df(k_6_mat)
r_6_cor_df <- make_cor_df(r_6_mat)
