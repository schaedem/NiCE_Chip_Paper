require(tidyverse)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")

#at this point, neither NiCE data nor 16S data has been corrected for dilution or standardized

nice_data <- read_csv("NiCE_clean_data.csv")
std_16S_data <- read_csv("processed_16S_qPCR.csv") %>%
  select(-c("timepoint", "block", "location", "treatment"))

head(nice_data)
head(std_16S_data)

no_dup_nice_data <- nice_data %>%
  group_by(sample, assay) %>%
  mutate(Ct = mean(Ct),
         log_conc = mean(log_conc),
         conc = mean(conc)) %>%
  ungroup() %>%
  unique()

#ready to join datasets

full_nice_data <- no_dup_nice_data %>%
  full_join(std_16S_data, by=c("sample")) %>%
  mutate(adj_conc = conc*10,
         adj_conc_16S = conc_16S*10,
         adj_log_conc = log10(adj_conc),
         adj_log_conc_16S = log10(adj_conc_16S))

weird1 <- full_nice_data %>%
  filter(log_conc > log_conc_16S) #1034

weird2 <- full_nice_data %>%
  filter(conc > conc_16S) #1043

weird3 <- full_nice_data %>%
  filter(adj_conc > adj_conc_16S) #1043

weird4 <- full_nice_data %>% #1043
  filter(adj_log_conc > adj_log_conc_16S)

#standardize to copies/g soil
full_adj_df <- full_nice_data %>%
  mutate(
    truesoil = 0.25-(0.25*gwc),
    conc_truesoil = (adj_conc*50)/truesoil, #10x dilution * 50ul elution volume from DNA extraction
    log_conc_truesoil = log10(conc_truesoil),
    conc_16S_truesoil = (adj_conc_16S*50)/truesoil, #unit: copies per gram soil
    log_conc_16S_truesoil = log10(conc_16S_truesoil)
  )

weird5 <- full_adj_df %>% #1043
  filter(log_conc_truesoil > log_conc_16S_truesoil)

weird6 <- full_adj_df %>% #1043
  filter(conc_truesoil > conc_16S_truesoil)

#include 16S archaea in abundance standardization

archaea_data <- full_adj_df %>%
  filter(assay == 2) %>%
  select(sample, adj_conc_16S, adj_conc, adj_log_conc_16S, adj_log_conc, location, timepoint, conc_truesoil, conc_16S_truesoil) %>%
  mutate(adj_conc_16S = as.numeric(adj_conc_16S),
         conc_16S_arch = as.numeric(adj_conc), 
         log_16S_arch = log10(conc_16S_arch),
         bact_arch_abund = adj_conc_16S + conc_16S_arch,
         log_bact_arch_abund = log10(bact_arch_abund),
         conc_total_16S_truesoil = conc_truesoil + conc_16S_truesoil,
         log_total_16S_truesoil = log10(conc_total_16S_truesoil)) %>%
  select(-c(adj_conc, adj_log_conc_16S, adj_conc_16S, conc_truesoil, conc_16S_truesoil, adj_log_conc))


#join archaea abund data back to main df

full_adj_df_final <- full_adj_df %>%
  filter(assay != 2) %>%
  full_join(archaea_data, by = c("sample", "location", "timepoint")) %>%
  mutate(std_conc = adj_conc/bact_arch_abund, 
         log_std_conc = log10(std_conc),  
         std_conc_truesoil = conc_truesoil/conc_total_16S_truesoil,
         log_std_conc_truesoil = log10(std_conc_truesoil))

names(full_adj_df_final)

summary(full_adj_df_final$log_std_conc) #20 NAs - sample 350
summary(full_adj_df_final$log_std_conc_truesoil)
summary(full_adj_df_final$log_conc_truesoil)

see <- full_adj_df_final %>% filter(is.na(log_std_conc))

weird7 <- full_adj_df_final %>%
  filter(log_conc_truesoil > log_total_16S_truesoil) %>% #744
  filter(target_group != "fungi") #527, ~7% of data

weird8 <- full_adj_df_final %>%
  filter(log_std_conc > 0) %>% #744
  filter(target_group != "fungi") #527, 6.9%  of data
  #mutate(adj_log_std_conc = log_bact_arch_abund - log_conc)
  
weird9 <- full_adj_df_final %>%
  filter(log_std_conc_truesoil>0) %>% #744
  filter(target_group != "fungi") #527, 1.2% of data

write.csv(weird9, "assay_samps_above_16S_conc.csv")

#plot histogram to see which standardization looks more normal

#visualize distributions of different standardization methods
ggplot(aes(x=log_std_conc_truesoil, color=location, fill=location), data=full_adj_df_final) +
  geom_histogram(bins=50) +
  #facet_wrap(~timepoint) +
  theme_classic() +
  xlab("") +
  ggtitle("Log (gene copies/g soil / total 16S copies/g soil)")

ggplot(aes(x=log_conc_truesoil, color=location, fill=location), data=full_adj_df_final)+
  #facet_wrap(~timepoint) +
  theme_classic() +
  ggtitle("Log(gene copies/g soil)")

#write final df to .csv
final_write_data <- full_adj_df_final 

write_csv(final_write_data, "final_nice_std_data.csv")

names(final_write_data)

#matrix with final unit: log(copies/g soil)

nice_final_mat <- final_write_data %>%
  select(sample, assay, log_conc_truesoil)

mat_nice_final <- as.matrix(nice_final_mat)

write.table(mat_nice_final, file="nice_matrix.Rdata")
