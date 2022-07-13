library(tidyverse)
library(vegan)
library(reshape2)
library(bioDist)
#source("cooccurrence_permanova.R")

setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/metadata")

nice_mat <- read.table(file="nice_matrix.Rdata") 
meta <- read_csv("SNC_metadata_2020.csv") %>%
  mutate(growth = case_when(timepoint==1 ~ "anthesis",
                            timepoint ==2 ~ "regrowth",
                            timepoint == 3 ~ "anthesis",
                            timepoint ==4 ~ "regrowth",
                            timepoint == 5 ~ "anthesis",
                            timepoint == 6 ~ "regrowth"))

nice_dist_mat <- nice_mat %>%
  arrange(sample) %>%
  pivot_wider(names_from=assay, values_from=log_conc_truesoil) %>%
  as.data.frame()

rownames(nice_dist_mat) <- nice_dist_mat$sample
nice_dist_mat <- nice_dist_mat[,2:21]

nice_dist_mat <- as.matrix(nice_dist_mat)

dist <- vegdist(nice_dist_mat, method="bray")

#permanova testing
adonis_df <- as.matrix(dist) %>%
  as.data.frame()

adonis_merge <- rownames_to_column(adonis_df, var="sample") %>% 
  as_tibble() %>%
  merge(meta, by="sample")

set.seed(123)
#timepoint_adonis <- adonis2(dist ~ timepoint, data=adonis_merge, permutations=9999)
#timepoint highly sig, R2=0.48792, p=0.0001

set.seed(123)
#timepoint_loc_adonis <- adonis2(dist ~ timepoint*location*block*treatment, 
  #                              data=adonis_merge, permutations=9999)
#timepoint_loc_adonis
#timepoint, location, and timepoint:location are sig

set.seed(123)
#growth_all_adonis <- adonis2(dist~as.factor(timepoint)*location*gwc*growth, 
 #                            data=adonis_merge, permutations=9999)
#growth_all_adonis

####
#pairwise comparisons to find sig differences between timepoints
# pairwise_df <- numeric() #initialize empty df to put p-values
# 
# pairwise_adonis <- function(val1, val2, adonis_merge) {
#   val1 <- val1
#   val2 <- val2
#   
#   filter_merge <- adonis_merge %>%
#     filter(timepoint == val1 | timepoint == val2)
#   
#   dist_filter <- filter_merge %>%
#     select(all_of(.[["sample"]])) %>%
#     as.dist()
#   
#   set.seed(123)
#   result <- adonis(dist_filter ~ timepoint, data=filter_merge)
#   
#   result_p <- result[["aov.tab"]][["Pr(>F)"]][1]
#   
#   return(result_p)
#   
# }
# 
# # pairwise_1_2 <- pairwise_adonis(1,2, adonis_merge)
# # pairwise_1_3 <- pairwise_adonis(1,3, adonis_merge) 
# # pairwise_1_4 <- pairwise_adonis(1,4, adonis_merge) 
# # pairwise_1_5 <- pairwise_adonis(1,5, adonis_merge) 
# # pairwise_1_6 <- pairwise_adonis(1,6, adonis_merge) 
# # pairwise_2_3 <- pairwise_adonis(2,3, adonis_merge) 
# # pairwise_2_4 <- pairwise_adonis(2,4,adonis_merge)
# #pairwise_2_5 <- pairwise_adonis(2,5, adonis_merge)
# #pairwise_2_6 <- pairwise_adonis(2,6, adonis_merge)
# #pairwise_3_4 <- pairwise_adonis(3,4, adonis_merge)
# #pairwise_3_5 <- pairwise_adonis(3,5, adonis_merge)
# #pairwise_3_6 <- pairwise_adonis(3,6, adonis_merge)
# #pairwise_4_5 <- pairwise_adonis(4,5, adonis_merge)
# #pairwise_4_6 <- pairwise_adonis(4,6,adonis_merge)
# #pairwise_5_6 <- pairwise_adonis(5,6,adonis_merge)
# 
# #pairwise_df[c("1-2", "1-3", "1-4", "1-5", "1-6", "2-3", "2-4", "2-5", "2-6",
#               "3-4", "3-5", "3-6", "4-5", "4-6", "5-6")] <- 
# #   c(pairwise_1_2, pairwise_1_3, pairwise_1_4,pairwise_1_5,pairwise_1_6,
# #     pairwise_2_3, pairwise_2_4 ,pairwise_2_5 ,pairwise_2_6 ,pairwise_3_4 ,pairwise_3_5 ,
# #     pairwise_3_6 ,pairwise_4_5 , pairwise_4_6 ,pairwise_5_6  )
# # 
# # pairwise_p <- p.adjust(pairwise_df, method="BH") #adjust for multiple comparisons
# # all(pairwise_p < 0.05)
# 
# ####
# set.seed(222)
# loc_adonis <- adonis2(dist ~ location, data=adonis_merge)
# #location sig, R2=0.028, p=0.001
# set.seed(222)
# gwc_adonis <- adonis2(dist~gwc, data=adonis_merge)
# #gwc sig, R2=0.14665, p=0.001
# set.seed(222)
# treat_adonis <- adonis2(dist~treatment, data=adonis_merge)
# #treatment not sig, R2=0.0047, P=0.993
# set.seed(222)
# growth_adonis <- adonis2(dist~growth, data=adonis_merge)
# #growth sig, R2=0.03028, P=0.002; likely conflated with other env variables
# 
# #NMDS plots
# set.seed(100192)
# nmds <- metaMDS(dist)
# 
# nmds_plot <- scores(nmds) %>%
#   as_tibble(rownames="sample") %>%
#   merge(meta, by="sample")
# 
# abund_nmds_growth <- 
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=growth)) + 
#   geom_point(alpha=0.8) + 
#   theme_classic() +
#   scale_color_brewer(palette="Set2") +
#   theme(legend.title=element_blank())
# abund_nmds_growth
# 
# abund_nmds_time <- 
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) + 
#   geom_point() + 
#   theme_classic() +
#   viridis::scale_color_viridis(discrete="TRUE", name="Timepoint")
# abund_nmds_time
# 
# abund_nmds_plant <-
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=treatment)) + 
#   geom_point(alpha=0.8) + 
#   theme_classic() +
#   scale_color_brewer(palette="Dark2") +
#   theme(legend.title=element_blank())
# abund_nmds_plant
# 
# abund_nmds_loc <- 
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=location)) + 
#   geom_point(alpha=0.8) + 
#   theme_classic() +
#   scale_color_brewer(palette="Set1") +
#   theme(legend.title=element_blank())
# abund_nmds_loc
# 
# abund_nmds_gwc <- 
#   ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=gwc)) +
#   geom_point(alpha=0.8) +
#   theme_classic() +
#   viridis::scale_color_viridis(option="B", name="GWC") 
# abund_nmds_gwc

