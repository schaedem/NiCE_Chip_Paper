library(tidyverse)
library(vegan)
library(reshape2)
library(bioDist)
source("permanova_community_differences.R")

dist_spearman <- spearman.dist(nice_dist_mat)

#permanova testing
set.seed(123)
spearman_loc_timepoint <- 
  adonis2(dist_spearman ~ timepoint*location, data=adonis_merge, permutations=9999)
#timepoint highly sig, R2=0.39559, p=0.0001
#location sig, R2=0.07318, p=0.0001
#interaction between timepoint and loc is sig, R2=0.015, p=5e-4

set.seed(123)
spearman_loc_timepoint_block<- 
  adonis2(dist ~ timepoint*location*block, 
                                data=adonis_merge, permutations=9999)
spearman_loc_timepoint_block
#block (p = 0.6098) and block interactions are not significant

set.seed(123)
spearman_all_adonis <- adonis2(dist~timepoint*location*treatment, 
                               data=adonis_merge, permutations=9999)
spearman_all_adonis
#treatment (p=0.86) and treatment interactions not significant

####
#pairwise comparisons to find sig differences between timepoints
pairwise_df_spearman <- numeric() #initialize empty df to put p-values

pairwise_spearman_adonis <- function(val1, val2, adonis_merge) {
  val1 <- val1
  val2 <- val2
  
  filter_merge <- adonis_merge %>%
    filter(timepoint == val1 | timepoint == val2)
  
  dist_filter <- filter_merge %>%
    select(all_of(.[["sample"]])) 
    #spearman.dist(as.matrix(.))
  
 # set.seed(123)
 # result <- adonis(dist_filter ~ timepoint, data=filter_merge)
  
 # result_p <- result[["aov.tab"]][["Pr(>F)"]][1]
  
  #return(result_p)
  
}

pairwise_1_2 <- pairwise_spearman_adonis(1,2, adonis_merge)
pairwise_1_3 <- pairwise_adonis(1,3, adonis_merge) 
pairwise_1_4 <- pairwise_adonis(1,4, adonis_merge) 
pairwise_1_5 <- pairwise_adonis(1,5, adonis_merge) 
pairwise_1_6 <- pairwise_adonis(1,6, adonis_merge) 
pairwise_2_3 <- pairwise_adonis(2,3, adonis_merge) 
pairwise_2_4 <- pairwise_adonis(2,4,adonis_merge)
pairwise_2_5 <- pairwise_adonis(2,5, adonis_merge)
pairwise_2_6 <- pairwise_adonis(2,6, adonis_merge)
pairwise_3_4 <- pairwise_adonis(3,4, adonis_merge)
pairwise_3_5 <- pairwise_adonis(3,5, adonis_merge)
pairwise_3_6 <- pairwise_adonis(3,6, adonis_merge)
pairwise_4_5 <- pairwise_adonis(4,5, adonis_merge)
pairwise_4_6 <- pairwise_adonis(4,6,adonis_merge)
pairwise_5_6 <- pairwise_adonis(5,6,adonis_merge)

pairwise_df[c("1-2", "1-3", "1-4", "1-5", "1-6", "2-3", "2-4", "2-5", "2-6",
              "3-4", "3-5", "3-6", "4-5", "4-6", "5-6")] <- 
  c(pairwise_1_2, pairwise_1_3, pairwise_1_4,pairwise_1_5,pairwise_1_6,
    pairwise_2_3, pairwise_2_4 ,pairwise_2_5 ,pairwise_2_6 ,pairwise_3_4 ,pairwise_3_5 ,
    pairwise_3_6 ,pairwise_4_5 , pairwise_4_6 ,pairwise_5_6  )

pairwise_p <- p.adjust(pairwise_df, method="BH") #adjust for multiple comparisons
all(pairwise_p < 0.05)


#NMDS plots
set.seed(100192)
nmds <- metaMDS(dist)

nmds_plot <- scores(nmds) %>%
  as_tibble(rownames="sample") %>%
  merge(meta, by="sample")

abund_nmds_growth <- 
  ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=growth)) + 
  geom_point(alpha=0.8) + 
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(legend.title=element_blank())
abund_nmds_growth

abund_nmds_time <- 
  ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) + 
  geom_point() + 
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint")
abund_nmds_time

abund_nmds_plant <-
  ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=treatment)) + 
  geom_point(alpha=0.8) + 
  theme_classic() +
  scale_color_brewer(palette="Dark2") +
  theme(legend.title=element_blank())
abund_nmds_plant

abund_nmds_loc <- 
  ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=location)) + 
  geom_point(alpha=0.8) + 
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(legend.title=element_blank())
abund_nmds_loc

abund_nmds_gwc <- 
  ggplot(data=nmds_plot, aes(x=NMDS1, y=NMDS2, color=gwc)) +
  geom_point(alpha=0.8) +
  theme_classic() +
  viridis::scale_color_viridis(option="B", name="GWC") 
abund_nmds_gwc

