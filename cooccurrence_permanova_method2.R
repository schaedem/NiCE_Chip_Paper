library(tidyverse)
library(vegan)
library(reshape2)
library(bioDist)
source("permanova_community_differences.R")

dist_spearman <- spearman.dist(nice_dist_mat)

spearman_dist_df <- as.matrix(dist_spearman) %>%
  as.data.frame()

spearman_merge <-  rownames_to_column(spearman_dist_df, var="sample") %>% 
  as_tibble() %>%
  merge(meta, by="sample")

#permanova testing
set.seed(123)
spearman_loc_timepoint <- 
  adonis2(dist_spearman ~ timepoint*location, data=spearman_merge, permutations=9999)
#timepoint highly sig, R2=0.39559, p=0.0001
#location sig, R2=0.07318, p=0.0001
#interaction between timepoint and loc is sig, R2=0.015, p=5e-4

set.seed(123)
# spearman_loc_timepoint_block<- 
#   adonis2(dist_spearman ~ timepoint*location*block, 
#                                 data=spearman_merge, permutations=9999)
# spearman_loc_timepoint_block
#block (p = 0.6098) and block interactions are not significant

set.seed(123)
# spearman_all_adonis <- adonis2(dist_spearman~timepoint*location*treatment, 
#                                dataspearman_merge, permutations=9999)
# spearman_all_adonis
#treatment (p=0.86) and treatment interactions not significant

#### pairwise comparisons to find sig differences between timepoints ###
pairwise_df_spearman <- numeric() #initialize empty df to put p-values

pairwise_spearman_adonis <- function(val1, val2, adonis_merge) {
  val1 <- val1
  val2 <- val2
  
  filter_merge <- adonis_merge %>%
    filter(timepoint == val1 | timepoint == val2)
  
  dist_filter <- filter_merge %>%
    select(all_of(.[["sample"]])) 
    #spearman.dist(as.matrix(.))
  
  set.seed(123)
  result <- adonis(dist_filter ~ timepoint, data=filter_merge)
  
  result_p <- result[["aov.tab"]][["Pr(>F)"]][1]
  
  return(result_p)
  
}

pairwise_1_2 <- pairwise_spearman_adonis(1,2, spearman_merge)
pairwise_1_3 <- pairwise_spearman_adonis(1,3, spearman_merge) 
pairwise_1_4 <- pairwise_spearman_adonis(1,4, spearman_merge) 
pairwise_1_5 <- pairwise_spearman_adonis(1,5, spearman_merge) 
pairwise_1_6 <- pairwise_spearman_adonis(1,6, spearman_merge) 
pairwise_2_3 <- pairwise_spearman_adonis(2,3, spearman_merge) 
pairwise_2_4 <- pairwise_spearman_adonis(2,4,spearman_merge)
pairwise_2_5 <- pairwise_spearman_adonis(2,5, spearman_merge)
pairwise_2_6 <- pairwise_spearman_adonis(2,6, spearman_merge)
pairwise_3_4 <- pairwise_spearman_adonis(3,4, spearman_merge)
pairwise_3_5 <- pairwise_spearman_adonis(3,5, spearman_merge)
pairwise_3_6 <- pairwise_spearman_adonis(3,6, spearman_merge)
pairwise_4_5 <- pairwise_spearman_adonis(4,5, spearman_merge)
pairwise_4_6 <- pairwise_spearman_adonis(4,6,spearman_merge)
pairwise_5_6 <- pairwise_spearman_adonis(5,6,spearman_merge)

pairwise_df_spearman[c("1-2", "1-3", "1-4", "1-5", "1-6", "2-3", "2-4", "2-5", "2-6",
              "3-4", "3-5", "3-6", "4-5", "4-6", "5-6")] <- 
  c(pairwise_1_2, pairwise_1_3, pairwise_1_4,pairwise_1_5,pairwise_1_6,
    pairwise_2_3, pairwise_2_4 ,pairwise_2_5 ,pairwise_2_6 ,pairwise_3_4 ,pairwise_3_5 ,
    pairwise_3_6 ,pairwise_4_5 , pairwise_4_6 ,pairwise_5_6  )

pairwise_p <- p.adjust(pairwise_df_spearman, method="BH") %>% #adjust for multiple comparisons
  as.data.frame() %>%
  rownames_to_colnames()
all(pairwise_p < 0.05)


#NMDS plots
set.seed(100192)
nmds_spear <- metaMDS(dist_spearman)

nmds_plot_spear <- scores(nmds_spear, display="sites") %>%
  as_tibble(rownames="sample") %>%
  merge(meta, by="sample")

spear_nmds_growth <- 
  ggplot(data=nmds_plot_spear, aes(x=NMDS1, y=NMDS2, color=growth)) + 
  geom_point(alpha=0.8) + 
  theme_classic() +
  scale_color_brewer(palette="Set2") +
  theme(legend.title=element_blank())
spear_nmds_growth

spear_nmds_time <- 
  ggplot(data=nmds_plot_spear, aes(x=NMDS1, y=NMDS2, color=as.factor(timepoint))) + 
  stat_ellipse(level=0.8, geom="polygon", alpha=0,show.legend=FALSE, 
               aes(color=as.factor(timepoint)))+
  geom_point() + 
  theme_classic() +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") 
spear_nmds_time
ggsave("spear_nmds_time.tiff", dpi=300)
# spear_nmds_plant <-
#   ggplot(data=nmds_plot_spear, aes(x=NMDS1, y=NMDS2, color=treatment)) + 
#   geom_point(alpha=0.8) + 
#   theme_classic() +
#   scale_color_brewer(palette="Dark2") +
#   theme(legend.title=element_blank())
# spear_nmds_plant

spear_nmds_loc <- 
  ggplot(data=nmds_plot_spear, aes(x=NMDS1, y=NMDS2, color=location)) + 
  stat_ellipse(level=0.9, geom="polygon", alpha=0,show.legend=FALSE, 
               aes(color=location)) +
  geom_point(alpha=0.8) + 
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  theme(legend.title=element_blank())
spear_nmds_loc

ggsave("spear_nmds_loc.tiff", dpi=300)

#collate for Fig 2

fig2 <- final_fig <-
  cowplot::plot_grid(spear_nmds_time, spear_nmds_loc,
                     labels="AUTO")


set.seed(123)
envfit(nmds_spear ~ location*as.factor(timepoint), spearman_merge, perm=999)
#location R2 = 0.0232, p=0.001
#timepoint R2 = 0.55, p=0.001