library(tidyverse)
library(vegan)
library(SNCAnalysis)
source("permanova_community_differences.R")

#without timepoint or location to get globally important soil factors

#format the dissimilarity response matrix
dist_spearman <- spearman.dist(nice_dist_mat) #input is standardized nice data in 335x20 matrix

spearman_dist_mat <- as.matrix(dist_spearman) 

spearman_dist_df <- as.matrix(dist_spearman) %>%
  as.data.frame()

spearman_merge <-  rownames_to_column(spearman_dist_df, var="sample") %>% 
  as_tibble() %>%
  merge(meta, by="sample")

#format the constraining environmental variable dataset
env_data <- full_dat %>%
  select(-c(auc_np, auc_dea, log_auc_dea, mmol_kg_hr))
colnames(env_data)

#merge dist mat with metadata to filter by timepoint
spearman_merge <-  rownames_to_column(spearman_dist_df, var="sample") %>%
  as_tibble() %>%
  merge(env_data, by="sample")

#variable selection
#null model
mod00 <- dbrda(spearman_dist_mat ~ 1, data=env_data, add=FALSE)

#full model
spear.dbrda.full2 <- dbrda(spearman_dist_mat ~  gwc + pH + treatment + log_dea +
                             mg_kg_NO3N + mg_kg_NH4N + NH4N_net + mg_kg_POXC + NP, data=env_data, add=FALSE)
anova(spear.dbrda.full2, permutations=how(nperm=999))
vif.cca(spear.dbrda.full2)

#forward selection
spear.dbrda.step2 <- ordiR2step(mod00, scope = formula(spear.dbrda.full2), direction="forward", permutations=how(nperm=999))

#final model
fullmod <- dbrda(spearman_dist_mat ~ mg_kg_NO3N + NH4N_net + pH + log_dea + gwc + 
                   mg_kg_POXC + mg_kg_NH4N, data=env_data, add=FALSE)
summary(fullmod)
anova(fullmod)
screeplot(fullmod)
vif.cca(fullmod) # variance inflation factors. Anything above 10 should be examined or avoided
coef(fullmod)
RsquareAdj(fullmod)$r.squared #0.3953
RsquareAdj(fullmod)$adj.r.squared #0.3823
anova(fullmod)

#plot dbRDA results
rda.scores <- scores(fullmod, display="sites") %>% #constant 5.437482
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample=as.numeric(sample)) %>%
  full_join(env_data, by="sample") %>%
  rename(Location = location)
head(rda.scores)

rda.vect <- as.matrix(scores(fullmod, display="bp", scaling="species")) %>% #no species scores yet - need these to plot vectors
  as.data.frame()
#rownames(rda.vect) <- c("Timepoint 1", "Timepoint 2", "Timepoint 3", "Timepoint 4", "Timepoint 5", "Location", "POXC")

plot_dbRDA <- ggplot() +
  geom_point(data=rda.scores, aes(x=dbRDA1, y=dbRDA2, color=as.factor(timepoint),shape=Location), alpha=0.6) +
  theme_classic() +
  xlim(-1,1) +
  ylim(-2,2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(data=rda.vect, aes(x=0,y=0, xend=dbRDA1, yend=dbRDA2), arrow=arrow(length=unit(0.2, "cm"))) +
  annotate("text", x=-0.8, y= -0.4, label = "NO[3]^'-'-N", parse=TRUE)+
  annotate("text", x = -0.33, y=-0.1, label = "pH") +
  annotate("text", x= -0.35, y=0.1, label = "NH[4]^'+'-N", parse=TRUE) +
  annotate("text", x=-0.1, y=-0.33, label = "POX-C") +
  annotate("text", x=0.25, y=-1, label = "PMN") +
  annotate("text", x= 0.63, y=0.15, label = "GWC") +
  annotate("text", x=0.6, y=0, label= "ln(DEA)")
plot_dbRDA

plot(fullmod, scaling=1)

ggsave("co-occurrence_dbRDA_full.tiff", dpi=300, width=8, height=6)

