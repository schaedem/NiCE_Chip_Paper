library(tidyverse)
library(vegan)
library(SNCAnalysis)
source("coccurrence_dbRDA_2.R")

#dbRDA within timepoint
t1_dat <- spearman_merge %>%
  filter(timepoint ==1)
t1_merge <- t1_dat[,2:56]

t2_dat <- spearman_merge %>%
  filter(timepoint==2)
t2_merge <- t2_dat[,57:112]

t3_dat <- spearman_merge %>%
  filter(timepoint==3)
t3_merge <- t3_dat[,113:168]

t4_dat <- spearman_merge %>%
  filter(timepoint==4)
t4_merge <- t4_dat[,169:224]

t5_dat <- spearman_merge %>%
  filter(timepoint==5)
t5_merge <- t5_dat[,225:280]

t6_dat <- spearman_merge %>%
  filter(timepoint==6) 
t6_merge <- t6_dat[,281:336]

#small function to return merged dfs back to matrix format
make_matrix <- function(merge_df) {
  merge_df <- merge_df
  rownames(merge_df) <- colnames(merge_df)
  new_mat <- as.matrix(merge_df)
  return(new_mat)
}

t1_distmat <- make_matrix(t1_merge)
t2_distmat <- make_matrix(t2_merge)
t3_distmat <- make_matrix(t3_merge)
t4_distmat <- make_matrix(t4_merge)
t5_distmat <- make_matrix(t5_merge)
t6_distmat <- make_matrix(t6_merge)

########## Within timepoints, both locations ############
#after accounting for main factors, is there any other variable that accounts for a sig portion of the variation?
#empty models
t1.dbrda.null <- dbrda(t1_distmat ~ 1, t1_dat)
t2.dbrda.null <- dbrda(t2_distmat ~ 1, t2_dat)
t3.dbrda.null <- dbrda(t3_distmat ~ 1, t3_dat)
t4.dbrda.null <- dbrda(t4_distmat ~ 1, t4_dat)
t5.dbrda.null <- dbrda(t5_distmat ~ 1, t5_dat)
t6.dbrda.null <- dbrda(t6_distmat ~ 1, t6_dat)

#full models
full_dbrda <- function(distmat, envdat) {
  distmat <- distmat
  envdat <- envdat
  
  dbrda.full <- dbrda(distmat ~  gwc + pH + treatment + log_dea + mg_kg_NO3N + mg_kg_NH4N
                      + NH4N_net + mg_kg_POXC + NP, envdat, add = FALSE)
  return(dbrda.full)
  
}

t1.dbrda.full <- full_dbrda(t1_distmat, t1_dat)
anova(t1.dbrda.full, permutations = how(nperm=999))
vif.cca(t1.dbrda.full) #all below 10

t2.dbrda.full <- full_dbrda(t2_distmat, t2_dat)
anova(t2.dbrda.full, permutations = how(nperm=999))
vif.cca(t2.dbrda.full) #all below 10

t3.dbrda.full <- full_dbrda(t3_distmat, t3_dat)
anova(t3.dbrda.full, permutations = how(nperm=999))
vif.cca(t3.dbrda.full) #all below 10

t4.dbrda.full <- full_dbrda(t4_distmat, t4_dat)
anova(t4.dbrda.full, permutations = how(nperm=999))
vif.cca(t4.dbrda.full) #all below 10

t5.dbrda.full <- full_dbrda(t5_distmat, t5_dat)
anova(t5.dbrda.full, permutations = how(nperm=999))
vif.cca(t5.dbrda.full) #all below 10

t6.dbrda.full <- full_dbrda(t6_distmat, t6_dat)
anova(t6.dbrda.full, permutations = how(nperm=999))
vif.cca(t6.dbrda.full) #all below 10

#variable selection
step_fun <- function(full_mod, null_mod) {
  null_mod <- null_mod
  full_mod <- full_mod
  spear.dbrda.step <- ordiR2step(null_mod, scope=formula(full_mod), direction="forward", 
                                 permutations=how(nperm=999))
  return(spear.dbrda.step)
}

t1.dbrda.step <- step_fun(t1.dbrda.full, t1.dbrda.null) #nothing significant
t1.mod2 <- dbrda(t1_distmat ~ 1, t1_dat)
t2.dbrda.step <- step_fun(t2.dbrda.full, t2.dbrda.null) #nothing significant
t2.mod2 <- dbrda(t1_distmat ~ 1, t2_dat)

t3.dbrda.step <- step_fun(t3.dbrda.full, t3.dbrda.null) #pH, mg_kg_NO3N; R2=0.458
t3.mod2 <- dbrda(t3_distmat ~ pH + mg_kg_NO3N, t3_dat)
vif.cca(t3.mod2)
screeplot(t3.mod2)
ordiplot(t3.mod2)
anova(t3.mod2, by='axis')

t4.dbrda.step <- step_fun(t4.dbrda.full, t4.dbrda.null) #gwc, pH, NP, POXC
t4.mod2 <- dbrda(t4_distmat ~ pH + gwc + NP + mg_kg_POXC, t4_dat)
summary(t4.mod2)
vif.cca(t4.mod2)
screeplot(t4.mod2)
ordiplot(t4.mod2)
anova(t4.mod2, by='axis')

t5.dbrda.step <- step_fun(t5.dbrda.full, t5.dbrda.null) #log_dea, NP, pH
t5.mod2 <- dbrda(t5_distmat ~ log_dea + NP + pH, t5_dat)
anova(t5.mod2, by='axis')
summary(t5.mod2)
vif.cca(t5.mod2)
screeplot(t5.mod2)
ordiplot(t5.mod2)

t6.dbrda.step <- step_fun(t6.dbrda.full, t6.dbrda.null) #error??
t6.mod2 <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NH4N + NP, t6_dat, add= FALSE)
vif.cca(t6.mod2)
screeplot(t6.mod2)
ordiplot(t6.mod2)
anova(t6.mod2, by='axis')

#####manual forward selection for t6  since ordiR2step is throwing an error (likely to do with a data object or other package) ######
t6.dbrda.null <- dbrda(t6_distmat ~ 1, t6_dat)
t6.dbrda.full <- dbrda(t6_distmat ~ gwc + pH + treatment + log_dea + mg_kg_NO3N + mg_kg_NH4N
                       + NH4N_net + mg_kg_POXC + NP, t6_dat, add = FALSE)
t6.dbrda.2var <- dbrda(t6_distmat ~ 1 + gwc, t6_dat, add=FALSE)
anova(t6.dbrda.null, t6.dbrda.2var) #keep gwc

t6.dbrda.3var <- dbrda(t6_distmat ~ gwc + pH, data=t6_dat, add=FALSE)
anova(t6.dbrda.2var, t6.dbrda.3var) #keep pH

t6.dbrda.4var <- dbrda(t6_distmat ~ gwc + pH + treatment, data=t6_dat, add=FALSE)
anova(t6.dbrda.3var, t6.dbrda.4var) #don't include treatment

t6.dbrda.5var <- dbrda(t6_distmat ~ gwc + pH + log_dea, t6_dat, add=FALSE)
anova(t6.dbrda.5var, t6.dbrda.3var) #don't include log_dea

t6.dbrda.6var <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NO3N, t6_dat, add=FALSE)
anova(t6.dbrda.6var, t6.dbrda.3var) #don't include NO3

t6.dbrda.7var <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NH4N, t6_dat, add=FALSE)
anova(t6.dbrda.7var, t6.dbrda.3var) #keep NH4N

t6.dbrda.8var <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NH4N + NH4N_net, t6_dat, add=FALSE)
anova(t6.dbrda.7var, t6.dbrda.8var) #don't keep PMN

t6.dbrda.9var <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NH4N + mg_kg_POXC, t6_dat, add=FALSE)
anova(t6.dbrda.7var, t6.dbrda.9var) #borderline, but don't keep POXC

t6.dbrda.10var <- dbrda(t6_distmat ~ gwc + pH + mg_kg_NH4N + NP, t6_dat, add= FALSE)
anova(t6.dbrda.7var, t6.dbrda.10var) #keep NP; final model

#Plot dbRDAs


#plot dbRDA results

plot_result <- function(model, data) {
  rda.scores <- scores(model, display="sites") %>% #constant 5.437482
    as.data.frame() %>%
    rownames_to_column("sample") %>%
   # mutate(sample=as.numeric(sample)) %>%
    full_join(data, by="sample") %>%
    rename(Location = location)
  
  rda.vect <- as.matrix(scores(model, display="bp", scaling="species")) %>% #no species scores yet - need these to plot vectors
    as.data.frame() %>%
    rownames_to_column(var="label")

  rda_plot <- ggplot() +
    geom_point(data=rda.scores, aes(x=dbRDA1, y=dbRDA2, shape=Location, color=Location)) +
    theme_classic() +
    scale_color_brewer(palette="Set1") +
    xlim(-1,1) +
    ylim(-2,2) +
    geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
    geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
    geom_segment(data=rda.vect, aes(x=0,y=0, xend=dbRDA1, yend=dbRDA2), arrow=arrow(length=unit(0.2, "cm"))) +
    geom_text(data = rda.vect, aes(x = dbRDA1, y=dbRDA2, label = label), nudge_y = -0.08)
  
  return(rda_plot)
  
}

t3_rda <- plot_result(t3.mod2, t3_dat) + ggtitle("Timepoint 3: mid rainy/ anthesis")
t4_rda <- plot_result(t4.mod2, t4_dat)+ ggtitle("Timepoint 4: mid rainy/ regrowth")
t5_rda <- plot_result(t5.mod2, t5_dat)+ ggtitle("Timepoint 5: late rainy/ anthesis")
t6_rda <- plot_result(t6.mod2, t6_dat)+ ggtitle("Timepoint 6: late rainy/ regrowth")

#### plot unconstrained ordination for timepoints 1 and 2 ####
rda.scores.t1 <- scores(t1.dbrda.null, display="sites") %>% #constant 5.437482
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  # mutate(sample=as.numeric(sample)) %>%
  full_join(t1_dat, by="sample") %>%
  rename(Location = location)

rda.scores.t2 <- scores(t2.dbrda.null, display="sites") %>% #constant 5.437482
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  # mutate(sample=as.numeric(sample)) %>%
  full_join(t2_dat, by="sample") %>%
  rename(Location = location)


t1_rda <- ggplot(data=rda.scores.t1, aes(x=MDS1, y=MDS2, shape=Location, color=Location), alpha=0.6) +
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  xlim(-1,1) +
  ylim(-2,2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  ggtitle("Timepoint 1: dry season / anthesis")
t1_rda


t2_rda <- ggplot(data=rda.scores.t2, aes(x=MDS1, y=MDS2, shape=Location, color=Location), alpha=0.6) +
  geom_point() +
  theme_classic() +
  scale_color_brewer(palette="Set1") +
  xlim(-1,1) +
  ylim(-2,2) +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  ggtitle("Timepoint 2: dry season / regrowth")
t2_rda

#collate dbRDAs into one figure
final_fig <-
  cowplot::plot_grid(t1_rda + theme(legend.position="none"), 
                   t2_rda + theme(legend.position="none"),
                   t3_rda + theme(legend.position="none"),
                   t4_rda + theme(legend.position="none"),
                   t5_rda + theme(legend.position="none"), 
                   t6_rda + theme(legend.position="none"),
                   labels="AUTO")

# extract the legend from one of the plots
legend <- cowplot::get_legend(
  # create some space to the left of the legend
  t1_rda + theme(legend.box.margin = margin(0, 0, 0, 12))
)

# add the legend to the row we made earlier. Give it one-third of 
# the width of one plot (via rel_widths).
cowplot::plot_grid(final_fig, legend, rel_widths = c(3, .4))

ggsave("timepoint_dbrda.tiff", dpi=300, width=12, height=6)
