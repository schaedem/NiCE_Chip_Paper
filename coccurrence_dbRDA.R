library(tidyverse)
library(vegan)
library(SNCAnalysis)
source("permanova_community_differences.R")

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

# #matrix of covariables (continuous environmental variables)
# covariables <- env_data %>%
#   select(-c(location, location, block, timepoint, sample, treatment))
# rownames(covariables) <- env_data$sample
# head(covariables)
# 
# covariables_mat <- covariables %>%
#   as.matrix()

#distance-based Redundancy Analysis (db-RDA)
#takes the spearman distance matrix as input
#multivariate multiple linear regression
#test the relative contributions of soil/environmental variables on gene co-occurrence
#db-RDA can be performed either with the capscale() or dbrda() functions in vegan

#capscale() computes a PCoA of the dissimilarity matrix and runs an RDA of the principle coordinates
# constrained by the explanatory variables

#dbrda() runs the RDA directly on the dissimilarity response matrix, avoiding the PCoA step
#this has benefits for multifactorial ANOVA designs, with correct type I error. 

# spear.dbrda <- dbrda(spearman_dist_mat ~ timepoint + location + timepoint:location + 
#                        Condition(covariables_mat), data=env_data, add=FALSE)
# anova(spear.dbrda, permutations=how(nperm=999))

########## All timepoints, both locations ############
#after accounting for main factors, is there any other variable that accounts for a sig portion of the variation?
#full model
spear.dbrda.full <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location + gwc + pH + treatment + log_dea +
                            mg_kg_NO3N + mg_kg_NH4N + NH4N_net + mg_kg_POXC + NP, data=env_data, add=FALSE)
anova(spear.dbrda.full, permutations=how(nperm=999))

#variable selection
mod0 <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location, data=env_data, add=FALSE)

spear.dbrda.step <- ordiR2step(mod0, scope=formula(spear.dbrda.full), direction="forward", permutations=how(nperm=999))

#with timepoint and location
mod2 <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location + mg_kg_POXC, data=env_data, add=FALSE)
anova(mod0, mod2, permutations=how(nperm=999))

mod3 <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location + mg_kg_POXC + NP, data=env_data, add=FALSE )

mod4 <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location + mg_kg_POXC + log_dea, data=env_data, add=FALSE )

mod5 <- dbrda(spearman_dist_mat ~ as.factor(timepoint) + location + mg_kg_POXC + NP + log_dea, data=env_data, add=FALSE )

anova(mod0, mod2, mod3, mod4, mod5) #model with only POXC as additional variable is preferred (mod2)

#check model summary and statistics
summary(mod2)
screeplot(mod2)
vif.cca(mod2) # variance inflation factors. Anything above 10 should be examined or avoided
vif.cca(spear.dbrda.full)
coef(mod2)
R2 <- RsquareAdj(mod2)$r.squared
R2 #0.6287
R2adj <- RsquareAdj(mod2)$adj.r.squared
R2adj #0.6207

#compare to R2 mod0
RsquareAdj(mod0)$r.squared
RsquareAdj(mod0)$adj.r.squared

#capscale for plotting species scores
spear.cap <- capscale(spearman_dist_mat ~ as.factor(timepoint) + location + mg_kg_POXC, data=env_data, add="lingoes")
plot(spear.cap, scaling=2)

#plot capscale results
cap.scores <- scores(spear.cap)
cap.sites <- cap.scores$sites%>%
  as.data.frame()%>%
  rownames_to_column("sample") %>%
  mutate(sample=as.numeric(sample)) %>%
  full_join(env_data, by="sample")
head(cap.sites)

cap.cent <- cap.scores$centroids %>%
  as.data.frame() %>%
  rownames_to_column("label")
cap.cent$label <- c("Timepoint 1", "Timepoint 2", "Timepoint 3", "Timepoint 4", "Timepoint 5", "Timepoint 6",
                    "Karama", "Rubona")

cap.arrow <- cap.scores$biplot %>%
  as.data.frame() %>%
  rownames_to_column("label") %>%
  filter(label=="mg_kg_POXC")

#plot capscale

plot_cap <- ggplot() +
  geom_point(data=cap.sites, aes(x=CAP1, y=CAP2, color=as.factor(timepoint), shape=location), alpha=0.6) +
  theme_classic() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(data=cap.arrow, aes(x=0,y=0, xend=CAP1*13.35134, yend=CAP2*13.35134), arrow=arrow(length=unit(0.2, "cm"))) +
  geom_text(data = cap.arrow, aes(x = CAP1*13.35134, y = CAP2*13.35134, label = label),
            nudge_x = 0.11, nudge_y=-0.08) +
  xlab("dbRDA1") +
  ylab("dbRDA2") +
  geom_point(data=cap.cent, aes(x=CAP1, y=CAP2), color="black", size=3, stroke = 2,shape=3) +
  geom_text(data = cap.cent, aes(x= CAP1, y=CAP2, label=label), nudge_x=0.3, cex=4, fontface="bold")

plot_cap #equivalent to scaling 2


#plot dbRDA results
rda.scores <- scores(mod2, display="sites") %>%
  as.data.frame() %>%
  rownames_to_column("sample") %>%
  mutate(sample=as.numeric(sample)) %>%
  full_join(env_data, by="sample")
head(rda.scores)

rda.vect <- as.matrix(scores(mod2, display="bp", scaling="species")) %>% #no species scores yet - need these to plot vectors
  as.data.frame()
rownames(rda.vect) <- c("Timepoint 1", "Timepoint 2", "Timepoint 3", "Timepoint 4", "Timepoint 5", "Location", "POXC")

plot_dbRDA <- ggplot() +
  geom_point(data=rda.scores, aes(x=dbRDA1, y=dbRDA2, color=as.factor(timepoint), shape=location)) +
  theme_classic() +
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +
  viridis::scale_color_viridis(discrete="TRUE", name="Timepoint") +
  geom_segment(data=rda.vect, aes(x=0,y=0, xend=dbRDA1, yend=dbRDA2), arrow=arrow(length=unit(0.2, "cm"))) +
  geom_text(data = rda.vect, aes(x = dbRDA1, y = dbRDA2, label = rownames(rda.vect)), nudge_x = 0.11, nudge_y=0.03) 

plot_dbRDA

ggsave("co-occurrence_dbRDA.tiff", dpi=300, width=9, height=5)

##########looking at distribution of env variables by time and location##########

ggplot(aes(x=as.factor(timepoint), y=NH4N_net), data=env_data) + #PMN decreased in early rainy season, highest in dry and late rainy
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location)

ggplot(aes(x=as.factor(timepoint), y=mg_kg_NO3N), data=env_data) + #NO3 decreased
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location)

ggplot(aes(x=as.factor(timepoint), y=mg_kg_NH4N), data=env_data) + #NH4 decreased (Rubona) or stayed the same
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location)

ggplot(aes(x=as.factor(timepoint), y=NP), data=env_data) + #np generally increased in rainy season
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location)

ggplot(aes(x=as.factor(timepoint), y=log_dea), data=env_data) + #dea generally increased in rainy season
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location)

ggplot(aes(x=as.factor(timepoint), y=mg_kg_POXC), data=env_data, fill=as.factor(timepoint)) + #poxc was potentially more related to growth stage; effects varied based on location 
  geom_boxplot() +
  theme_classic() +
  facet_wrap(~location) 


ggplot(aes(x=mg_kg_POXC, y=log_dea, group=location), data=env_data) +
  geom_point(aes(fill=as.factor(timepoint))) +
  facet_wrap(~location) +
  theme_classic() +
  geom_smooth(method="lm") +
  viridis::scale_fill_viridis(discrete="TRUE", name="Timepoint") 

ggplot(aes(x=mg_kg_POXC, y=NP), data=env_data) +
  geom_point() +
  facet_wrap(~location) +
  theme_classic() +
  geom_smooth(method="lm")

#POXC relates to potential activity, but the strength of the effect varies

ggplot(aes(x=mg_kg_POXC, y=log_dea, color=location), data=env_data) +
  geom_point() +
  facet_wrap(~timepoint) +
  theme_classic() +
  geom_smooth(method="lm", se=FALSE) +
  ylab(expression(~"ln DEA("~"nmol"~N[2]~"O"~""~kg^-1~"soil"~hr^-1~")")) +
  xlab(expression(~"mg"~"POX-C"~kg^-1~"soil")) +
  scale_color_brewer(palette="Set1") +
  theme(legend.title = element_blank())+
  ggpmisc::stat_poly_eq(formula=y~x)

ggplot(aes(x=mg_kg_POXC, y=NP, color=location), data=env_data) +
  geom_point() +
  facet_wrap(~timepoint) +
  theme_classic() +
  geom_smooth(method="lm", se=FALSE) +
  ylab(expression(~"NP"~"(mg"~NO[3]^-~""~kg^-1~"soil)")) +
  xlab(expression(~"mg"~"POX-C"~kg^-1~"soil")) +
  scale_color_brewer(palette="Set1") +
  theme(legend.title = element_blank()) +
  ggpmisc::stat_poly_eq(formula=y~x)
