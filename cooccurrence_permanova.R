library(bioDist)
library(vegan)
library(mvtnorm)
library(ggplot2)
library(igraph)
library(SNCAnalysis)
library(reshape2)
source("correlation_dataframes.R")
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/metadata")
meta <- read_csv("SNC_metadata_2020.csv") 
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")

set.seed(344)

env_data <- SNCAnalysis::full_dat %>%
  subset(select= c(timepoint, location, sample))

nice_all <- rbind(r_1, k_1, r_2, k_2, r_3, k_3, r_4, k_4, r_5, k_5, r_6, k_6) %>%
  select(-assay) %>%
  merge(full_dat, by=c("sample", "timepoint", "location"))

nice_spread <- nice_all %>%
  pivot_wider(id_cols = c(sample, timepoint, treatment, location, block),  
              names_from=primer_pair, values_from=log_conc_truesoil)

nice_melt <-melt(nice_spread, id=c("sample","timepoint", "treatment", "location", "block")) %>% 
  rename(primer=variable)

#matrices are made for each ecosystem type and recombined for the analysis so that there is a
#column of species id's and samples are the y variables

t1<-dcast(subset(nice_melt, timepoint=="1"), timepoint+primer~sample, value="value") %>%
  mutate("246"="NA")
t1 <- t1[,c(1:47,58,48:57)]
t2<-dcast(subset(nice_melt, timepoint=="2"), timepoint+primer~sample, value="value")
t3<-dcast(subset(nice_melt, timepoint=="3"), timepoint+primer~sample, value="value")
t4<-dcast(subset(nice_melt, timepoint=="4"), timepoint+primer~sample, value="value")
t5<-dcast(subset(nice_melt, timepoint=="5"), timepoint+primer~sample, value="value")
t6<-dcast(subset(nice_melt, timepoint=="6"), timepoint+primer~sample, value="value")

colnames(t6) <- colnames(t2)
colnames(t5) <- colnames(t2)
colnames(t4) <- colnames(t2)
colnames(t3) <- colnames(t2)
colnames(t1) <- colnames(t2)

#to be used in spearman.dist(), which calculates correlational distance for all *rows* of a matrix
#this means that the microbes have to be x's and the samples are the y's (columns)
nice_cocur<-data.frame(rbind(t1,t2, t3, t4, t5, t6))

trts<-as.vector((unique((nice_all$timepoint))))

#making the spearman's distance matrix
nice_cocur_dist<-spearman.dist(data.matrix(nice_cocur[,-c(1:2)]))

#PERMANOVA
set.seed(222)
test <- adonis(nice_cocur_dist~nice_cocur$timepoint, permutations=9999) 
summary(test)#timepoint highly sig (p=1e-4; R2=0.13768)

#NMDS
set.seed(213)
nice_mds<-monoMDS(nice_cocur_dist, model="local")
plot(nice_mds)
stressplot(nice_mds)

#Repeat same process to test whether co-occurrence relationships differ by location
Rub_dist <-dcast(subset(nice_melt, location=="Rubona"), location+primer~sample, value="value") 
Kar_dist <- dcast(subset(nice_melt, location=="Karama"), location+primer~sample, value="value")%>%
  mutate("246"="NA") #missing data point

Kar_dist <- Kar_dist[,c(1:19,170,20:169)]

colnames(Kar_dist) <- colnames(Rub_dist)

loc_cocur <- data.frame(rbind(Rub_dist, Kar_dist))

locs<-as.vector((unique((nice_all$location))))

#making the spearman's distance matrix
loc_cocur_dist<-spearman.dist(data.matrix(loc_cocur[,-c(1:2)]))

#PERMANOVA
set.seed(223)
adonis(loc_cocur_dist~loc_cocur$location, permutations=9999) 
#loc highly sig (p=1e-4; R2=0.45433)

#NMDS
set.seed(213)
nice_mds<-monoMDS(loc_cocur_dist, model="local")
plot(nice_mds)
stressplot(nice_mds)

#Repeat to include info for both timepoint and location
t1_r <- dcast(subset(nice_melt, timepoint=="1" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t1_k <- dcast(subset(nice_melt, timepoint=="1" & location=="Karama"), timepoint+primer+location~sample, value="value") %>%
  mutate('246'="NA")
t1_k <- t1_k[,c(1:20,31,21:30)]

t2_r<-dcast(subset(nice_melt, timepoint=="2" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t2_k <- dcast(subset(nice_melt, timepoint=="2" & location=="Karama"), timepoint+primer+location~sample, value="value")

t3_r <- dcast(subset(nice_melt, timepoint=="3" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t3_k <- dcast(subset(nice_melt, timepoint=="3" & location=="Karama"), timepoint+primer+location~sample, value="value")

t4_r <- dcast(subset(nice_melt, timepoint=="4" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t4_k <- dcast(subset(nice_melt, timepoint=="4" & location=="Karama"), timepoint+primer+location~sample, value="value")

t5_r <- dcast(subset(nice_melt, timepoint=="5" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t5_k <- dcast(subset(nice_melt, timepoint=="5" & location=="Karama"), timepoint+primer+location~sample, value="value")

t6_r <- dcast(subset(nice_melt, timepoint=="6" & location=="Rubona"), timepoint+primer+location~sample, value="value")
t6_k <- dcast(subset(nice_melt, timepoint=="6" & location=="Karama"), timepoint+primer+location~sample, value="value")

colnames(t1_r) <- colnames(t1_k)
colnames(t2_r) <- colnames(t1_k)
colnames(t2_k) <- colnames(t1_k)
colnames(t3_r)<- colnames(t1_k)
colnames(t3_k)<- colnames(t1_k)
colnames(t4_r)<- colnames(t1_k)
colnames(t4_k)<- colnames(t1_k)
colnames(t5_r)<- colnames(t1_k)
colnames(t5_k)<- colnames(t1_k)
colnames(t6_r)<- colnames(t1_k)
colnames(t6_k)<- colnames(t1_k)

loc_time_cocur<-data.frame(rbind(t1_r, t1_k, t2_r, t2_k, t3_r, t3_k, 
                                 t4_r, t4_k, t5_r, t5_k, t6_r, t6_k))

#making the spearman's distance matrix
loc_time_cocur_dist<-spearman.dist(data.matrix(loc_time_cocur[,-c(1:3)]))

#PERMANOVA
set.seed(222)
adonis(loc_time_cocur_dist~timepoint*location,data=loc_time_cocur, permutations=9999) 
#loc and timepoint both highly sig (p=1e-4)
#suggests that locations should be analyzed separately for each timepoint

#NMDS
set.seed(213)
loc_time_mds<-monoMDS(loc_time_cocur_dist, model="local")
plot(loc_time_mds)
stressplot(loc_time_mds)
