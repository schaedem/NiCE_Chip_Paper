library(tidyverse)
library(readxl)
library(ggtext)
library(viridis)
source("cooccurrence_permanova.R")
d<-SNCAnalysis::full_dat

#extract species scores from monoMDS
data.scores <- scores(nice_mds) %>%
  as_tibble(rownames="sample") %>%
  mutate(timepoint = nice_cocur$timepoint,
         season = ifelse(timepoint == "1", "dry", ifelse(timepoint=="2", "dry", "rainy")))

#create tibble for plotting centroids
centroids <- data.scores %>%
  select(MDS1, MDS2, timepoint) %>%
  group_by(timepoint) %>%
  summarise(MDS1 = mean(MDS1),
            MDS2 = mean(MDS2)) 

#new tibble for creating a star plot
starplot_meta <- data.scores %>%
  select(MDS1, MDS2, timepoint) %>%
  group_by(timepoint) %>%
  mutate(centroid1 = mean(MDS1), centroid2=mean(MDS2)) %>%
  ungroup()

#adding timepoint to MDS scores
data.scores$timepoint = nice_cocur$timepoint

#en_coord_cont = as.data.frame(scores(en, "vectors")) * ordiArrowMul(en)
#en_coord_cat = as.data.frame(scores(en, "factors")) * ordiArrowMul(en)

#simple ordination plot with centroids overlaid
centroid_plot <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2)) + 
  geom_point(data = data.scores, aes(colour = as.factor(timepoint)), size = 2, alpha = 0.4) + 
  #scale_colour_manual(values = c("orange", "steelblue")) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30")) +
  labs(colour = "Timepoint") +
  geom_point(data=centroids, 
             mapping=aes(x=MDS1, y=MDS2, color=as.factor(timepoint)), shape=18, size=4.5, show.legend=FALSE)
centroid_plot

ggsave("timepoint_nmds.tiff", width=8, height=5, dpi=300)


#ordination plot with shaded ellipses
#stat_ellipse level controls confidence interval - lower confidence = smaller ellipses

my_legend <- tibble(x=c(0.6, 1.8, -1.2, 1.3, -1.2, 0.65),
                    y=c(-1.25, 0.9, 1.1, -0.5, -0.7, 1.4),
                    color = c("1", "2", "3", "4", "5", "6"),
                    label= c("1. late dry\nanthesis", "2. late dry\nregrowth", 
                             "3. mid rainy\nanthesis",
                             "4. mid rainy\nregrowth", "5. late rainy\nanthesis",
                             "6. late rainy\nregrowth"))

ellipse_plot <- ggplot(data = data.scores, aes(x = MDS1, y = MDS2, color=as.factor(timepoint))) + 
  stat_ellipse(level=0.8, geom="polygon", alpha=0.2,show.legend=FALSE, aes(fill=as.factor(timepoint)))+
  geom_point(size = 2, alpha = 0.4) + 
  theme(axis.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        panel.background = element_blank(), panel.border = element_rect(fill = NA, colour = "grey30"), 
        axis.ticks = element_blank(), axis.text = element_blank(), legend.key = element_blank(), 
        legend.title = element_text(size = 10, face = "bold", colour = "grey30"), 
        legend.text = element_text(size = 9, colour = "grey30"),
        legend.position = "none") +
  labs(colour = "Timepoint") +
  scale_color_viridis(discrete=TRUE, option="D")+
  scale_fill_viridis(discrete=TRUE, option="D") +
  theme_classic() +
  theme(legend.position="none") +
  geom_text(data=my_legend, 
          aes(x=x, y=y, color=color, label=label),
          inherit.aes=FALSE, show.legend=FALSE, fontface="bold", lineheight=0.8, hjust=1) 

ellipse_plot

ggsave("timepoint_nmds_ellipse_viridis.tiff", width=8, height=5, dpi=300)

##repeat for loc-time MDS ####

#extract species scores from monoMDS
data.scores.loc.time <- scores(loc_time_mds) %>%
  as_tibble(rownames="sample") %>%
  mutate(timepoint = loc_time_cocur$timepoint,
         location = loc_time_cocur$location)

#create tibble for plotting centroids
centroids.time<- data.scores.loc.time %>%
  select(MDS1, MDS2, timepoint, location) %>%
  group_by(timepoint) %>%
  summarise(MDS1_time = mean(MDS1),
            MDS2_time = mean(MDS2)) %>%
  ungroup()

#adding timepoint to MDS scores
data.scores.loc.time$timepoint = loc_time_cocur$timepoint

#simple ordination plot with centroids overlaid
centroid_plot <- ggplot() + 
  geom_point(data = data.scores.loc.time, 
             aes(x = MDS1, y = MDS2, color=location), size = 2, alpha = 0.4, stat="identity") + 
  scale_color_manual(values=c("#E41A1C", "#377EB8")) +
  theme_classic() +
  geom_point(data=centroids.time, 
             aes(x=MDS1_time, y=MDS2_time, fill=as.factor(timepoint)),
              size=4.5, inherit.aes=FALSE) +
  viridis::scale_fill_viridis(option="D", discrete=TRUE)
centroid_plot

RColorBrewer::brewer.pal(4,"Set1")

#ordination plot with shaded ellipses
#stat_ellipse level controls confidence interval - lower confidence = smaller ellipses

my_legend <- tibble(x=c(0.6, 1.8, -1.2, 1.3, -1.2, 0.65),
                    y=c(-1.25, 0.9, 1.1, -0.5, -0.7, 1.4),
                    color = c("1", "2", "3", "4", "5", "6"),
                    label= c("1. late dry\nanthesis", "2. late dry\nregrowth", 
                             "3. mid rainy\nanthesis",
                             "4. mid rainy\nregrowth", "5. late rainy\nanthesis",
                             "6. late rainy\nregrowth"))

ellipse_plot <- ggplot() + 
  geom_point(data = data.scores.loc.time, 
             aes(x = MDS1, y = MDS2, color=location), 
             size = 2, alpha = 0.8) + 
  scale_color_manual(values=c("#E41A1C", "#377EB8")) +
  theme(legend.key=element_rect(fill=NA, color=NA)) +
  theme_classic() +
  stat_ellipse(data=data.scores.loc.time,level=0.8, geom="polygon", alpha=0.3,
               show.legend=TRUE, aes(x=MDS1, y=MDS2,fill=as.factor(timepoint)))+
  labs(fill = "Timepoint", color="") +
  viridis::scale_fill_viridis(discrete=TRUE, option="D") 

  #geom_text(data=my_legend, 
           # aes(x=x, y=y, color=color, label=label),
           # inherit.aes=FALSE, show.legend=FALSE, fontface="bold", lineheight=0.8, hjust=1) 

ellipse_plot

ggsave("loc_time_nmds_ellipse_viridis.tiff", width=8, height=5, dpi=300)

