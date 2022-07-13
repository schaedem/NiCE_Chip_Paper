library(tidyverse)
library(igraph)
source('intersection_networks.R')

degree_df <- function(graph) {
  degree_mat <- degree(graph) %>%
    data.frame() %>%
    rownames_to_column("acronym") %>%
    rename(degree = '.') %>%
    merge(assay_list, by="acronym")
}

between_df <- function(graph) {
  between_mat <- betweenness(graph) %>%
    data.frame() %>%
    rownames_to_column("acronym") %>%
    rename(betweenness = '.') %>%
    merge(assay_list, by="acronym")
}

closeness_df <- function(graph) {
  closeness_mat <- closeness(graph) %>%
    data.frame() %>%
    rownames_to_column("acronym") %>%
    rename(closeness_centrality = '.')
  
}

#make df with node degree
degree1 <- degree_df(graph_t_1_int) %>% mutate(timepoint="1")
degree2 <- degree_df(graph_t_2_int) %>% mutate(timepoint="2")
degree3 <- degree_df(graph_t_3_int)%>% mutate(timepoint="3")
degree4 <- degree_df(graph_t_4_int)%>% mutate(timepoint="4")
degree5 <- degree_df(graph_t_5_int)%>% mutate(timepoint="5")
degree6 <- degree_df(graph_t_6_int)%>% mutate(timepoint="6")
degree_all <- rbind(degree1, degree2, degree3, degree4, degree5, degree6)

#make df with closeness scores
close1 <- closeness_df(graph_t_1_int) %>% mutate(timepoint="1")
close2 <- closeness_df(graph_t_2_int) %>% mutate(timepoint= "2")
close3 <- closeness_df(graph_t_3_int) %>% mutate(timepoint = "3")
close4 <- closeness_df(graph_t_4_int) %>% mutate(timepoint = "4")
close5 <- closeness_df(graph_t_5_int) %>% mutate(timepoint = "5")
close6 <- closeness_df(graph_t_6_int) %>% mutate(timepoint = "6")
close_all <- rbind(close1, close2, close3, close4, close5, close6)

#make df with betweenness scores and metadata for graphing
between_t_1_int <- between_df(graph_t_1_int) %>% mutate(timepoint="1", date=lubridate::mdy("9/16/2020"))
between_t_2_int <- between_df(graph_t_2_int) %>% mutate(timepoint="2", date=lubridate::mdy("10/06/2020"))
between_t_3_int <- between_df(graph_t_3_int) %>% mutate(timepoint= "3", date=lubridate::mdy("11/18/2020"))
between_t_4_int <- between_df(graph_t_4_int) %>% mutate(timepoint="4", date=lubridate::mdy("12/07/2020"))
between_t_5_int <- between_df(graph_t_5_int) %>% mutate(timepoint="5", date=lubridate::mdy("01/19/2021"))
between_t_6_int <- between_df(graph_t_6_int) %>% mutate(timepoint="6", date=lubridate::mdy("02/10/2021"))

between_boxplot <- function(between_graph) {
 graph <- between_graph
  plot <- 
   ggplot(data=graph, aes(x=reorder(acronym, betweenness), y=betweenness, 
                                 fill=acronym)) +
  geom_bar(stat="identity") +
  theme_classic() +
  coord_flip() +
  theme(legend.position = "none") +
  scale_fill_manual(values= graph$colorcode)+
  ylab("Betweenness") +
  xlab("")
 
 return(plot)

}

between_plot_1 <- between_boxplot(between_t_1_int)  + ggtitle("1")
between_plot_2 <- between_boxplot(between_t_2_int)  + ggtitle("2")
between_plot_3 <- between_boxplot(between_t_3_int)  + ggtitle("3")
between_plot_4 <- between_boxplot(between_t_4_int)  + ggtitle("4")
between_plot_5 <- between_boxplot(between_t_5_int)  + ggtitle("5")
between_plot_6 <- between_boxplot(between_t_6_int)  + ggtitle("6")


plot.list <- lapply(list(between_plot_1, between_plot_2, between_plot_3, between_plot_4,
                         between_plot_5, between_plot_6), 
                    function(p) p + theme(plot.background = element_rect(color = "white")))
ggarrange(plotlist = plot.list)

#alternative to faceted box plots: time/line plot
all_between <- rbind(between_t_1_int, between_t_2_int, between_t_3_int,
                     between_t_4_int, between_t_5_int, between_t_6_int) %>%
  full_join(close_all, by=c('acronym', 'timepoint')) %>%
  full_join(degree_all, by=c('acronym', 'timepoint', 'gene_org', 'pathway', 'process')) %>%
  filter(pathway== "nitrification" | pathway == "denitrification") %>%
  mutate(process = factor(process, levels =c("ammonia_oxidation", "hydroxylamine_oxidation",
                                             "nitrite_oxidation", "nitrite_reduction", "nitric_oxide_reduction", 
                                             "nitrous_oxide_reduction")),
       #  primer_pair = factor(primer_pair, levels=c("amoa_F1 / amoA_2R", "Arch-amoAFA / Arch-amoAR", "Arch-amoAFB / Arch-amoAR",
        #                                            "Gamo172_F1 / Gamo172_F1_R", "Gamo172_F1 / Gamo172_F1_R2", "Gamo172_F2 / Gamo172_F2_R1",
         #                                           "haoF4 / haoR2", "hzocl1F1 / hzocl1R2", "NxrB169F / NxrB638R", "NxrB1F / NxrB1R",
          #                                          "FlaCu / R3Cu", "nirK876 / nirK1040", "nirKfF / nirKfR", "norB2 / norB6",
           #                                         "qnorB2F / qnorB5R", "nosZ1F / nosZ1R", "NosZ912F / NosZ853R")),
         acronym = factor(acronym, levels = c("Beta_amoA", "Arch_amoA_FA", "Arch_amoA_FB", "Gamma_amoA_F1R1", "Gamma_amoA_F1R2", "Gamma_amoA_F2R1",
                                              "Proteo_hao", "Annamox_hzocl", "Nitrospira_nxrB", "Nitrobacter_nxrB", "nirK_FlaCu", "nirK_876",
                                              "nirK_Fungi", "norB_2", "qnorB_2F-5R", "nosZ_1F", "nosZ_912F")),
         pathway = factor(pathway, levels=c("nitrification", "denitrification"))) %>%
  filter(!is.na(acronym)) 
levels(all_between$pathway) <- c("Nitrification", "Denitrification")

pal <-c("#AD2831", "#640D14","#640D14", "#F38375", "#EF6351", "#EF6351", "#FBC3BC",
        "#F7A399", "#D00000", "#DC2F02", "#01497C", "#01497C", "#80CDC1",
        "#012A4A", "#012A4A", "#61A5C2", "#89C2D9")

#need to correctly align primer pairs with their colors!
betweenness_time_plot <-
  ggplot(data=all_between, aes(x=date, y=betweenness, color=acronym, 
                               group=acronym, reorder(process, acronym))) +
  geom_point(size=3) +
  geom_path() +
  theme_classic() +
  scale_color_manual(values=pal) +
  facet_wrap(~pathway) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  ylab("Betweenness") +
  xlab("") 
betweenness_time_plot
ggsave("betweenness_plot.tiff", plot=betweenness_time_plot, height=6, width=10)


closeness_time_plot <-
  ggplot(data=all_between, aes(x=date, y=closeness_centrality, color=acronym, 
                               group=acronym, reorder(process, acronym))) +
  geom_point(size=3) +
  geom_path() +
  theme_classic() +
  scale_color_manual(values=pal) +
  facet_wrap(~pathway) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  ylab("Closeness Centrality") +
  xlab("") 
closeness_time_plot

head(all_between)

#find keystones based on betweenness, centrality, degree
keystone_df <- all_between %>%
  select(timepoint, acronym, closeness_centrality, betweenness, pathway, degree)

sum_df <- keystone_df %>%
  group_by(pathway) %>%
  summarise(mean(betweenness), (sd(betweenness)/length(pathway)))

keystone_nit_between <- keystone_df %>%
  filter(pathway=="Nitrification") %>%
  group_by(timepoint) %>%
  top_n(1, betweenness)
 
keystone_nit_close <- keystone_df %>%
  filter(pathway=="Nitrification") %>%
  group_by(timepoint) %>%
  top_n(1, closeness_centrality)

keystone_nit_degree <- keystone_df %>%
  filter(pathway=="Nitrification") %>%
  group_by(timepoint) %>%
  top_n(1, degree)

keystone_denit_between <- keystone_df %>%
  filter(pathway=="Denitrification") %>%
  group_by(timepoint) %>%
  top_n(1, betweenness)

keystone_denit_close <- keystone_df %>%
  filter(pathway=="Denitrification") %>%
  group_by(timepoint) %>%
  top_n(1, closeness_centrality)

keystone_denit_degree <- keystone_df %>%
  filter(pathway=="Denitrification") %>%
  group_by(timepoint) %>%
  top_n(1, degree)

keystone_nit <- keystone_nit_between %>%
  rbind(keystone_nit_close)%>%
  rbind(keystone_nit_degree) %>%
  unique()
keystone_denit <- keystone_denit_between %>%
  rbind(keystone_denit_close) %>%
  rbind(keystone_denit_degree) %>%
  unique()

all_keystone <- keystone_nit %>%
  full_join(keystone_denit)

key_primers <- all_keystone[,2] %>%
  unique()

all_keystone_time <- all_between %>%
  semi_join(key_primers) %>%
  as_tibble() %>%
  mutate(acronym = droplevels(acronym),
         process=droplevels(process),
         pathway=droplevels(pathway))
levels(all_keystone_time$pathway)

#simplified graphs using only potential keystone genes
pal2 <- all_keystone_time %>%
  select(acronym, colorcode.x) %>%
  distinct(acronym, .keep_all=TRUE)%>%
  arrange(acronym)

pal3 <- pal2[[2]]

keystone_closeness_time_plot <-
  ggplot(data=all_keystone_time, aes(x=date, y=closeness_centrality, color=acronym, 
                               group=acronym, reorder(process,acronym))) +
  geom_point(size=3) +
  geom_line() +
  theme_classic() +
  scale_color_manual(values=pal3) +
  facet_wrap(~pathway) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  ylab("Closeness Centrality") +
  xlab("") 
keystone_closeness_time_plot 

ggsave("closeness_plot_simplified.tiff", plot=keystone_closeness_time_plot, height=6, width=10)


keystone_betweenness_time_plot <-
  ggplot(data=all_keystone_time, aes(x=date, y=betweenness, color=acronym, 
                                     group=acronym, reorder(process, acronym))) +
  geom_point(size=3) +
  geom_line(na.rm=T) +
  theme_classic() +
  scale_color_manual(values=pal3) +
  facet_wrap(~pathway) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  ylab("Betweenness Centrality") +
  xlab("") 
keystone_betweenness_time_plot
ggsave("betweenness_plot_simplified.tiff", plot=keystone_betweenness_time_plot, height=6, width=10)


#	Gamo172_F1 / Gamo172_F1_R2 & haoF4 / haoR2
#t5: nirK876 / nirK1040
#t2:norB2 / norB6

summ <- all_keystone_time %>%
  group_by(pathway) %>%
  summarise(mean(betweenness))

keystone_degree_time_plot <-
  ggplot(data=all_keystone_time, aes(x=date, y=degree, color=acronym, 
                                     group=acronym, reorder(acronym))) +
  geom_point(size=3) +
  geom_line(na.rm=TRUE) +
  theme_classic() +
  scale_color_manual(values=pal3) +
  facet_wrap(~pathway) +
  theme(legend.title = element_blank(), legend.position="bottom") +
  ylab("Normalized node degree") +
  xlab("") 
keystone_degree_time_plot
ggsave("degree_plot_simplified.tiff", plot=keystone_degree_time_plot, height=6, width=10)

