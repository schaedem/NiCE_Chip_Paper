library(tidyverse)
library(igraph)
source('intersection_networks.R')

between_df <- function(graph) {
  between_mat <- betweenness(graph) %>%
    data.frame() %>%
    rownames_to_column("primer_pair") %>%
    rename(betweenness = '.') %>%
    merge(assay_list, by="primer_pair")
}

between_t_1_int <- between_df(graph_t_1_int) %>% mutate(timepoint="1", date=lubridate::mdy("9/16/2020"))
between_t_2_int <- between_df(graph_t_2_int) %>% mutate(timepoint="2", date=lubridate::mdy("10/06/2020"))
between_t_3_int <- between_df(graph_t_3_int) %>% mutate(timepoint= "3", date=lubridate::mdy("11/18/2020"))
between_t_4_int <- between_df(graph_t_4_int) %>% mutate(timepoint="4", date=lubridate::mdy("12/07/2020"))
between_t_5_int <- between_df(graph_t_5_int) %>% mutate(timepoint="5", date=lubridate::mdy("01/19/2021"))
between_t_6_int <- between_df(graph_t_6_int) %>% mutate(timepoint="6", date=lubridate::mdy("02/10/2021"))

between_boxplot <- function(between_graph) {
 graph <- between_graph
  plot <- 
   ggplot(data=graph, aes(x=reorder(primer_pair, betweenness), y=betweenness, 
                                 fill=primer_pair)) +
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
  filter(pathway== "nitrification" | pathway == "denitrification") %>%
  mutate(process = factor(process, levels =c("ammonia_oxidation", "hydroxylamine_oxidation",
                                             "nitrite_oxidation", "nitrite_reduction", "nitric_oxide_reduction", 
                                             "nitrous_oxide_reduction")),
         primer_pair = factor(primer_pair, levels=c("amoa_F1 / amoA_2R", "Arch-amoAFA / Arch-amoAR", "Arch-amoAFB / Arch-amoAR",
                                                    "Gamo172_F1 / Gamo172_F1_R", "Gamo172_F1 / Gamo172_F1_R2", "Gamo172_F2 / Gamo172_F2_R1",
                                                    "haoF4 / haoR2", "hzocl1F1 / hzocl1R2", "NxrB169F / NxrB638R", "NxrB1F / NxrB1R",
                                                    "FlaCu / R3Cu", "nirK876 / nirK1040", "nirKfF / nirKfR", "norB2 / norB6",
                                                    "qnorB2F / qnorB5R", "nosZ1F / nosZ1R", "NosZ912F / NosZ853R")),
         pathway = factor(pathway, levels=c("nitrification", "denitrification"))) %>%
  filter(!is.na(primer_pair))
levels(all_between$pathway) <- c("Nitrification", "Denitrification")

pal_df <- all_between %>% select(colorcode, process, primer_pair) %>% unique()
pal <-c("#AD2831", "#640D14","#640D14", "#F38375", "#EF6351", "#EF6351", "#FBC3BC",
        "#F7A399", "#D00000", "#DC2F02", "#01497C", "#01497C", "#80CDC1",
        "#012A4A", "#012A4A", "#61A5C2", "#89C2D9")
#need to correctly align primer pairs with their colors!
betweenness_time_plot <-
  ggplot(data=all_between, aes(x=date, y=betweenness, color=primer_pair, 
                               group=primer_pair, reorder(process, primer_pair))) +
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

#Modules/Modularity
mod_t_1_int <- edge.betweenness.community(graph_t_1_int, directed=FALSE) #mod=0.00028; 4 groups
mod_t_2_int <- edge.betweenness.community(graph_t_2_int, directed=FALSE) #mod=0.19; 4 groups
mod_t_3_int <- edge.betweenness.community(graph_t_3_int, directed=FALSE) #mod=0.18; 3 groups
mod_t_4_int <- edge.betweenness.community(graph_t_4_int, directed=FALSE) #mod=0.019; 5 groups
mod_t_5_int <- edge.betweenness.community(graph_t_5_int, directed=FALSE) #mod=0.24; 4 groups
mod_t_6_int <- edge.betweenness.community(graph_t_6_int, directed=FALSE) #mod=0; 1 group

#Degree vs betweenness
deg_betw_t_1_int <- nodes_t_1_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)

deg_betw_t_2_int <- nodes_t_2_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)

deg_betw_t_3_int <- nodes_t_3_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)

deg_betw_t_4_int <- nodes_t_4_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)

deg_betw_t_5_int <- nodes_t_5_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)

deg_betw_t_6_int <- nodes_t_6_int %>% 
  select(primer_pair, degrees) %>%
  merge(between_t_1_int)


ggplot(deg_betw_t_1_int, aes(x=degrees, y=betweenness)) +
  geom_point() +
  theme_classic() +
  geom_smooth(se=FALSE, method='lm') +
  ggtitle("1")
