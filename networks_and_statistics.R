library(tidyverse)
library(igraph)
source("intersection_networks.R")
source("network_graphs.R")

#Rubona node attribute df
r_1_nodes <- make_nodes_df(r_1_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=1, location="Rubona") 

r_2_nodes <- make_nodes_df(r_2_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=2, location="Rubona") 

r_3_nodes <- make_nodes_df(r_3_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=3, location="Rubona") 

r_4_nodes <- make_nodes_df(r_4_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=4, location="Rubona") 

r_5_nodes <- make_nodes_df(r_5_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=5, location="Rubona") 

r_6_nodes <- make_nodes_df(r_6_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint=6, location="Rubona") 

#Karama node attribute df
k_1_nodes <- make_nodes_df(k_1_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 1, location = "Karama")

k_2_nodes <- make_nodes_df(k_2_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 2, location = "Karama")

k_3_nodes <- make_nodes_df(k_3_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 3, location = "Karama")

k_4_nodes <- make_nodes_df(k_4_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 4, location = "Karama")

k_5_nodes <- make_nodes_df(k_5_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 5, location = "Karama")

k_6_nodes <- make_nodes_df(k_6_graph, assay_list) %>%
  select(primer_pair, degree, betweenness, closeness_centrality, pathway) %>%
  mutate(timepoint = 6, location = "Karama")

#bind all node data together
all_nodes <- rbind(r_1_nodes, r_2_nodes, r_3_nodes, r_4_nodes, r_5_nodes, r_6_nodes,
                   k_1_nodes, k_2_nodes, k_3_nodes, k_4_nodes, k_5_nodes, k_6_nodes)

#rename pathway levels for plotting
all_nodes$pathway <- recode(all_nodes$pathway, comammox = 'Comammox', 
                        denitrification = 'Denitrification',
                        n_fixation = 'N Fixation', nitrification = 'Nitrification')
all_nodes <- all_nodes %>%
  mutate(growth = case_when(timepoint==1 ~ "anthesis",
                            timepoint==2 ~ "regrowth",
                            timepoint==3 ~"anthesis",
                            timepoint==4 ~ "regrowth",
                            timepoint==5 ~ "anthesis", 
                            timepoint ==6 ~ "regrowth"))

#Investigating relationships between network statistics - linear relationships?
summ <- all_nodes %>%
  group_by(timepoint, growth) %>%
  mutate(mean_degree = mean(degree)) %>%
  ungroup()

ggplot(summ, aes(x=timepoint, y=mean_degree, color=growth)) +
  geom_point() + 
  geom_line() +
  theme_classic() +
  ylab("Mean Degree") +
  theme(legend.title=element_blank())

degree_closeness <-
  ggplot(all_nodes, aes(x=degree, y=closeness_centrality, color=pathway)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("#C18C5D", "#3A68AE", "#66A61E", "#C03830")) +
 facet_wrap(~pathway) +
  theme(legend.position = "none") +
  geom_smooth(method="lm") +
  ylab("Closeness Centrality") +
  xlab("Normalized Node Degree")

degree_between <- 
  ggplot(all_nodes, aes(x=degree, y=betweenness, color=pathway)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("#C18C5D", "#3A68AE", "#66A61E", "#C03830")) +
  facet_wrap(~pathway) +
  theme(legend.position = "none") +
  geom_smooth(method="lm") +
  ylab("Betweenness Centrality") +
  xlab("Normalized Node Degree")

close_between <- 
  ggplot(all_nodes, aes(x=closeness_centrality, y=betweenness, color=pathway)) +
  geom_point() +
  theme_classic() +
  scale_color_manual(values=c("#C18C5D", "#3A68AE", "#66A61E", "#C03830")) +
  facet_wrap(~pathway) +
  theme(legend.position="none") +
  geom_smooth(method="lm") +
  ylab("Betweenness Centrality") +
  xlab("Closeness Centrality")

plot.list <- list(degree_closeness,degree_between, close_between)
      
plot_arrange <- ggarrange(plotlist = plot.list)

ggsave(plot=plot_arrange, "network_statistics.tiff", dpi=300, width=10, height=10)


