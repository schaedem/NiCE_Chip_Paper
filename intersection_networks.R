#create graphs by loc and timepoint 
#create dfs for each locxtime with node degree, timepoint, loc, and betweenness centrality
library(tidyverse)
library(igraph)
source("correlation_dataframes.R")
source("network_graphs.R")

#Networks by loc x time
k_1_graph <- k_1_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_1_graph <- r_1_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

k_2_graph <- k_2_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_2_graph <- r_2_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

k_3_graph <- k_3_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_3_graph <- r_3_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

k_4_graph <- k_4_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_4_graph <- r_4_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

k_5_graph <- k_5_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_5_graph <- r_5_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

k_6_graph <- k_6_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

r_6_graph <- r_6_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

#intersection graphs within same timepoint

graph_t_1_int <- graph.intersection(k_1_graph, r_1_graph, keep.all.vertices=FALSE) %>% simplify()
graph_t_2_int <- graph.intersection(k_2_graph, r_2_graph, keep.all.vertices=FALSE) %>% simplify()
graph_t_3_int <- graph.intersection(k_3_graph, r_3_graph,keep.all.vertices=FALSE)%>% simplify()
graph_t_4_int <- graph.intersection(k_4_graph, r_4_graph,keep.all.vertices=FALSE)%>% simplify()
graph_t_5_int <- graph.intersection(k_5_graph, r_5_graph,keep.all.vertices=FALSE)%>% simplify()
graph_t_6_int <- graph.intersection(k_6_graph, r_6_graph,keep.all.vertices=FALSE)

#get edgelists for intersection graphs
edgelist_t1_int <- get.edgelist(graph_t_1_int)
edgelist_t2_int <- get.edgelist(graph_t_2_int)
edgelist_t3_int <- get.edgelist(graph_t_3_int)
edgelist_t4_int <- get.edgelist(graph_t_4_int)
edgelist_t5_int <- get.edgelist(graph_t_5_int)
edgelist_t6_int <- get.edgelist(graph_t_6_int)

#intersection for all timepoints x locations
graph_t_all_int<- graph.intersection(graph_t_1_int, graph_t_2_int, graph_t_3_int,
                                graph_t_4_int, graph_t_5_int, graph_t_6_int)
edgelist_t_all_int <- get.edgelist(graph_t_all_int)
#graph_t_all_int <- graph_from_edgelist(edgelist_t_all_int, directed=FALSE)

#anthesis intersection graph
graph_anth_int <- graph.intersection(graph_t_1_int, graph_t_3_int, graph_t_5_int)
#edgelist_anth_int <- get.edgelist(anthesis_intersect)
#graph_anth_int <- graph_from_edgelist(edgelist_anth_int, directed=FALSE)

#regrowth intersection graph
graph_regrow_int <- graph.intersection(graph_t_2_int, graph_t_4_int, graph_t_6_int)
#edgelist_regrow_int <- get.edgelist(regrowth_intersect)
#graph_regrow_int <- graph_from_edgelist(edgelist_regrow_int, directed=FALSE)

#make node list for graphing
make_node_list <- function(intersect_df, assay_list) {
  nodes_df <- make_nodes_df(intersect_df, assay_list)%>%
    mutate(process = factor(process, levels = c("ammonia_oxidation", "hydroxylamine_oxidation", "nitrite_oxidation",
                                                "nitrite_reduction", "nitric_oxide_reduction", "nitrous_oxide_reduction", "n_fixation")),
           pathway = factor(pathway, levels=c("nitrification", "denitrification", "comammox", "n_fixation")),
           pathway = ordered(pathway)) %>%
    as_tibble()  %>%
    select(-id) %>%
    filter(degree > 0)
  
  return(nodes_df)
}

nodes_t_1_int <- make_node_list(graph_t_1_int, assay_list)
nodes_t_2_int <- make_node_list(graph_t_2_int, assay_list)
nodes_t_3_int <- make_node_list(graph_t_3_int, assay_list) 
nodes_t_4_int <- make_node_list(graph_t_4_int, assay_list) 
nodes_t_5_int <- make_node_list(graph_t_5_int, assay_list) 
nodes_t_6_int <- make_node_list(graph_t_6_int, assay_list) 
nodes_t_all_int <- make_node_list(graph_t_all_int, assay_list)
nodes_anth_int <- make_node_list(graph_anth_int, assay_list)
nodes_regrow_int <- make_node_list(graph_regrow_int, assay_list)

#Plot networks
t_1_int_graph <- make_network_graph(nodes_t_1_int, edgelist_t1_int) + ggtitle("A")
t_2_int_graph <- make_network_graph(nodes_t_2_int, edgelist_t2_int) + ggtitle("B")
t_3_int_graph <- make_network_graph(nodes_t_3_int, edgelist_t3_int) + ggtitle("C")
t_4_int_graph <- make_network_graph(nodes_t_4_int, edgelist_t4_int) + ggtitle("D")
t_5_int_graph <- make_network_graph(nodes_t_5_int, edgelist_t5_int) + ggtitle("E")
t_6_int_graph <- make_network_graph(nodes_t_6_int, edgelist_t6_int) + ggtitle("F")
t_all_int_graph <- make_network_graph(nodes_t_all_int, edgelist_t_all_int)

#Arrange network graphs in one plot

plot.list <- lapply(list(t_1_int_graph,t_3_int_graph, t_5_int_graph, t_2_int_graph,t_4_int_graph,t_6_int_graph),
                    function(p) p + theme(plot.background = element_rect(color = "black")))

ggarrange(plotlist = plot.list)
ggsave(width = 19, height=11, filename="int_graphs.tiff", dpi=300)
#w=15; h=11
####Network-level statistics####

#Node number - didn't use V() because it included nodes with degree 0
nodes1 <- length(unique(nodes_t_1_int$primer_pair))
nodes2 <- length(unique(nodes_t_2_int$primer_pair))
nodes3 <- length(unique(nodes_t_3_int$primer_pair))
nodes4 <- length(unique(nodes_t_4_int$primer_pair))
nodes5 <- length(unique(nodes_t_5_int$primer_pair))
nodes6 <- length(unique(nodes_t_6_int$primer_pair))

#Edge number
edgenum1 <- length(E(graph_t_1_int))
edgenum2 <- length(E(graph_t_2_int))
edgenum3 <- length(E(graph_t_3_int))
edgenum4 <- length(E(graph_t_4_int))
edgenum5 <- length(E(graph_t_5_int))
edgenum6 <- length(E(graph_t_6_int))

#Modules/Modularity
mod_t_1_int <- edge.betweenness.community(graph_t_1_int, directed=FALSE) #mod=0.00028
mod_t_2_int <- edge.betweenness.community(graph_t_2_int, directed=FALSE) #mod=0.19
mod_t_3_int <- edge.betweenness.community(graph_t_3_int, directed=FALSE) #mod=0.18
mod_t_4_int <- edge.betweenness.community(graph_t_4_int, directed=FALSE) #mod=0.019
mod_t_5_int <- edge.betweenness.community(graph_t_5_int, directed=FALSE) #mod=0.24
mod_t_6_int <- edge.betweenness.community(graph_t_6_int, directed=FALSE) #mod=0

mod_all_int <- data.frame(modularity = c(0, 0.19, 0.18, 0.02, 0.24, 0), timepoint = c(1:6),
                            growth = c("anth", "regrow", "anth", "regrow", "anth", "regrow"))

#Transitivity
trans_t_1_int <- transitivity(graph_t_1_int) #0.8851894
trans_t_2_int <- transitivity(graph_t_2_int) #0.6157895
trans_t_3_int <- transitivity(graph_t_3_int) #0.7317073
trans_t_4_int <- transitivity(graph_t_4_int) #0.6830357
trans_t_5_int <- transitivity(graph_t_5_int) # 0.4636364
trans_t_6_int <- transitivity(graph_t_6_int) #0.375

trans_all_int <- data.frame(transitivity = c(trans_t_1_int, trans_t_2_int, trans_t_3_int,
                       trans_t_4_int, trans_t_5_int, trans_t_6_int), timepoint = c(1:6),
                       growth = c("anth", "regrow", "anth", "regrow", "anth", "regrow"))

#Density
dens_t_1_int <- edge_density(graph_t_1_int) #0.83
dens_t_2_int <- edge_density(graph_t_2_int) #0.28
dens_t_3_int <- edge_density(graph_t_3_int) #0.46
dens_t_4_int <- edge_density(graph_t_4_int) #0.42
dens_t_5_int <- edge_density(graph_t_5_int) #0.25
dens_t_6_int <- edge_density(graph_t_6_int) #0.38

#testing network differences for plant growth stage
anth <- c(dens_t_1_int, dens_t_3_int, dens_t_5_int)
regr <- c(dens_t_2_int, dens_t_4_int, dens_t_6_int)
#mean
mean(anth)
mean(regr)
#SE
sd(anth)/sqrt(3)
sd(regr)/sqrt(3)

#Average path length
path_t_1_int <- average.path.length(graph_t_1_int)
path_t_2_int <- average.path.length(graph_t_2_int)
path_t_3_int <- average.path.length(graph_t_3_int)
path_t_4_int <- average.path.length(graph_t_4_int)
path_t_5_int <- average.path.length(graph_t_5_int)
path_t_6_int <- average.path.length(graph_t_6_int)
path_all_int <- data.frame(path_length = c(path_t_1_int, path_t_2_int, path_t_3_int,
                                            path_t_4_int, path_t_5_int, path_t_6_int), timepoint = c(1:6))

dens_all_int <- data.frame(density = c(dens_t_1_int, dens_t_2_int, dens_t_3_int,
                                             dens_t_4_int, dens_t_5_int, dens_t_6_int), timepoint = c(1:6),
                            growth = c("anth", "regrow", "anth", "regrow", "anth", "regrow"))

graph_stats <- full_join(dens_all_int, trans_all_int, by=c("growth", "timepoint")) %>%
  merge(mod_all_int) %>%
  merge(path_all_int)

ggplot(graph_stats, aes(x=transitivity, y=density)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method="lm", se=FALSE)

ggplot(graph_stats, aes(x=transitivity, y=modularity)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method="lm", se=FALSE)

ggplot(graph_stats, aes(x=density, y=modularity)) +
  geom_point() +
  theme_classic() +
  geom_smooth(method="lm")

ggplot(graph_stats, aes(x=timepoint, y=modularity)) +
  geom_point() +
  theme_classic() +
  geom_line()

ggplot(graph_stats, aes(x=timepoint, y=density)) +
  geom_point() +
  theme_classic() +
  geom_line()

ggplot(graph_stats, aes(x=timepoint, y=transitivity)) +
  geom_point() +
  geom_line() +
  theme_classic() 


#Small world vs free-scale

deg_dist1 <- degree.distribution(graph_t_1_int)
plot(deg_dist1,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph1",log="xy")

deg_dist2 <- degree.distribution(graph_t_2_int)
plot(deg_dist2,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph2",log="xy")

deg_dist3 <- degree.distribution(graph_t_3_int)
plot(deg_dist3,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph3",log="xy")

deg_dist4 <- degree.distribution(graph_t_4_int)
plot(deg_dist4,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph4",log="xy")

deg_dist5 <- degree.distribution(graph_t_5_int)
plot(deg_dist5,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph5",log="xy")

deg_dist6 <- degree.distribution(graph_t_6_int)
plot(deg_dist6,pch=20,xlab="k",ylab="P(k)",las=1,main="Graph6",log="xy")



#anth_int_graph <- make_network_graph(nodes_anth_int, edgelist_anth_int)
#regrow_int_graph <- make_network_graph(nodes_regrow_int, edgelist_regrow_int)

#Save graphs
#ggsave("graph1_int.tiff", plot=t_1_int_graph, height=8, width=8)
#ggsave("graph2_int.tiff", plot=t_2_int_graph, height=8, width=8)
#ggsave("graph3_int.tiff", plot=t_3_int_graph, height=8, width=8)
#ggsave("graph4_int.tiff", plot=t_4_int_graph, height=8, width=8)
#ggsave("graph5_int.tiff", plot=t_5_int_graph, height=8, width=8)
#ggsave("graph6_int.tiff", plot=t_6_int_graph, height=8, width=8)
#ggsave(plot= t_all_int_graph, "t_loc_consensus_network.tiff", dpi=300, height=4, width=4)

#ggsave(plot= anth_int_graph, "anthesis_consensus_network.tiff", dpi=300, height=6, width=6)
#ggsave(plot= regrow_int_graph, "regrowth_consensus_network.tiff", dpi=300, height=6, width=6)

