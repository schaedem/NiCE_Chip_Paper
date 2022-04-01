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

graph_t_1_int <- graph.intersection(k_1_graph, r_1_graph) %>% simplify()
graph_t_2_int <- graph.intersection(k_2_graph, r_2_graph) %>% simplify()
graph_t_3_int <- graph.intersection(k_3_graph, r_3_graph)%>% simplify()
graph_t_4_int <- graph.intersection(k_4_graph, r_4_graph)%>% simplify()
graph_t_5_int <- graph.intersection(k_5_graph, r_5_graph)%>% simplify()
graph_t_6_int <- graph.intersection(k_6_graph, r_6_graph)%>% simplify()

#get edgelists for intersection graphs
edgelist_t1_int <- get.edgelist(graph_t_1_int)
edgelist_t2_int <- get.edgelist(graph_t_2_int)
edgelist_t3_int <- get.edgelist(graph_t_3_int)
edgelist_t4_int <- get.edgelist(graph_t_4_int)
edgelist_t5_int <- get.edgelist(graph_t_5_int)
edgelist_t6_int <- get.edgelist(graph_t_6_int)

#make graph object
#graph_t_1_int <- graph_from_edgelist(edgelist_t1_int, directed=FALSE)
#graph_t_2_int <- graph_from_edgelist(edgelist_t2_int, directed=FALSE)
#graph_t_3_int <- graph_from_edgelist(edgelist_t3_int, directed=FALSE)
#graph_t_4_int <- graph_from_edgelist(edgelist_t4_int, directed=FALSE)
#graph_t_5_int <- graph_from_edgelist(edgelist_t5_int, directed=FALSE)
#graph_t_6_int <- graph_from_edgelist(edgelist_t6_int, directed=FALSE)

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
t_1_int_graph <- make_network_graph(nodes_t_1_int, edgelist_t1_int)
t_2_int_graph <- make_network_graph(nodes_t_2_int, edgelist_t2_int)
t_3_int_graph <- make_network_graph(nodes_t_3_int, edgelist_t3_int)
t_4_int_graph <- make_network_graph(nodes_t_4_int, edgelist_t4_int)
t_5_int_graph <- make_network_graph(nodes_t_5_int, edgelist_t5_int)
t_6_int_graph <- make_network_graph(nodes_t_6_int, edgelist_t6_int)
t_all_int_graph <- make_network_graph(nodes_t_all_int, edgelist_t_all_int)

####Network-level statistics####

#Modules/Modularity
mod_t_1_int <- edge.betweenness.community(graph_t_1_int, directed=FALSE) #mod=0.00028
mod_t_2_int <- cluster_edge_betweenness(graph_t_2_int, directed=FALSE) #mod=0.19
mod_t_3_int <- edge.betweenness.community(graph_t_3_int, directed=FALSE) #mod=0.18
mod_t_4_int <- edge.betweenness.community(graph_t_4_int, directed=FALSE) #mod=0.019
mod_t_5_int <- edge.betweenness.community(graph_t_5_int, directed=FALSE) #mod=0.24
mod_t_6_int <- edge.betweenness.community(graph_t_6_int, directed=FALSE) #mod=0


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

ggplot(trans_all_int, aes(x=timepoint, y=transitivity, color=growth)) +
  geom_point() +
  geom_line() +
  theme_classic()

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

