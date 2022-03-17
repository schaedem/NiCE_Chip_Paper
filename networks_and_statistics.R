library(tidyverse)
library(igraph)
library(Hmisc)
library(corrr)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")
source("correlation_dataframes.R")

#Networks for each timepoint
nice_1_graph <- nice_1_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE) 

nice_2_graph <- nice_2_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p < 0.01) %>%
  graph_from_data_frame(directed=FALSE)

nice_3_graph <- nice_3_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p<0.01) %>%
  graph_from_data_frame(directed=FALSE)

nice_4_graph <- nice_4_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p<0.01) %>%
  graph_from_data_frame(directed=FALSE)

nice_5_graph <- nice_5_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p<0.01) %>%
  graph_from_data_frame(directed=FALSE)

nice_6_graph <- nice_6_cor_df %>%
  filter(abs(r) > 0.75) %>%
  filter(p<0.01) %>%
  graph_from_data_frame(directed=FALSE)

#Network statistics for each network
#network statistics
#get edgelist
edgelist_1 = get.edgelist(nice_1_graph)
edgelist_2 = get.edgelist(nice_2_graph)
edgeList_3 = get.edgelist(nice_3_graph)
edgelist_4 = get.edgelist(nice_4_graph)
edgelist_5 = get.edgelist(nice_5_graph)
edgelist_6 = get.edgelist(nice_6_graph)

#get vertex labels
#vlabels_1 = get.vertex.attribute(nice_1_graph)


#get vertex fill color
#vfill = get.vertex.attribute(nice_1_graph)

#get vertext border color
#vborders = get.vertex.attribute(nice_1_graph, "border")

#get vertex degree
degrees_1 = degree(nice_1_graph)
degrees_2 = degree(nice_2_graph)
degrees_3 = degree(nice_3_graph)
degrees_4 = degree(nice_4_graph)
degrees_5 = degree(nice_5_graph)
degrees_6 = degree(nice_6_graph)

#get edge value
evalues_1 = get.edge.attribute(nice_1_graph, "r")
evalues_2 = get.edge.attribute(nice_2_graph, "r")
evalues_3 = get.edge.attribute(nice_3_graph, "r")
evalues_4 = get.edge.attribute(nice_4_graph, "r")
evalues_5 = get.edge.attribute(nice_5_graph, "r")
evalues_6 = get.edge.attribute(nice_6_graph, "r")

#betweenness scores
between_1 = betweenness(nice_1_graph, directed=FALSE)
between_2 = betweenness(nice_2_graph, directed=FALSE)
between_3 = betweenness(nice_3_graph, directed = FALSE)
between_4 = betweenness(nice_4_graph, directed = FALSE)
between_5 = betweenness(nice_5_graph, directed=FALSE)
between_6 = betweenness(nice_6_graph, directed = FALSE)

#get density 
nice_1_dens <- edge_density(nice_1_graph, loops=FALSE) #0.908 -harvest 
nice_2_dens <- edge_density(nice_2_graph, loops=FALSE) #0.544 -early growth
nice_3_dens <- edge_density(nice_3_graph, loops=FALSE) #0.346 - harvest
nice_4_dens <- edge_density(nice_4_graph, loops=FALSE) #0.949 - early growth
nice_5_dens <- edge_density(nice_5_graph, loops=FALSE) #0.477 - harvest
nice_6_dens <- edge_density(nice_6_graph, loops=FALSE) #0.345 - early growth

