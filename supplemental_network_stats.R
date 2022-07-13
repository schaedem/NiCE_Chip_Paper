library(tidyverse)
library(igraph)
source("intersection_networks.R")

#Density
#edge_density

edge_density(k_1_graph)
edge_density(r_1_graph)
edge_density(k_2_graph)
edge_density(r_2_graph)
edge_density(k_3_graph)
edge_density(r_3_graph)
edge_density(k_4_graph)
edge_density(r_4_graph)
edge_density(k_5_graph)
edge_density(r_5_graph)
edge_density(k_6_graph)
edge_density(r_6_graph)

#Clustering coefficient
#transitivity

transitivity(k_1_graph)
transitivity(r_1_graph)
transitivity(k_2_graph)
transitivity(r_2_graph)
transitivity(k_3_graph)
transitivity(r_3_graph)
transitivity(k_4_graph)
transitivity(r_4_graph)
transitivity(k_5_graph)
transitivity(r_5_graph)
transitivity(k_6_graph)
transitivity(r_6_graph)

#Path length
#mean_distance
mean_distance(k_1_graph)
mean_distance(r_1_graph)
mean_distance(k_2_graph)
mean_distance(r_2_graph)
mean_distance(k_3_graph)
mean_distance(r_3_graph)
mean_distance(k_4_graph)
mean_distance(r_4_graph)
mean_distance(k_5_graph)
mean_distance(r_5_graph)
mean_distance(k_6_graph)
mean_distance(r_6_graph)
