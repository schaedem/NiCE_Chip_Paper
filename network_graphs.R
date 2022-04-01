library(ggraph)
library(igraph)
library(arcdiagram)
library(tidyverse)
library(tidygraph)
library(ggrepel)
library(gridExtra)
library(cowplot)
library(ggpubr)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")
source("correlation_dataframes.R")
#source("networks_and_statistics.R")
options(warnings=-1)

#Read in assay list: primer pairs, gene/org designation, pathway, process (function), unique hexcodes
assay_list <- read_csv("assay_list_2.csv") %>%
  select(c(gene_org, colorcode, pathway, process, primer_pair)) 
assay_list$value <- assay_list$primer_pair
assay_list <- assay_list %>%
  unique()

make_nodes_df <- function(graphobject, assay_list){
  
  graph <- graphobject
  assays <- assay_list
  
  graph_labels <- as_tibble(get.vertex.attribute(graph, "name")) %>%
    merge(assays)
  
  #matching vertices to N cycle processes
  graph <- set.vertex.attribute(graph, "process", V(graph), value=graph_labels$process)
  graph <- set.vertex.attribute(graph, "pathway", V(graph), value=graph_labels$pathway)
  graph <- set.vertex.attribute(graph, "gene_org", V(graph), value=graph_labels$gene_org)
  graph <- set.vertex.attribute(graph, "colorcode", V(graph), value=graph_labels$colorcode)
  
  process <- get.vertex.attribute(graph, "process")
  pathway <- get.vertex.attribute(graph, "pathway")
  gene_org <- get.vertex.attribute(graph, "gene_org")
  color <- get.vertex.attribute(graph, "colorcode")
  
  edges <- get.edgelist(graph) %>% as_tibble()
  
  degrees <- degree(graph) %>% as_tibble_col() %>% rename(degree = value)
  degreenames <- degree(graph) %>% as_tibble_row() 
  names <- colnames(degreenames) %>% unlist() %>% as_tibble_col() %>% rename(primer_pair=value)

  betweenness_dat <- betweenness(graph) %>%
    as_tibble_col() %>%
    rename(betweenness=value)
  
  betw_names <- betweenness(graph) %>%
    as_tibble_row()
  
  betw_names <- colnames(betw_names) %>%
    unlist() %>%
    as_tibble_col() %>%
    rename(primer_pair=value)
  
  closeness_dat <- closeness(graph) %>%
    as_tibble_col() %>% rename(closeness_centrality = value)
  
  close_names <- closeness(graph) %>% as_tibble_row() 
  close_names <- colnames(close_names) %>% as_tibble_col() %>% rename(primer_pair=value)
  
  degree_df <- bind_cols(degrees, names) 
  between_df <- bind_cols(betweenness_dat, betw_names)
  closeness_df <- bind_cols(closeness_dat, close_names)
  all <- merge(degree_df, between_df) %>% merge(closeness_df)

  # data frame with groups, degree, labels and id
  nodes <- data.frame(process, pathway, graph_labels, color, gene_org, id=1:vcount(graph)) %>%
    merge(all, by="primer_pair")
  
  return(nodes)
  
}


make_network_graph <- function(nodes, edgelist) {

  edgelist <- edgelist %>%
    as_tibble() %>%
    rename(value = V1)
  
  nodes <- nodes %>%
    mutate(pathway = ordered(pathway))
  
  #getting colors for the network graph
  colordf <- nodes %>%
    select(colorcode, gene_org, pathway) %>%
    arrange(pathway) %>%
    select(-pathway) %>%
    unique()
  pal <- colordf$colorcode
  
  nodes_g <- nodes %>% 
    mutate(gene_org = factor(gene_org, levels=colordf$gene_org),
           gene_org = ordered(gene_org))
  
  #normalize degree size; code from blog.schoschastics.net

  normalise <- function (x, from = range(x), to = c(0, 1)) {
    x <- (x - from[1])/(from[2] - from[1])
    if (!identical(to, c(0, 1))) {
      x <- x * (to[2] - to[1]) + to[1]
    }
    x
  }
  
  # map to the range you want
  nodes$degree <- normalise(nodes$degree, to = c(3, 11))
  
  # prepare data for edges
  graph.tidy <- tbl_graph(nodes = nodes_g, edges = as_tibble(edgelist), 
                          directed = FALSE, node_key = "value") 
  
  graph_out <- 
    ggraph(graph.tidy, layout= "linear", circular=TRUE, sort=process) + 
    geom_edge_arc(alpha=0.1, fold=TRUE) + 
    scale_edge_width(range = c(0.2, 2)) +
    geom_node_point(aes(size = I(degree), color=gene_org)) +
    scale_color_manual(values=pal) + 
    geom_node_text(aes(label = value), angle = 0, hjust = 0.5, 
                   nudge_y = -0.15, size = 2, fontface="bold") +   
    coord_cartesian(clip = "off") + 
    theme_graph() +
    theme(legend.position = "none") 
  
  return(graph_out)

}


#arrange networks in one graph
# 
# plot_grid(
#   # row 1
#   plot_grid(graph1, graph3, graph5, nrow = 1, labels = c('A', 'C', 'E')) +
#     theme(plot.background = element_rect(color = "black")),
#   
#   # row 2
#   plot_grid(graph2, graph4, graph6, nrow = 1, labels = c('B', 'D', 'F')) +
#     theme(plot.background = element_rect(color = "black")), 
#   
#   nrow = 2)
# 
# plot.list <- lapply(list(graph1,graph3, graph5, graph2,graph4,  graph6), 
#                     function(p) p + theme(plot.background = element_rect(color = "black")))
# 
# ggarrange(plotlist = plot.list)
