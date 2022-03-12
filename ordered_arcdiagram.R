library(ggraph)
library(igraph)
library(arcdiagram)
library(tidyverse)
library(tidygraph)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")
source("correlation_dataframes.R")


assay_list <- read_csv("assay_list_2.csv") %>%
  select(c(gene_org, colorcode, pathway, process, primer_pair)) 
assay_list$value <- assay_list$primer_pair
assay_list <- assay_list %>%
  unique()

test <- nice_1_graph 
#test <- set.vertex.attribute(test, "group", V(test), value=LETTERS[1:10])
test_labels <- as_tibble(get.vertex.attribute(test, "name")) %>%
  merge(assay_list)

#matching vertices to N cycle processes
test <- set.vertex.attribute(test, "process", V(test), value=test_labels$process)
test <- set.vertex.attribute(test, "pathway", V(test), value=test_labels$pathway)
test <- set.vertex.attribute(test, "gene_org", V(test), value=test_labels$gene_org)
test <- set.vertex.attribute(test, "colorcode", V(test), value=test_labels$colorcode)

process <- get.vertex.attribute(test, "process")
pathway <- get.vertex.attribute(test, "pathway")
gene_org <- get.vertex.attribute(test, "gene_org")
color <- get.vertex.attribute(test, "colorcode")

edges <- get.edgelist(test) %>% as_tibble()

degrees <- degree(test)

# data frame with groups, degree, labels and id
nodes <- data.frame(process, degrees, pathway, test_labels, color, gene_org, id=1:vcount(test)) %>%
  mutate(process = factor(process, levels = c("ammonia_oxidation", "hydroxylamine_oxidation", "nitrite_oxidation", "nitrite_reduction", "nitric_oxide_reduction", "nitrous_oxide_reduction", "n_fixation")),
         pathway = factor(pathway, levels=c("nitrification", "denitrification", "comammox", "n_fixation")),
         pathway = ordered(pathway)) %>%
  as_tibble() 

#getting colors for the network graph
colordf <- nodes %>%
  select(colorcode, gene_org, pathway) %>%
  arrange(pathway) %>%
  select(-pathway) %>%
  unique()
pal <- colordf$colorcode

nodes <- nodes %>% 
  mutate(gene_org = factor(gene_org, levels=colordf$gene_org),
         gene_org = ordered(gene_org))

#ggraph arc diagram

# prepare data for edges
test.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = "value") 
#dev.off()
ggraph(test.tidy, layout = "linear", sort=gene_org) + 
  geom_edge_arc(alpha=0.3, fold=TRUE) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_point(aes(size = degrees, color=gene_org)) +
  scale_color_manual(values=pal) + 
  geom_node_text(aes(label = value), angle = 30, hjust = 1, nudge_y = -0.2, size = 2.5) +   
  coord_cartesian(clip = "off") + 
  theme_graph() +
  theme(legend.position = "top") 

#auto layout
test.tidy <- tbl_graph(nodes = nodes, edges = edges, directed = FALSE, node_key = "value") 
#dev.off()
ggraph(test.tidy) + 
  geom_edge_arc(alpha=0.1, fold=TRUE) + 
  scale_edge_width(range = c(0.2, 2)) +
  geom_node_point(aes(size = degrees, color=gene_org)) +
  scale_color_manual(values=pal) + 
  geom_node_text(aes(label = value), angle = 0, hjust = 0.1, nudge_y = -0.1, size = 3) +   
  coord_cartesian(clip = "off") + 
  theme_graph() +
  theme(legend.position = "top") 


#900 height; 650 width or as landscape PDF
