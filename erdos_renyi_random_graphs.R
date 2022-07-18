library(tidyverse)
library(igraph)
library(SingleCaseES)

#This script:
#1.generates random networks with the same number of nodes and edges as the
#observed networks acccording to the Erdos-Renyi model.
#2. calculates network statistics for random graphs: 
#density, transitivity (clustering coefficient), modularity, average path length
#3. calculates the log response ratio between the observed and random network statistics

#graph 1:18 nodes, 127 edges
#graph 2:16, 38
#graph 3:14, 42
#graph 4:13, 38
#graph 5:15, 27
#graph 6:6, 8

#1. erdos.renyi.game(n,p.or.m,type = c("gnm"),directed = FALSE,loops = FALSE)

generate_random_graphs <- function(nodes, edgenum) {
  
  nodes<- nodes
  edgenum <- edgenum
  
  datalist <- list()
  
  for(i in 1:100) {
    random <- erdos.renyi.game(nodes, edgenum, type= "gnm", directed=FALSE)
    edges <- edge_density(random)
    trans <- transitivity(random)
    paths <- average.path.length(random)
    dat <- data.frame(edges, trans, paths, shortest)
    datalist[[i]] <- dat
  }
  
  all_data = do.call(rbind, datalist)
  return(all_data)
  
}

graph1 <- generate_random_graphs(18,127)
graph2 <- generate_random_graphs(16, 38)
graph3 <- generate_random_graphs(14, 42)
graph4 <- generate_random_graphs(13, 38)
graph5 <- generate_random_graphs(15, 27)
graph6 <- generate_random_graphs(6,8)

#2. Summarize random network stats
rsum_1 <- graph1 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths)) %>% mutate(i = 1)

rsum_2 <- graph2 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths))%>% mutate(i = 2)

rsum_3 <- graph3 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths))%>% mutate(i = 3)
rsum_4 <- graph4 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths))%>% mutate(i = 4)

rsum_5 <- graph5 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths))%>% mutate(i = 5)

rsum_6 <- graph6 %>%
  summarize(mean(edges), mean(trans), mean(paths),
            sd(edges), sd(trans), sd(paths))%>% mutate(i = 6)

all_rsum <- rbind(rsum_1, rsum_2, rsum_3, rsum_4, rsum_5, rsum_6)
lrr <- c(0.2, 0.73, 0.51, 0.37, 0.69, -0.11)
dens <- c(0.83, 0.28, 0.46, 0.42, 0.25, 0.38)

