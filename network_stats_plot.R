library(tidyverse)
library(RColorBrewer)
timepoint <- 1:6
season <- c("late dry", "late dry", "mid rainy", "mid rainy",
            "late rainy", "late rainy")
density <- c(0.74, 0.2, 0.22, 0.22, 0.16, 0.08)
transitivity <- c(0.89, 0.62, 0.73, 0.68, 0.46, 0.38)
modularity <- c(0, 0.19, 0.18, 0.02, 0.24, 0)

data <- cbind(timepoint, season, density, transitivity, modularity) %>%
  as_tibble() %>%
  mutate(density = as.numeric(density),
         transitivity = as.numeric(transitivity),
         modularity = as.numeric(modularity))

ggplot(data, aes(x=timepoint, y=density, group=1)) +
  scale_y_continuous(limits=c(0,1), 
                     breaks=c(0, 0.2, 0.4, 0.6, 0.8, 1)) +
  geom_point(size=3, alpha=0.9, color="#66C2A5") +
  geom_line(size=2, alpha=0.8, color="#66C2A5") +
  theme_classic() +
  geom_point(inherit.aes=FALSE, aes(x=timepoint, y=transitivity), color="#FC8D62", size=3, alpha=0.8) +
  geom_line(inherit.aes=FALSE, aes(x=timepoint, y=transitivity, group=1), color="#FC8D62", size=2, alpha=0.7) +
  geom_point(inherit.aes=FALSE, aes(x=timepoint, y=modularity), color="#8DA0CB", size=3, alpha=0.8) +
  geom_line(inherit.aes=FALSE, aes(x=timepoint, y=modularity, group=1), color="#8DA0CB", size=2, alpha=0.7) +
  ylab("") +
  annotate("text", x=2.5, y=.35, label="Density", color="#66C2A5", size=6) +
  annotate("text", x=2.5, y=0.85, label="Transitivity", color="#FC8D62", size=6) +
  annotate("text", x=2.5, y=0.08, label= "Modularity", color="#8DA0CB", size=6) +
  xlab("")


