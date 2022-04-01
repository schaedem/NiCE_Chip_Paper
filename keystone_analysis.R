library(tidyverse)
library(lme4)
library(stats)
library(pscl)
source("networks_and_statistics.R")

#abandon ship on modeling endeavor using betweenness; Works with closeness centrality!
#little to no identifiable relationship between degree and node; predictions from
#zero-inflated model do not make sense

all_nodes$between <- all_nodes$betweenness
all_nodes <- all_nodes %>%
  select(-betweenness) %>%
  mutate(log_between = log10(between +1),
         log_degree = log10(degree),
         count_between = round(between))

mod1 <- lmer(closeness_centrality ~ log_degree + between+ (1|location:timepoint), data=all_nodes)
mod2 <- lmer(closeness_centrality ~ log_degree +(1|location:timepoint), data=all_nodes)
summary(mod1)
anova(mod1, mod2) #mod2 is preferred

#mod2 <- glm(between ~ log_degree + location:timepoint, family="quasipoisson", data=all_nodes )
#mod3 <- glm(between ~ degree + location:timepoint, family= "gaussian", data=all_nodes)
#mod4 <- zeroinfl(count_between ~ degree |location:timepoint, data=all_nodes)
#mod5 <- zeroinfl(count_between ~ degree + location:timepoint, data=all_nodes)
#mod6 <- zeroinfl(count_between ~ log_degree +location:timepoint, data=all_nodes)

#18 is maximum node degree in network
#log10=1.255273
dat <- data.frame(log_degree=1.255273, timepoint=6, location="Rubona")
preds <- predict(mod2,dat)
max(preds) #max is for 11 degrees


