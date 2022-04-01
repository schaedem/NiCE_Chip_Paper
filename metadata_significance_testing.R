library(tidyverse)
library(lubridate)
library(SNCAnalysis)
library(lme4)
head(full_dat)

full_dat <- full_dat %>%
  mutate(timepoint = as.factor(timepoint),
         date = case_when(timepoint == 1 ~ mdy("9/16/2020"),
                          timepoint == 2 ~ mdy("10/6/2020"),
                          timepoint ==3 ~ mdy("11/18/2020"),
                          timepoint==4 ~ mdy("12/7/2020"),
                          timepoint==5 ~ mdy("01/19/2021"),
                          timepoint==6~ mdy("02/10/2021"))) %>%
  mutate(date = as.factor(date))
  

ggplot(data=full_dat, aes(log_dea, NP, color=as.factor(date))) +
  geom_point() +
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_color_viridis(discrete=TRUE) +
  labs(fill = "Timepoint", color="") 

dea_plot <- 
  ggplot(data=full_dat, aes(x=date, y=log_dea, fill=as.factor(date))) +
  geom_boxplot() +
  geom_jitter(alpha=0.3)+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_fill_viridis(discrete=TRUE, alpha=0.8) +
  labs(fill = "Timepoint", fill="") +
  ylab("log DEA") +
  xlab("")


np_plot <- 
  ggplot(data=full_dat, aes(x=date, y=NP, fill=as.factor(date))) +
  geom_boxplot() +
  geom_jitter(alpha=0.3)+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_fill_viridis(discrete=TRUE, alpha=0.8) +
  labs(fill = "Timepoint", fill="") +
  ylab("NP") +
  xlab("")

#Significance testing for time effects
np_mod <- lmer(NP~date+treatment+(1|date:location:block), data=full_dat)
np_mod2 <- lmer(NP~date+(1|date:location:block), data=full_dat)
summary(np_mod2)
anova(np_mod, np_mod2)

dea_mod <- lmer(log_dea~date+treatment+(1|date:location:block), data=full_dat)
dea_mod2 <- lmer(log_dea~date+(1|date:location:block), data=full_dat)
summary(dea_mod)
anova(dea_mod, dea_mod2)

dea_lm <- aov(log_dea ~ date + location:block, data=full_dat)
summary(dea_lm)
TukeyHSD(dea_lm, which="date")
#first two timepoints sig different from all others, including themselves; 
#last four timepoints are not sig different from themselves

np_lm <- aov(NP ~ date + location:block, data=full_dat)
summary(np_lm)
TukeyHSD(np_lm, which="date")
#first two timepoints sig different from all others, but not themselves; 
#last four timepoints are not sig different from themselves
