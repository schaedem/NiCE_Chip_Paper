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
  labs(fill = "Timepoint", color="") +
  geom_smooth(se=FALSE)

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
summary(dea_mod2)
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

#16S vs timepoint and GWC
std_16S_data <- read_csv("processed_16S_qPCR.csv")
head(std_16S_data)

ggplot(data=std_16S_data, aes(x=as.factor(timepoint), y=log_conc_16S, fill=as.factor(timepoint))) +
  geom_boxplot() +
  geom_jitter(alpha=0.3)+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_fill_viridis(discrete=TRUE, alpha=0.8) +
  labs(fill = "Timepoint", fill="") +
  ylab("log(16S) copies") +
  xlab("")

#with 16S copies/g soil

full_16S_data <- read_csv("full_adj_df.csv") %>%
  select(timepoint, location, treatment, gwc, log_conc_16S_truesoil, sample) %>%
  unique()
head(full_16S_data)

ggplot(full_16S_data, aes(x=as.factor(timepoint), y=log_conc_16S_truesoil, fill=as.factor(timepoint)))+
geom_boxplot() +
  geom_jitter(alpha=0.3)+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_fill_viridis(discrete=TRUE, alpha=0.8) +
  labs(fill = "Timepoint", fill="") +
  ylab("log(16S/g soil)") +
  xlab("")

ggplot(full_16S_data, aes(x=gwc, y=log_conc_16S_truesoil, color=timepoint)) +
  geom_jitter(alpha=0.8, width=0.001)+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_color_viridis(alpha=0.8) +
  labs(fill = "Timepoint", fill="") +
  ylab("log(16S/g soil)") +
  xlab("Gravimetric Water Content") +
  geom_smooth(se=FALSE, method="lm")

#16S vs NP and DEA
full_16S_data_meta <- full_16S_data %>%
  merge(full_dat, by=c("sample", "location", "timepoint", "treatment")) %>%
  rename(gwc = gwc.x)
head(full_16S_data_meta)


ggplot(full_16S_data_meta, aes(x=NP, y=log_conc_16S_truesoil)) +
  geom_jitter(alpha=0.8, width=0.001, aes(color=as.factor(timepoint)))+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_color_viridis(alpha=0.8, discrete=TRUE) +
  guides(color=guide_legend("Timepoint")) +
  ylab("log(16S/g soil)") +
  xlab("Nitrification Potential") +
  geom_smooth(se=FALSE, method="lm")

ggplot(full_16S_data_meta, aes(x=log_dea, y=log_conc_16S_truesoil)) +
  geom_jitter(alpha=0.8, width=0.001, aes(color=as.factor(timepoint)))+
  theme_classic() +
  facet_wrap(~location) +
  viridis::scale_color_viridis(alpha=0.8, discrete=TRUE) +
  guides(color=guide_legend("Timepoint")) +
  ylab("log(16S/g soil)") +
  xlab("log(Denitrification Enzyme Activity)") +
  geom_smooth(se=FALSE, method="lm")
