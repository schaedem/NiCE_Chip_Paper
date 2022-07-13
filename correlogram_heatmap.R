library(tidyverse)
library(lubridate)
library(SNCAnalysis)
library(viridis)

#Correlogram with NiCE abundances vs soil data
assay <- read_csv("assay_list_2.csv") %>%
  select(primer_pair, acronym)

only_16S <- read_csv("16S_data.csv")[,2:4] %>%
  mutate(acronym = ifelse(grepl("Arch_16S", acronym), "16S_Arch", acronym))

final_nice <- read_csv("final_nice_std_data.csv") %>%
  select(sample, log_conc_16S_truesoil, log_conc_truesoil, pH,
         gwc, timepoint, location, primer_F, primer_R, assay) %>%
  mutate(pH = ifelse(is.na(pH), 5.14, pH),
         primer_pair = paste(primer_F, primer_R, sep=" / ")) %>%
  select(-c(primer_F, primer_R)) %>%
  merge(assay, by= "primer_pair")

primers <- final_nice %>%
  select(sample, log_conc_truesoil, acronym) %>%
  rbind(only_16S)

meta_dat <- full_dat %>%
  select(location, timepoint, sample)

long_data <- primers %>% 
  rbind(only_16S) %>%
  merge(meta_dat, by="sample")

spread_primers <- spread(primers, key=acronym, value=log_conc_truesoil) 

# Heatmap
nice.heatmap <- ggplot(data = long_data, mapping = aes(x = acronym, y = timepoint, fill = log_conc_truesoil))  + 
  geom_tile() + 
  ylab("Timepoint") +
  xlab("") +
  facet_grid(location~.) +  
  theme_classic() +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) +
  scale_fill_viridis(option="B",name = expression(~"log(copies"~g^-1~"soil)"))
nice.heatmap
ggsave("nice_heatmap.tiff", dpi=300, width=12, height = 7)

#wide dataframe for correlogram
full_spread <- spread_primers %>%
  merge(full_dat, by="sample") %>%
  select(-c(block, timepoint, location, mmol_kg_hr, 
            log_auc_dea, auc_dea, treatment, auc_np, sample)) %>%
  mutate(pH = ifelse(is.na(pH), 5.14, pH))

#plot correlogram -- having difficulty with margins
library(corrgram)
library(GGally)
library(ggstatsplot)

#instead of traditional correlogram, we can compute all pairwise correlations and 
#only show the significant ones
#devtools::install_github("laresbernardo/lares")
library(lares)
np<- corr_var(full_spread, NP, method = "spearman", max_pvalue = 0.05, top=200) + ggtitle("NP")
dea <- corr_var(full_spread, log_dea, method="spearman", max_pvalue=0.05, top=200) + ggtitle("DEA")
gwc <- corr_var(full_spread, gwc, method="spearman", max_pvalue=0.05, top=200) + ggtitle("GWC")
pox <- corr_var(full_spread, mg_kg_POXC, method="spearman", max_pvalue=0.05, top=200) + ggtitle("POXC")
pmn <- corr_var(full_spread, NH4N_net, method="spearman", max_pvalue=0.05, top=200) + ggtitle("PMN")
nit <- corr_var(full_spread, mg_kg_NO3N, method="spearman", max_pvalue=0.05, top=200) + ggtitle("NO3-N")
amm <- corr_var(full_spread, mg_kg_NH4N, method="spearman", max_pvalue=0.05, top=200) + ggtitle("NH4-N")

#correlogram as a tile plot
cormat <- round(cor(full_spread, method="spearman", use="complete.obs"),2)
melted_cormat <- melt(cormat)

melted_tile <- melted_cormat %>%
  filter(Var1 == "mg_kg_POXC" | Var1 == "mg_kg_NH4N" | Var1 == 
           "mg_kg_NO3N" | Var1 == "NH4N_net" | Var1 == "gwc" | Var1 == "pH" | 
           Var1 == "NP" | Var1 == "log_dea") %>%
  subset(!(Var2 %in% c("mg_kg_POXC", "mg_kg_NH4N", "mg_kg_NO3N", "NH4N_net",
                       "gwc", "pH", "NP", "log_dea")))
  na.omit()

ggplot(data = melted_tile, aes(Var2, Var1, fill = value))+
  geom_tile(color = "white")+
  scale_fill_gradient2(low = "darkred", high = "steelblue", mid = "white", 
                       midpoint = 0, limit = c(-1,1), space = "Lab", 
                       name="Spearman\nCorrelation") +
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed() +
  xlab("") + ylab("")

citation(package = "vegan")

#correlogram just for primers, showing significance

spread_final <- spread_primers %>%
  select(-sample)

cor_plot <- 
  ggstatsplot::ggcorrmat(spread_final,
                       type="nonparametric", 
                       colors=c("darkred", "white", "steelblue")) 
ggsave("primers_correlogram.tiff", width=30, height=30, dpi=300)


ggstatsplot::ggcorrmat(full_spread,
                         type="nonparametric", 
                         colors=c("darkred", "white", "steelblue")) 
