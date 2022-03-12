library(tidyverse)
library(viridis)
library(ggpubr)
library(grid)
library(cowplot)
library(lubridate)
library(ggrepel)
library(directlabels)
library(gridExtra)
library(grid)
library(lattice)

nice <- read_csv("final_nice_std_data.csv") %>%
  mutate(primer_pair = paste(primer_F, primer_R, sep=" / "),
         timepoint=as.factor(timepoint))
assay <- read_csv("assay_list_2.csv")
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/metadata")
meta <- read_csv("SNC_metadata_2020.csv") %>%
  select(sample, gwc)
setwd("/Volumes/Backup_1/Marie/Thesis/Rwanda/Lab work/03_NiCE_Chip/final_data")


#log_conc_truesoil = log10(gene copies/dry soil (g))

  sum_dat <- nice %>%
    select(location, primer_pair, log_conc_truesoil, timepoint, gwc) %>%
    group_by(primer_pair, timepoint, location) %>%
    summarise(mean(log_conc_truesoil), mean(gwc)) %>% 
    rename('mean_val' = 'mean(log_conc_truesoil)',
           'mean_gwc' = 'mean(gwc)') %>%
    merge(assay, by="primer_pair") %>%
    mutate(primer_pair = factor(primer_pair, 
          levels=c("IGK3 / DVV_correct","nifHF / nifHR", "Arch-amoAFA / Arch-amoAR","Arch-amoAFB / Arch-amoAR","Arch-amoA-for / Arch-amoA-rev",
                   "Arch-amoAF / Arch-amoAR", "amoa_F1 / amoA_2R", "Gamo172_F1 / Gamo172_F1_R2", "Gamo172_F2 / Gamo172_F2_R1",
                   "Gamo172_F1 / Gamo172_F1_R", "hzocl1F1 / hzocl1R2", "haoF4 / haoR2", "NxrB169F / NxrB638R", "NxrB1F / NxrB1R",
                   "comaA-244F / comaA-659R", "comaB-244F / comaB-659R","V66 / V67","V17m / napA4R","narG1960F / narG2650r",
                   "W9F / t38R", "nirK876 / nirK1040","nirKC1F / nirKC1R","FlaCu / R3Cu", "nirKC2F / nirKC2R", "nirKC4F / nirKC4R", 
                   "nirK_166F / nirK_665R", "nirKfF / nirKfR","nirSC2F / nirSC2R","nirSC1F / nirSC1R", "nirSCd3aF / nirSR3cd", 
                   "nirSC3F / nirSC3R","qnorB2F / qnorB7R", "qnorB2F / qnorB5R", "norB2 / norB6", "cnorB-2F / cnorB-6R", "nosZ1F / nosZ1R",
                   "nosZ-F-1181 / nosZR-1880","NosZ912F / NosZ853R"," nosZ-II-F / nosZ-II-R_")))

#xmin <- c(0, 1.5, 12.5, 13.5)
#xmax <- c(1.5, 12.5,13.5, Inf)
#color <- c('#E0F4DA', '#F9E0DE', '#EDD9BA', '#DFEDFA')    
#group <- c("n_fixation", "nitrification", "comammox", "denitrification")
  
#rects <- cbind(xmin, xmax, color, group) %>% as_tibble()

  plot <- ggplot(data=sum_dat, aes(x=primer_pair, y=mean_val, 
                                   group=interaction(timepoint, primer_pair),fill=as.factor(timepoint))) +
    geom_bar(stat="identity", position=position_dodge()) +
    scale_fill_viridis(discrete=TRUE, option="D") +
    theme_classic() +
    #facet_wrap(~location, nrow=2) +
    xlab("") +
    ylab(expression('log(copies'~""~g^~'-1'~'soil)')) +
    theme(
          legend.position= "right",
          legend.title=element_text(size=10),
          axis.text.x = element_text(angle=60, vjust=1, hjust=1, size=6.5)) +
    guides(fill=guide_legend(title="Timepoint"))
    

plot1 <-
  plot  + annotate("rect",xmin=-Inf, xmax=1.5, ymin=-Inf, ymax=Inf, alpha=0.7, fill='#E0F4DA') +
  annotate("rect", xmin=1.5, xmax=12.5, ymin=-Inf, ymax=10, alpha=0.2, fill='#C03830', color="darkgrey") +
  annotate("rect", xmin=12.5, xmax=13.5, ymin=-Inf, ymax=10, alpha=0.35, fill='#C18C5D')+
  annotate("rect", xmin=13.5, xmax=Inf, ymin=-Inf, ymax=Inf, alpha=0.7, fill='#DFEDFA',color="darkgrey") +
  geom_bar(data=sum_dat, aes(x=primer_pair, y=mean_val, 
                             group=interaction(timepoint, primer_pair),fill=as.factor(timepoint)), 
           stat="identity", position=position_dodge()) +
  scale_y_continuous(expand=c(0,0))
  


plot2 <- add_sub(plot1, c("N Fixation","Nitrification","Comammox", "Denitrification"), 
                 x=c(-0.05,0.25,0.55,0.76), hjust=c(0,0,0,0), vjust=c(-1.2,-1.2,-1.2,-1.2), size=12,
                 color = c("#44A043", '#C03830', '#C18C5D','#3A68AE'))

final_plot <- ggdraw(plot2)

ggsave("timepoint_primer_abund_both_locs.tiff", dpi=300, height=5, width=8.5)

#trendline graph
date_df <- sum_dat %>%
  mutate(date = case_when(timepoint==1 ~ mdy("09/17/2020"),
                          timepoint==2 ~ mdy("10/07/2020"),
                          timepoint==3 ~ mdy("11/19/2020"),
                          timepoint==4 ~ mdy("12/07/2020"),
                          timepoint==5 ~ mdy("01/20/2021"),
                          timepoint==6 ~ mdy("02/11/2021"))) %>%
  group_by(pathway, date) %>%
  summarise(mean(mean_val)) %>%
  rename(mean_val = "mean(mean_val)") %>%
  as_tibble()
 

date_df$pathway <- as.factor(date_df$pathway)
levels(date_df$pathway) <- c("Comammox", "Denitrification", "N fixation", "Nitrification")

#line plot
colors <- c('#C18C5D', '#3A68AE',"#44A043",'#C03830')

line_plot <- 
  ggplot(data=date_df, aes(x=date, y=mean_val, color=pathway)) +
  expand_limits(x=c(mdy("09/15/2020"), mdy("03/01/2021")))+
  geom_line(alpha=0.8, size=1) +
  theme_classic() +
  xlab("") +
  ylab(expression('log(copies'~""~g^~'-1'~'soil)')) +
  scale_color_manual(values=colors) +
  theme(legend.position = "none") +
 # geom_dl(aes(label=pathway), method="last.points")
  geom_dl(aes(label=pathway), method=list(cex=3, dl.trans(x=x+0.2, y=y+0.1), "last.qp"))

#egg::ggarrange(plot1, line_plot, nrow=2, heights=c(1,0.5), widths=c(1,0.5))

#cowplot::plot_grid(plot1, line_plot,  nrow=2, rel_heights = c(5/8, 3/8),
#                   rel_widths=c(1, 1/2))

#gridExtra::grid.arrange(arrangeGrob(plot1, line_plot, ncol=1, nrow=2),
                       #  heights=c(3,1), widths=c(2,1))
ggdraw(xlim=c(0,1), ylim=c(-1.5,1)) +
  draw_plot(plot1, x=0.01, y=-0.65, width=1, height=1.6) +
  draw_plot(line_plot, x=0.01, y=-1.5, width=0.7, height=0.9)

ggsave("timepoint_two_panels.tiff", dpi=300, height=10, width=8.5)


gwc_plot <- 
  ggplot(nice, aes(x=gwc, y=log_conc_truesoil, color=timepoint)) +
  #geom_point(alpha=0.8) +
  geom_jitter(width=0.009, alpha=0.8)+
  scale_color_viridis(discrete=TRUE, option="D")  +
  theme_classic() +
 facet_wrap(~pathway) +
  geom_line(aes(x=gwc, y=log_conc_truesoil), stat="smooth", method="lm", size=2, 
            color="darkred", alpha=0.6, inherit.aes=FALSE) +
   stat_cor(aes(x=gwc, y=log_conc_truesoil),method="spearman", label.y.npc = "bottom", inherit.aes=FALSE) +
  xlab("Gravimetric Water Content") +
  ylab(expression('log(copies'~""~g^~'-1'~'soil)')) +
  theme(legend.position="top", legend.title=element_blank())
gwc_plot


