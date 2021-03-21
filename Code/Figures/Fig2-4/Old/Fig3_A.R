library(data.table)
library(tidyverse)
library(caret)
library(extrafont)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021//")

d = fread("Figures/Figure_3/Models_performance_allCohorts.tsv")

d = filter(d, Cohort != Test)

g1 = ggplot(d, aes(x = AUC))+
  geom_histogram(stat = "bin", fill = "grey", color = "black", bins = 20)+
  geom_vline(xintercept = 0.7, lty = 2, col = "grey")+
  coord_flip()+
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.92, 0.9),
    plot.margin = margin(0,0,0,0),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA)
  )+
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0,1))+
  ylab("#Models")+
  facet_grid(cols = vars(paste("All models\n", sep = "\n")), scales = "free", space = "free")
g1

d = filter(d, Cohort != Test)%>%mutate(ID = paste(Gene, Cohort, Test, sep =", "))%>%
  filter(AUC >= 0.7)

g2 = ggplot(d, aes(x = reorder(paste(Test, ", ", nfg, "/", nbg, sep =""), -AUC), y = AUC))+
  geom_histogram(stat = "identity", width = 0.2, fill = "grey", color = "black")+
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = c(0.92, 0.9),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA)
  )+
  facet_grid(cols = vars(paste(Gene, Cohort, sep = "\n")), scales = "free", space = "free")+
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0,1))+
  xlab("Test cohort")+
  ylab("")+
  geom_hline(yintercept = 0.7, lty = 2, col = "grey")

g2 + g1

library(patchwork)
layout = "
AAAAAB
"

out = g2+g1+
  plot_annotation(tag_levels = "a",
                  title = 'Figure 3',
                  subtitle = '',
                  caption = 'February 8th, 2021'
  ) +
  patchwork::plot_layout(design = layout) 
out

ggsave("Figures/Figure_3/Fig3.pdf", device = "pdf", width = 12, height = 6)

