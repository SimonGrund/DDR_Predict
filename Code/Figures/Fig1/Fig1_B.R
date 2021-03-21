######
# Figure 1 — Summary figure 
# - 1. Overview of samples per study, cancertype and metastasis
######
library(data.table)
library(tidyverse)
library(extrafont)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021/Github//")

######
# Load data
######
d = fread("Data/Generated_dataframes/all_samples.tsv")
d$Study = "PCAWG"
######
# Make figure
####

d2 = d%>%distinct(Sample_ID, .keep_all = T)%>%
  group_by(primaryTumorLocation)%>%mutate(n_sort = n())%>%ungroup()%>%
  group_by(primaryTumorLocation, Study)%>%mutate(n = n())%>%ungroup()%>%
  distinct(primaryTumorLocation, Study, .keep_all = T)

gA <<- ggplot(na.omit(d2), aes(y = n, x = reorder(primaryTumorLocation, -n_sort), fill = Study))+
  geom_histogram(stat = "identity", color = "black", width = 0.6, size = 0.1)+
  xlab("")+
  ylab("#Donors")+
  scale_fill_manual(values = c("brown3", "darkblue"))+
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        legend.position = c(0.9, 0.8),
        axis.title = element_text(size = 7),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA))+
  guides(fill = guide_legend(title = ""))

gA

ggsave(plot = gA, filename = "Output/Figures/Figure1B.pdf", device = "pdf", width = 4.5, height = 3)