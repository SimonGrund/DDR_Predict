library(data.table)
library(tidyverse)
library(caret)
library(extrafont)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021//")

d = fread("Figures/Figure_3/Models_performance_allCohorts.tsv")
mods = fread("Figures/Figure_2/Models_booty.tsv")
mods = filter(mods, significant == T)
#d = mutate(d, ID = paste(Gene, Cohort, sep ="|"))%>%
 # filter(ID %in% mods$ID)

# for(i in 1:1e4){
#   set.seed(i)
#   tmp = d[sample(1:nrow(d), size = nrow(d), replace = T), ]
#   tmp_out = data.frame(AUC = mean(tmp$AUC))
#   if(i == 1){
#     out = tmp_out
#   }else{
#     out = bind_rows(out, tmp_out)
#   }
# }

out = fread("Figures/Figure_3/null_distribution.tsv")
out$"AUC" = out$AUC_v
quantile(out$AUC)
hist(out$AUC)

qAUC = quantile(out$AUC, 0.95)

g1 = ggplot(out, aes(x = AUC))+
  geom_histogram(stat = "bin", fill = "grey", color = "black", bins = 150)+
  geom_vline(xintercept = qAUC, lty = 2, col = "grey")+
  annotate(geom = "text", size = 2, family = "Arial", 
           label = "Top 5%", x = qAUC+0.1, y = 7e3, angle = 0)+
  coord_flip()+
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    strip.text = element_text(family = "Arial", size =4 ),
    legend.position = c(0.92, 0.9),
    plot.margin = margin(0,0,0,0),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA)
  )+
  scale_x_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0,1))+
  ylab("#Bootstraps")+
  facet_grid(cols = vars(paste("Bootstrap\nAUC", sep = "\n")), scales = "free", space = "free")
g1


# d = filter(d, Cohort != Test)%>%mutate(ID = paste(Gene, Cohort, Test, sep =", "))%>%
#   filter(AUC >= 0.7)

mods = fread("Figures/Figure_2/Models_booty.tsv")
mods = filter(mods, significant == T)
d2 = mutate(d, ID = paste(Gene, Cohort, sep ="|"))%>%
  filter(
     ID %in% mods$ID,
         AUC >= qAUC
        )%>%
  rowwise()%>%
  mutate(
    p = sum(out$AUC>AUC)/1e5+0.00001,
    fill = ifelse(AUC > qAUC, "Above cutoff", "Below cutoff")
  )
d2$q = p.adjust(d2$p, method = "bonferroni", n = nrow(d))

g2 = ggplot(d2, aes(x = reorder(paste(Test, ", ", nfg, "/", nbg, sep =""), -AUC), y = AUC))+
  geom_histogram(stat = "identity", width = 0.2, fill = "grey30", color = "black")+
  geom_text(aes(label = paste("q=", round(q,3))), angle = 90, hjust = 0, vjust = 0.5, nudge_y = 0.02, size =2.5, family = "Arial")+
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
 #   strip.text = element_text(family = "Arial", size =4 ),
   # legend.position = c(0.92, 0.9),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA)
  )+
  facet_grid(cols = vars(paste(Gene, Cohort, sep = "\n")), scales = "free", space = "free")+
  scale_y_continuous(breaks = c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1), limits = c(0,1))+
  scale_fill_manual(values= c("darkblue", "grey"), na.value = "grey")+
  xlab("Test cohort")+
  ylab("AUC")+
  geom_hline(yintercept = qAUC, lty = 2, col = "grey")+
  guides(fill = guide_legend(title = ""))
g2

g2 + g1

library(patchwork)
layout = "
AAAAAAAAAAB
"

pout = g2+g1+
  plot_annotation(
                  title = 'Figure 3',
                  subtitle = 'a',
                  caption = ''
  ) +
  patchwork::plot_layout(design = layout, guides = "collect") &theme(legend.position = "bottom",
                                                                     text = element_text("Arial", size = 9))
pout 

ggsave(plot = pout, "Figures/Figure_3/Fig3_A.pdf", device = "pdf", width = 7.5, height = 4)

write.table(d2, "Figures/Figure_3/cross_ct_models.tsv", sep = "\t", col.names = T, row.names = F)
