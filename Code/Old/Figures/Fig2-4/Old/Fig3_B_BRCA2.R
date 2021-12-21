library(data.table)
library(tidyverse)
library(caret)
library(patchwork)  
library(extrafont)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021//")

#BRCA
c = readRDS("Figures/Figure_3/BRCA_model.Rdata")

t = as.data.frame(summary(c)$coefficients)
Features = data.frame(Gene = "BRCA1/2", 
                      Feature = rownames(t), 
                      Estimate = t$Estimate,
                      stderror = t$`Std. Error`,
                      t_val = t$`t value`, 
                      p_val = t$`Pr(>|t|)`)
Features$Feature = str_replace_all(Features$Feature, "\\__._1", ">1")
Features$Feature = str_replace_all(Features$Feature, "b\\.", "b-")
Features$Feature = str_replace_all(Features$Feature, "1.1", "1-1")
Features$Feature = str_replace_all(Features$Feature, "\\.", " ")
Features$Feature = str_replace_all(Features$Feature, "del", "del.")
Features$Feature = str_replace_all(Features$Feature, "rep", "rep.")
Features$Feature = str_replace_all(Features$Feature, "mh", "Microhomology")
Features$Feature = str_replace_all(Features$Feature, "\\_", " ")

Feature_plot = ggplot(Features, aes(x= Estimate, y = Feature))+
  geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
  facet_wrap(~paste(Gene, sep = ", "), scales = "free", nrow = 6)+
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = c(0.9, 0.1),
        panel.border = element_rect(color = "black", fill = NA),
        axis.line = element_line(color = NA)) +
  geom_vline(xintercept = 0, lty = 2, col = "grey")

Feature_plot

#Make ROC
CDK12 = fread("Figures/Figure_3/BRCA_predicted_data.tsv")
cp = cutpointr(data = CDK12, x = LOF, class = obs, pos_class = "LOF",
               method = maximize_metric, metric = sum_sens_spec,  direction = ">=", subgroup = primaryTumorLocation)
cp

t1 = as.data.frame(cp$roc_curve[[1]])
t1$"ct" = "Breast"; t1$AUC = cp$AUC[[1]]
t2 = as.data.frame(cp$roc_curve[[2]])
t2$"ct" = "Prostate"; t2$AUC = cp$AUC[[2]]
t3 = as.data.frame(cp$roc_curve[[3]])
t3$"ct" = "Ovary"; t3$AUC = cp$AUC[[3]]
t4 = as.data.frame(cp$roc_curve[[4]])
t4$"ct" = "Pancreas"; t4$AUC = cp$AUC[[4]]

tt = bind_rows(bind_rows(t1, t4), bind_rows(t2,t3))
tt = mutate(tt, col = paste(ct, ", AUC=", round(AUC, 2), sep = ""))

ROC = ggplot(tt, aes(x = 1-tnr, y = tpr, color = col))+
  geom_line(show.legend = T)+
  theme_bw(base_size = 12, base_family = "Arial")+
  theme(
    legend.position = c(0.7, 0.3),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = NA)
  )+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  scale_color_manual(values = c("cornflowerblue", "brown", "goldenrod", "purple"))+
  guides(color = guide_legend(title = ""))
ROC


p = fread("Figures/Figure_3/BRCA_predicted_data.tsv")
table(duplicated(p$Sample_ID))

#Plot
p$obs = p$Mutation_level
p =  mutate(p, 
            infered_state = ifelse(obs == 7, "Clinvar Pathogenic", "Uncertain significance"),
            infered_state = ifelse(obs == 6, "Pathogenic VUS (CADD > 30)", infered_state),
            infered_state = ifelse(obs == 5, "VUS (CADD > 20)", infered_state),
            infered_state = ifelse(obs == 4, "VUS (CADD > 15)", infered_state),
            infered_state = ifelse(obs == 3, "VUS (CADD > 10)", infered_state),
            infered_state = ifelse(obs == 2, "Benign VUS (CADD < 10)", infered_state),
            infered_state = ifelse(obs == 1, "Clinvar Benign", infered_state),
            infered_state = ifelse(obs == 0, "No variants", infered_state)
)

p$infered_state = factor(p$infered_state, levels = c( "No variants",
                                                      "Clinvar Benign", "Benign VUS (CADD < 10)",
                                                      "VUS (CADD > 10)", "VUS (CADD > 15)",
                                                      "VUS (CADD > 20)", "Pathogenic VUS (CADD > 30)",
                                                      "Clinvar Pathogenic"))

p = p%>%group_by(infered_state)%>%mutate(n= n(), label = paste(infered_state, ", n=" ,n, sep =""))

#hist(p$LOF, breaks = 1e2)
g1 <<- ggplot(p, aes(x = reorder(label, obs), y = LOF))+
  geom_boxplot(outlier.shape = NA, width = 0.3, fill = "grey")+
  geom_jitter(width = 0.05, size = 0.2, aes(color = primaryTumorLocation), show.legend = F)+
  #geom_hline(yintercept = cp$optimal_cutpoint, lty = 2, col = "grey")+
  theme_bw(base_size = 12, base_family = "Arial") +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "bottom",
    legend.direction = "vertical",
   # legend.position = c(0.92, 0.9),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA),
    plot.margin = margin(0,0,0,20)
  )+
  ylab("Posterior probability of CDK12 mutation")+
  xlab("")+ 
  scale_x_discrete(drop=FALSE)+
  scale_color_manual(values = c("cornflowerblue", "brown", "goldenrod", "purple"))+
  guides(color = guide_legend(title = ""))
g1

g2 <<- ggplot(data = p, aes(x = LOF, fill = primaryTumorLocation))+
  geom_histogram(stat = "bin",, color = "black", bins = 20)+
  coord_flip()+
  theme_bw(base_size = 12, base_family = "Arial") +
  geom_vline(xintercept = cp$optimal_cutpoint[[1]], lty = 2, col = "grey")+
  annotate(geom = "text", y = 250, x = cp$optimal_cutpoint[[1]], label = "Breast cutoff", family = "Arial", size  =3.2, vjust = 0)+
  geom_vline(xintercept = cp$optimal_cutpoint[[2]], lty = 2, col = "grey")+
  annotate(geom = "text", y = 250, x = cp$optimal_cutpoint[[2]], label = "Prostate cutoff", family = "Arial", size  =3.2, vjust = 0)+
  geom_vline(xintercept = cp$optimal_cutpoint[[3]], lty = 2, col = "grey")+
  annotate(geom = "text", y = 250, x = cp$optimal_cutpoint[[3]], label = "Ovary cutoff", family = "Arial", size  =3.2,vjust = 0)+
  geom_vline(xintercept = cp$optimal_cutpoint[[4]], lty = 2, col = "grey")+
  annotate(geom = "text", y = 250, x = cp$optimal_cutpoint[[4]], label = "Pancreas cutoff", family = "Arial", size  =3.2,vjust = 0)+
  theme(
    axis.text.y = element_blank(),
    axis.title.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid.major = element_blank(),
    axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
    panel.grid.minor = element_blank(),
    legend.position = c(0.92, 0.9),
    plot.margin = margin(0,0,0,0),
    panel.border = element_rect(color = "black", fill = NA),
    legend.background = element_rect(fill =NA),
    axis.line = element_line(color = NA)
  )+
  ylab("#Patients")+
  scale_fill_manual(values = c("cornflowerblue", "brown", "goldenrod", "purple"))+
  guides(fill = guide_legend(title = ""))+
  scale_y_continuous(limits = c(0,400))
# annotate(geom = "text", x = 0.4496159, y = 200, label = "Determined cutoff", size = 2.8, family = "Arial")
g2

g1+g2

library(patchwork)
layout = "
AABBCCD
"

out = Feature_plot + ROC + g1 + g2 +
  plot_annotation(title = '',
                  subtitle = 'c',
                  
  ) +
  patchwork::plot_layout(design = layout, guides = "collect")
out

ggsave("Figures/Figure_3/BRCA.pdf", device = "pdf", width  = 15, height= 6)

