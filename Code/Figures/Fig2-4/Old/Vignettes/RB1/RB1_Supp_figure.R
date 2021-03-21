library(data.table)
library(tidyverse)
library(scales)
library(extrafont)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/")

data = fread("Data/LOF_all.tsv")
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
data = filter(data, Sample_ID %in% features$Sample_ID)
table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

scientific_10 <- function(x) {
  parse(text=gsub("1e\\+*", "10^", scales::scientific_format()(x)))
}

#Increase in SV overall
source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment_SUMMED_FEAT.R")
RB1 =feat_enrich_summedFeat(gene = "RB1",cohort = c("Urinary tract", "Uterus"), 
                  feat = grep("inv", colnames(features), value = T)
)
RB1 =feat_enrich_summedFeat(gene = "RB1",cohort = c("Urinary tract", "Uterus"), 
                            feat = "non.clustered_inv_100Kb.1Mb"
)
RB1 = d2


d2 = RB1%>%na.omit()%>%
  rename(Feature = name)

for(gene in unique(d2$Gene)){
  for(ct in unique(d2$primaryTumorLocation)){
    try({
    x = filter(d2, Gene == gene & primaryTumorLocation == ct & CDK12 == "LOF")
    y = filter(d2, Gene == gene & primaryTumorLocation == ct & CDK12 == "WT")
    d2$p[d2$Gene == gene & d2$primaryTumorLocation == ct] = wilcox.test(x$value, y$value)$p.value
    })
  }
}


d2$Feature = "Summed #SVs per patient"
d2 = d2%>%group_by(CDK12, Gene)%>%mutate(n = n())%>%ungroup()
d2 = filter(d2, CDK12 != "VUS")
d2 = mutate(d2, CDK12 = ifelse(CDK12 == "LOF", "RB1-d", "RB1-wt"))
#d2$gene = factor(d2$Gene, levels = c("IDH1+ATRX", "ATRX", "IDH1", "BRCA2"))
SV = ggplot(filter(d2, CDK12 != "VUS"), aes(x = CDK12, y = value+0.1))+
  geom_boxplot(outlier.shape = NA, width= 0.4, size = 0.3)+
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.5)+
  ggpubr::stat_compare_means(comparisons = list(c("RB1-d","RB1-wt")), paired = F, 
                             method = "wilcox.test", size = 2.5, family = "Arial")+
  facet_wrap(primaryTumorLocation~Gene, scales = "free_x", ncol = 3)+
  #scale_y_log10(label=scientific_10, limits = c(0.03,1.5e4))+
  xlab("")+
  ylab("#SVs")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme( 
    panel.grid = element_blank(),
    axis.title = element_text(size = 7)
  )+
  scale_fill_manual(values = c("goldenrod", "skyblue"))+
  scale_color_manual(values = c("grey", "brown"))+
  guides(fill = guide_legend(title = "Status"))
  geom_text(data = distinct(d2, n, .keep_all = T), size = 2.5, family = "Arial",
            aes(x = CDK12, y = 0.05, label = paste("n=",n, sep= "")))

SV

#Increase in MMR signatures summed ATRX and IDH1
#source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment_SUMMED_FEAT.R")
ATRX =feat_enrich_summedFeat(gene = c("ATRX"),cohort = c("CNS"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
ATRX2 = d2
IDH1 =feat_enrich_summedFeat(gene = c("IDH1"), cohort = c("CNS"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
IDH12 = d2
d2 = bind_rows(ATRX2, IDH12)%>%na.omit()%>%
  rename(Feature = name)

#Annotate patients with mutation in both IDH1 and ATRX
ATRX = filter(d2, Gene == "ATRX", CDK12 == "LOF")
IDH1 = filter(d2, Gene == "IDH1", CDK12 == "LOF")
both_lof = filter(ATRX, Sample_ID %in% IDH1$Sample_ID)%>%mutate(Gene="ATRX+IDH1")

#WT
ATRX = filter(d2, Gene == "ATRX", CDK12 == "WT")
IDH1 = filter(d2, Gene == "IDH1", CDK12 == "WT")
both_wt = filter(ATRX, Sample_ID %in% IDH1$Sample_ID)%>%mutate(Gene="ATRX+IDH1")

#Remove ATRX and IDH1 mutated samples with mutations in both
single_lof = filter(d2, Sample_ID %in% both_lof$Sample_ID & Gene %in% c("ATRX", "IDH1"))
d2 = bind_rows(d2, both_lof, both_wt)
d2 = anti_join(d2, single_lof)

for(gene in unique(d2$Gene)){
  for(ct in unique(d2$primaryTumorLocation)){
    try({
      x = filter(d2, Gene == gene & primaryTumorLocation == ct & CDK12 == "LOF")
      y = filter(d2, Gene == gene & primaryTumorLocation == ct & CDK12 == "WT")
      d2$p[d2$Gene == gene & d2$primaryTumorLocation == ct] = wilcox.test(x$value, y$value)$p.value
    })
  }
}


d2 = filter(d2, CDK12 != "VUS")
d2 = d2%>%group_by(CDK12, Gene)%>%mutate(n = n())%>%ungroup()
d2 = mutate(d2, CDK12 = ifelse(CDK12 == "LOF", "RB1-d", "RB1-wt"))
MMR = ggplot(filter(d2, CDK12 != "VUS"), aes(x = CDK12, y = value+0.1))+
  geom_boxplot(outlier.shape = NA, width= 0.4, size = 0.3)+
  geom_jitter(width = 0.1, size = 0.3, alpha = 0.5)+
  ggpubr::stat_compare_means(comparisons = list(c("RB1-d","RB1-wt")), paired = F, 
                             method = "wilcox.test", size = 2.5, family = "Arial")+
  facet_wrap(~Gene, scales = "free_x", ncol = 3)+
  scale_y_log10(label=scientific_10, limits = c(0.03,1.5e4))+
  xlab("")+
  ylab("MMR signatures\n(Signature 26 + Signature MMR2)")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme( 
    panel.grid = element_blank(),
    axis.title = element_text(size = 7)
  )+
  scale_fill_manual(values = c("goldenrod", "skyblue"))+
  scale_color_manual(values = c("grey", "brown"))+
  guides(fill = guide_legend(title = "Status"))+
  geom_text(data = distinct(d2, n, .keep_all = T), size = 2.5, family = "Arial",
            aes(x = CDK12, y = 0.05, label = paste("n=",n, sep= "")))
MMR

library(patchwork)
layout = "
A
B
"
out = SV+MMR + patchwork::plot_layout(design = layout)
#ggsave(plot = out, filename = "Supp_figures/ATRX.pdf", device = "pdf", height = 5, width = 4)


###Three samples with co-mutation with POLQ
t = filter(var, Gene == "ATRX")
t2 = filter(var, Gene == "POLQ")
t3 = filter(t2, Sample_ID %in% t$Sample_ID)
ID = c("CPCT02050396T", "CPCT02070311T")

d2$sig = F
d2$sig[d2$Sample_ID %in% ID] = T
