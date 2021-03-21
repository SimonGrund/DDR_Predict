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

#IDH1 patients
source("Figures/Figure_4_new/Vignettes/GF/Extract_Variants.R")
IDH1_donors = Variant_extract(gene = "IDH1", cohort = unique(all_samples$primaryTumorLocation))
ATRX_donors = Variant_extract(gene = "ATRX", cohort = unique(all_samples$primaryTumorLocation))
both= filter(ATRX_donors, Sample_ID %in% IDH1_donors$Sample_ID)

#Increase in SV overall
source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment_SUMMED_FEAT.R")
ATRX =feat_enrich_summedFeat(gene = "ATRX",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                  feat = grep("cluster", colnames(features), value = T)
)
ATRX2 = d2
IDH1 =feat_enrich_summedFeat(gene = "IDH1",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = grep("cluster", colnames(features), value = T)
)
IDH12 = d2
BRCA =feat_enrich_summedFeat(gene = "BRCA2",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = grep("cluster", colnames(features), value = T)
)
BRCA2 = d2
d = bind_rows(ATRX, IDH1, BRCA)%>%na.omit()%>%
  rename(Feature = name)%>%
  pivot_longer(cols = c(LOF, WT))

ggplot(d, aes(x = Gene, y = value, fill = name))+
  geom_histogram(width= 0.4, stat = "identity", position = "dodge")+
  facet_wrap(~primaryTumorLocation, scales = "free")

d2 = bind_rows(ATRX2, IDH12, BRCA2)%>%na.omit()%>%
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


d2$Feature = "Summed #SVs per patient"
#d2$gene = factor(d2$Gene, levels = c("IDH1+ATRX", "ATRX", "IDH1", "BRCA2"))
ggplot(filter(d2, CDK12 != "VUS"), aes(x = paste(Gene, "\np<", round(p,5)+0.000001, sep = ""), y = value+0.1, fill = CDK12))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 0.7, aes(group = CDK12),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.66))+
  facet_wrap(primaryTumorLocation~Feature, scales = "free", ncol = 2)+
  scale_y_log10()  +
  xlab("")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme( 
    panel.grid = element_blank()
  )+
  scale_fill_manual(values = c("goldenrod", "skyblue"))+
  scale_color_manual(values = c("grey", "brown"))+
  guides(fill = guide_legend(title = "Status"))



#BRCA Feature enrichment
source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment.R")

#Comparing with BRCA
ATRX =feat_enrich(gene = "ATRX",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                  feat = c("del.mh", "Signature.3")
)
ATRX2 = d2
BRCA =feat_enrich(gene = "BRCA2",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                  feat = c("del.mh","Signature.3"))

BRCA2 = d2

IDH1 =feat_enrich(gene = "IDH1",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                  feat = c("del.mh", "Signature.3"))
IDH12 = d2


d2 = bind_rows(ATRX2, IDH12, BRCA2)%>%na.omit()%>%
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

ggplot(filter(d2, CDK12 != "VUS"), aes(x = paste(Gene, "\np<", round(p,5)+0.0000001, sep = ""), y = value, fill = CDK12))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 0.7, aes(group = CDK12),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.66))+
  facet_wrap(primaryTumorLocation~Feature, scales = "free", ncol = 2)+
  scale_y_log10()  +
  xlab("")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme( 
    panel.grid = element_blank()
  )+
  scale_fill_manual(values = c("goldenrod", "skyblue"))+
  scale_color_manual(values = c("grey", "brown"))+
  guides(fill = guide_legend(title = "Status"))


#Increase in MMR signatures overall
source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment_SUMMED_FEAT.R")
ATRX =feat_enrich_summedFeat(gene = "ATRX",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
ATRX2 = d2
IDH1 =feat_enrich_summedFeat(gene = "IDH1",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
IDH12 = d2
BRCA =feat_enrich_summedFeat(gene = "BRCA2",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
BRCA2 = d2

d2 = bind_rows(ATRX2, IDH12, BRCA2)%>%na.omit()%>%
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

d2$Feature = "Summed #MMR2 + Sig 26 per patient"

ggplot(filter(d2, CDK12 != "VUS"), aes(x = paste(Gene, "\np<", round(p,5)+0.00001, sep = ""), y = value+0.1, fill = CDK12))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 0.7, aes(group = CDK12),
             position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.66))+
  facet_wrap(primaryTumorLocation~Feature, scales = "free", ncol = 2)+
  scale_y_log10()  +
  xlab("")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme( 
    panel.grid = element_blank()
  )+
  scale_fill_manual(values = c("goldenrod", "skyblue"))+
  scale_color_manual(values = c("grey", "brown"))+
  guides(fill = guide_legend(title = "Status"))



#Increase in MMR signatures summed ATRX and IDH1
source("Figures/Figure_4_new/Vignettes/GF/Feature_enrichment_SUMMED_FEAT.R")
ATRX =feat_enrich_summedFeat(gene = c("ATRX", "IDH1"),cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
ATRX2 = d2
IDH1 =feat_enrich_summedFeat(gene = c("ATRX", "IDH1"), cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
IDH12 = d2
BRCA =feat_enrich_summedFeat(gene = c("ATRX", "IDH1"), cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), 
                             feat = c("Signature.MMR2", "Signature.26")
)
BRCA2 = d2
d = bind_rows(ATRX, IDH1, BRCA)%>%na.omit()%>%
  rename(Feature = name)%>%
  pivot_longer(cols = c(LOF, WT))

ggplot(d, aes(x = Gene, y = value, fill = name))+
  geom_histogram(width= 0.4, stat = "identity", position = "dodge")+
  facet_wrap(~primaryTumorLocation, scales = "free")

d2 = bind_rows(ATRX2, IDH12, BRCA2)%>%na.omit()%>%
  rename(Feature = name)
dt = dplyr::select(d, -value)%>%filter(name == "LOF")%>%distinct()%>%dplyr::select(-name)
d2 = left_join(d2, dt, by = c("primaryTumorLocation", "Feature","Gene"))
d2$Feature = "Summed #MMR2 + Sig 26 per patient"

ggplot(filter(d2, CDK12 != "VUS"), aes(x = paste(Gene, "\np<", round(p,5), sep = ""), y = value+.01, fill = CDK12))+
  geom_boxplot(outlier.shape = NA)+
  geom_point(size = 0.1, alpha = 0.5, position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75))+
  facet_wrap(primaryTumorLocation~Feature, scales = "free")+
  xlab("")+
  scale_y_log10()+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme(
    panel.grid = element_blank()
  )+
  scale_fill_brewer()+
  guides(fill = guide_legend(title = "Status"))

#Investigate S.ind enrichment—not present.
msi = fread("//Volumes/GenomeDK/PCAWG/simong_SV/subprojects/MMR_MSI/Data/MSI_status_all_samples.tsv")

d = both_lof%>%
  distinct(Sample_ID)

d = left_join(d, msi)

d = filter(var, Gene == "BLM")%>%
  distinct(Sample_ID, .keep_all = T)
d = left_join(d, msi, by  = "Sample_ID")
