library(data.table)
library(tidyverse)
library(caret)
library(glmnet)
library(patchwork)
library(scales)
library(extrafont)
#library(gridExtra)
library(cutpointr)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/")

data = fread("Data/LOF_all.tsv")
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
data = filter(data, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)
table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

feat_enrich = function(gene = "CDK12", cohort = c("Ovary", "Breast"), feat = "Signature.13"){
#gene = "CDK12"; cohort = c("Ovary", "Breast", "Prostate", "Skin"); feat = c("Signature.13")
  
  #Filter data set
  d = filter(data, GENE %in% gene)
  d2 = anti_join(all_samples, d)%>%
    mutate(LOF = 0)
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation %in% cohort)%>%
    dplyr::select(-cancertype)
  d$CDK12[d$LOF>= 6]= "LOF" 
  d$CDK12[d$LOF <3] = "WT"
  d$CDK12[is.na(d$CDK12)] = "VUS"
  
 #  #Remove BRCA1/2 mutated
 #  brca2 = filter(data, GENE %in% c("BRCA2", "BRCA1"))
 #  brca2 = filter(brca2, LOF >= 6)
 #  
 # # d = filter(d, !Sample_ID %in% brca2$Sample_ID)
  
  d = dplyr::select(d, Sample_ID, primaryTumorLocation, CDK12)
  d = left_join(d, features)%>%
    distinct()
  
  d = dplyr::select(d, Sample_ID, primaryTumorLocation, CDK12, all_of(feat))

  d2 = pivot_longer(d, cols = all_of(feat))
  
 g<<- ggplot(d2, aes(x = name, y = value, fill = CDK12))+
   geom_boxplot()+
   #coord_cartesian(ylim = c(0,1e4))+
   facet_wrap(primaryTumorLocation~name, scales ="free")+
   ggpubr::stat_compare_means(comparisons = list(c(1,2)))
 
d3 = d2%>%
  group_by(CDK12, name, primaryTumorLocation)%>%
  mutate(med = round(median(value),0),
         CI_up = round(quantile(value, 0.875),0),
         CI_down = round(quantile(value, 0.125),0))%>%
  distinct(CDK12, name, med, CI_up, CI_down)%>%
  filter(CDK12 != "VUS")%>%
  pivot_wider(names_from = CDK12, values_from = c(med, CI_up, CI_down))%>%
  mutate(Fold_change = med_LOF/med_WT)

for(ct in distinct(d, primaryTumorLocation)$primaryTumorLocation){
  for(ff in feat){
    x = dplyr::filter(d, CDK12 == "LOF", primaryTumorLocation == ct)%>%dplyr::select(all_of(ff))
    y = dplyr::filter(d, CDK12 == "WT", primaryTumorLocation == ct)%>%dplyr::select(all_of(ff))
    if(nrow(y)<3 | nrow(x)<3){next}
    d3$"p"[d3$primaryTumorLocation == ct & d3$name == ff] = 
      wilcox.test(x = unlist(x[,1]), 
                y = unlist(y[,1]))$p.value
    d3$"n_lOF"[d3$primaryTumorLocation == ct & d3$name == ff] = nrow(x)
    d3$"n_wt"[d3$primaryTumorLocation == ct & d3$name == ff] = nrow(y)
  }
}

d3$"Gene"  = gene
d2$"Gene" = gene
d2 <<- d2
d3
}

# feat_enrich()
 t =feat_enrich(cohort = c("Ovary", "Breast"), feat = "Signature.13")
# t =feat_enrich(cohort = c("Ovary", "Breast", "Prostate", "Skin"), feat = "non.clustered_inv_100Kb.1Mb")
# 
 t =feat_enrich(gene = "CDK12", cohort = c("Ovary", "Breast", "Prostate", "Skin"), feat = "non.clustered_inv_100Kb.1Mb")
# t =feat_enrich(gene = "CDK12", cohort = c("Ovary", "Breast", "Prostate", "Skin"), feat = "non.clustered_inv_10.100Kb")
# t =feat_enrich(gene = "CDK12", cohort = c("Ovary", "Breast", "Prostate", "Skin"), feat = c("non.clustered_inv_1.10Kb","non.clustered_inv_10.100Kb"))
# t =feat_enrich(gene = "CDK12", cohort = c("Ovary", "Breast", "Prostate", "Skin"), feat = "del.mh")
# g
# t
# t =feat_enrich(gene = "BRCA1",cohort = c("Ovary", "Breast"), feat = "Signature.3")
# 
# t =feat_enrich(gene = "BRCA1",cohort = c("Ovary", "Breast", "Prostate"), feat = "non.clustered_inv_10.100Kb")
# t =feat_enrich(gene = "BRCA1",cohort = c("Ovary", "Breast", "Prostate"), feat = "Signature.8")
# t =feat_enrich(gene = "BRCA1",cohort = c("Ovary", "Breast", "Prostate"), feat = "non.clustered_inv_1.10Kb")
# g
# t
# 
# #Arid1A
# t =feat_enrich(gene = "ARID1A",cohort = c("Biliary", "Breast", "Kidney", "Stomach", "Bone/Soft tissue", 
#                                           "Esophagus", "NET", "Unknown"), feat = "Signature.2")
 t =feat_enrich(gene = "ARID1A",cohort = c("Kidney", "Prostate"), feat = "non.clustered_inv_1.10Kb")
# g
# t
# 
# #ATRX
# t =feat_enrich(gene = "ATRX",cohort = c("Bone/Soft tissue", "CNS"), feat = "del.rep")
# g
# t
# 
# #ATRX
# t =feat_enrich(gene = "ATRX",cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"), feat = "del.mh")
# g
# t
# 
# 
# #RB1
 t =feat_enrich(gene = "RB1",cohort = c("Urinary tract", "Uterus"), feat =  c("non.clustered_inv_100Kb.1Mb", "Signature.7"))
# g
# t
# 
# #MTOR
# t =feat_enrich(gene = "MTOR",cohort = c("Lung", "Skin"), feat = "del.none")
# g
# t
