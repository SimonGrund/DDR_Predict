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

  gene = "ARID1A"; cohort = c("Billiary", "Breast", "Kidney", "Stomach", "Bone/Soft tissue", 
                              "Esophagus", "NET", "Unknown"); Cutoff = 25
  
  #Filter data set
  d = filter(data, GENE == gene)
  d2 = anti_join(all_samples, d)%>%
    mutate(LOF = 0, GENE = gene, Germline = F, Double = F, MAX_CADD= 0)
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation %in% cohort)%>%
    dplyr::select(-cancertype)
  d$CDK12[d$LOF>= 6]= "LOF" 
  d$CDK12[d$LOF <3] = "WT"
  d$CDK12[is.na(d$CDK12)] = "VUS"
  
  #Remove BRCA1/2 mutated
  brca2 = filter(data, GENE %in% c("BRCA2", "BRCA1"))
  brca2 = filter(brca2, LOF >= 6)
  
  #d = filter(d, !Sample_ID %in% brca2$Sample_ID)
  
  
  d = dplyr::select(d, Sample_ID, CDK12)
  d = left_join(d, features)
  
  d = dplyr::select(d, Sample_ID, CDK12, Signature.2, Signature.13, Signature.3, Signature.8, del.mh, non.clustered_inv_1.10Kb, 
                    non.clustered_inv_10.100Kb, non.clustered_inv_1Mb.10Mb)
  
  d2 = pivot_longer(d, cols = c(Signature.2, Signature.13, Signature.3, Signature.8, del.mh, non.clustered_inv_1.10Kb, 
                                non.clustered_inv_10.100Kb, non.clustered_inv_1Mb.10Mb))
  
 ggplot(d2, aes(x = name, y = value, fill = CDK12))+
   geom_boxplot()+
   #coord_cartesian(ylim = c(0,1e4))+
   facet_wrap(~name, scales ="free")
 
d3 = d2%>%
  group_by(CDK12, name)%>%
  mutate(med = median(value))%>%
  distinct(CDK12, name, med)%>%
  filter(CDK12 != "VUS")%>%
  pivot_wider(names_from = CDK12, values_from = med)%>%
  mutate(Fold_change = LOF/WT)

t = cor(d[,3:ncol(d)])

 
 