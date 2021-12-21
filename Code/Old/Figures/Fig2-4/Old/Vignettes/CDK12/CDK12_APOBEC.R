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

  gene = "CDK12"; cohort = c("Ovary", "Breast"); Cutoff = 25
  
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
  
 # d = filter(d, !Sample_ID %in% brca2$Sample_ID)
  
  
  d = dplyr::select(d, Sample_ID, primaryTumorLocation, CDK12)
  d = left_join(d, features)%>%
    distinct()
  
  d = dplyr::select(d, Sample_ID, primaryTumorLocation,CDK12, Signature.13)
  
  d2 = pivot_longer(d, cols = c(Signature.13))
  
 ggplot(d2, aes(x = name, y = value, fill = CDK12))+
   geom_boxplot()+
   #coord_cartesian(ylim = c(0,1e4))+
   facet_wrap(primaryTumorLocation~name, scales ="free")
 
#Get fold change 
d3 = d2%>%
  group_by(CDK12, primaryTumorLocation, name)%>%
  mutate(med = median(value))%>%
  distinct(CDK12, name, med)%>%
  filter(CDK12 != "VUS")%>%
  pivot_wider(names_from = CDK12, values_from = med)%>%
  mutate(Fold_change = LOF/WT)

#p-value
d3$"p"[d3$primaryTumorLocation == "Breast"] = wilcox.test(x = filter(d, CDK12 == "LOF", primaryTumorLocation == "Breast")$Signature.13, 
            y = filter(d, CDK12 == "WT", primaryTumorLocation == "Breast")$Signature.13)$p.value
 
d3$p[d3$primaryTumorLocation == "Ovary"] = wilcox.test(x = filter(d, CDK12 == "LOF", primaryTumorLocation == "Ovary")$Signature.13, 
                                                        y = filter(d, CDK12 == "WT", primaryTumorLocation == "Ovary")$Signature.13)$p.value



