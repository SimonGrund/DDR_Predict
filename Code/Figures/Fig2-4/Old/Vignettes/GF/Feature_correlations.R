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


feat_correlate = function(gene = "BRCA2", cohort = c("Ovary", "Breast"), feat = "Signature.13"){
#  gene = "BRCA2"; cohort = c("Ovary", "Breast"); Cutoff = 25; feat = "del.mh"
  Cutoff = 25
  #Filter data set
  d = filter(data, GENE == gene)
  d2 = anti_join(all_samples, d)%>%
    mutate(LOF = 0, GENE = gene, Germline = F, Double = F, MAX_CADD= 0)
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation %in% cohort)%>%
    dplyr::select(-cancertype)
  d = d[d$MAX_CADD>Cutoff | d$MAX_CADD<5 |LOF == 7|LOF == 1 |LOF == 0, ]
  d = mutate(d, y = ifelse(MAX_CADD >= Cutoff|LOF == 7, 1,0))
  d$y[d$LOF < 2] = 0
  
  # #Remove BRCA2
  # brca2 = filter(data, GENE == "BRCA2")
  # brca2 = filter(brca2, LOF >= 6)
  # 
  # d = filter(d, !Sample_ID %in% brca2$Sample_ID)
  
  table(d$y)
  
  #Add features
 # d = all_samples
  d = left_join(d, features)

  cols_to_remove = c("MAX_CADD","Donor_ID", "VAF", "Sample_ID", "LOF","Double", "GENE", "y", "Germline", "LOF", "Study", "primaryTumorLocation")

x = (dplyr::select(d, -(cols_to_remove)))  
x = na.omit(x)
x = dplyr::select(x, all_of(colnames(x)[colSums(x)>0]))


#Get correlation of each
t = as.data.frame(cor(x))
t$feature = rownames(t)
del2 = dplyr::select(t, feature, all_of(feat))


#p-value for correlations
for(col in colnames(x)){
  wx = unlist(dplyr::select(x, all_of(feat)))
  wy = unlist(dplyr::select(x, all_of(col)))
  p = cor.test(wx, wy)$p.value
  tmp = data.frame(feature = col, p = p)
  if(col == colnames(x)[1]){
    out = tmp
  }else{
    out = bind_rows(out, tmp)
  }
}

#Combine
del2 = dplyr::select(t, feature, all_of(feat))
del2 = left_join(del2, out)
del2$q = p.adjust(del2$p, method = "bonferroni")
#del2 = left_join(del2, dplyr::select(pval, feature, p_val = del.mh), by = "feature")
#del2$q = p.adjust(del2$p_val, method = "bonferroni")
# 
# ggplot(del2, aes(x = feature, y = -log(p_val)))+
#   geom_point(data = filter(del2, q< 0.05), aes(col = "q<0.05"))+
#   geom_point(data = filter(del2, q> 0.05), col = "grey", alpha = 0.6)+
#   geom_text(data = filter(del2, q< 0.05), aes(col = "red", label = paste(feature, sep ="")))

del2$correlation = del2[,feat]
g <<- ggplot(del2, aes(x = reorder(feature, -correlation), y = correlation))+
  geom_histogram(stat = "identity")+
  #geom_text(aes(label = paste("q=", round(q, 10), sep ="")))+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme(
    axis.text.x = element_text(angle = 45, hjust =1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(10,10,10,55)
  )+
  ggtitle(paste("Correlations with ", feat, sep = ""))+
  xlab("")+
  ylab("Pearson correlation")

del2
}

t =  feat_correlate(feat = "Signature.7",  cohort = c("Urinary tract", "Uterus"), gene = "RB1")

t = feat_correlate(feat = "del.mh",  cohort = c("Breast", "Ovary", "Prostate", "Pancreas"), gene = "BRCA1")
t = feat_correlate(gene = "BRCA1", cohort = c("Breast", "Ovary", "Prostate"), feat = "non.clustered_inv_1.10Kb")
t
g
#t
t = feat_correlate(gene = "CDK12", cohort = c("Breast", "Ovary", "Prostate", "Skin"), feat = c("Signature.19", "Signature.7","non.clustered_inv_100Kb.1Mb"))


t
g

t = feat_correlate(gene = "MTOR", cohort = c("Lung", "Skin"), feat = "non.clustered_inv_100Kb.1Mb")
t
g

t = feat_correlate(gene = "RB1", cohort = c("Urinary tract", "Uterus"), feat = "non.clustered_inv_100Kb.1Mb")
t
g