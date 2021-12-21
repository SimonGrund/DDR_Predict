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

  gene = "BRCA1"; cohort = c("Ovary", "Breast"); Cutoff = 25
  
  #Filter data set
  d = filter(data, GENE == "BRCA1")
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
#x = dplyr::select(x, Signature.3, Signature.8, del.mh, clustered_inv_10.100Kb, clustered_inv_1.10Kb)
#x = x[,colSums(x)>0]
  
t = as.data.frame(cor(x))
#corrplot::corrplot(t, order = "hclust")
#library(psych)  
#p =  rownames(t)
t$feature = rownames(t)
t = filter(t, del.mh<1)
del2 = dplyr::select(t, del.mh)
ggplot(t, aes(x = reorder(feature, -del.mh), y = del.mh))+
  geom_histogram(stat = "identity")+
  theme(
    axis.text.x = element_text(angle = 45, hjust =1, vjust = 1)
  )

ggplot(t, aes(x = "Del.mh", y = del.mh))+
  geom_boxplot()+
  geom_point()+
  geom_text(aes(label = feature))
  
theme(
    axis.text.x = element_text(angle = 45, hjust =1, vjust = 1)
  )



# # K-Means Clustering with 5 clusters
# fit <- kmeans(x, 5)
# 
# # Cluster Plot against 1st 2 principal components
# 
# # vary parameters for most readable graph
# library(cluster) 
# clusplot(x, fit$cluster, color=TRUE, shade=TRUE, 
#          labels=2, lines=0)
# 
# # Centroid Plot against 1st 2 discriminant functions
# library(fpc)
# plotcluster(x, fit$cluster)
#     
#     
