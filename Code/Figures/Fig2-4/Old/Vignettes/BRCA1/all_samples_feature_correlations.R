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

all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")

all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)
#table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

  #Add features
  d = all_samples
  d = left_join(d, features)

  #cols_to_remove = c("MAX_CADD","Donor_ID", "VAF", "Sample_ID", "LOF","Double", "GENE", "y", "Germline", "LOF", "Study", "primaryTumorLocation")

x = (dplyr::select(d, -Donor_ID, -Sample_ID, -primaryTumorLocation, -Study))  
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
del = dplyr::select(t, del.mh)
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
