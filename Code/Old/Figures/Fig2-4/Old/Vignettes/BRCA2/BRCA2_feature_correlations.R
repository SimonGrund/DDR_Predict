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

  gene = "BRCA2"; cohort = c("Ovary", "Breast"); Cutoff = 25
  
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
#x = dplyr::select(x, Signature.3, Signature.8, del.mh, clustered_inv_10.100Kb, clustered_inv_1.10Kb)
#x = x[,colSums(x)>0]

#Get correlation of each
t = as.data.frame(cor(x))

#Get a p-val for each — Doesnt work!
pval = matrix(NA, ncol(x), ncol(x),
              dimnames = list(colnames(x),
                              colnames(x)))
xm = as.matrix(x)
for (row in 1:ncols(x)) {
  for (col in  colnames(x)) {
    pval[i, j] <- cor.test(xm[i, ], xm[j, ])$p.value
  }
}
pval = as.data.frame(pval)
pval$"feature" =row.names(pval)

#corrplot::corrplot(t, order = "hclust")
#library(psych)  
#p =  rownames(t)
t$feature = rownames(t)
t = filter(t, del.mh<1)
del2 = dplyr::select(t, feature, del.mh)
del2 = left_join(del2, dplyr::select(pval, feature, p_val = del.mh), by = "feature")
del2$q = p.adjust(del2$p_val, method = "bonferroni")
# 
# ggplot(del2, aes(x = feature, y = -log(p_val)))+
#   geom_point(data = filter(del2, q< 0.05), aes(col = "q<0.05"))+
#   geom_point(data = filter(del2, q> 0.05), col = "grey", alpha = 0.6)+
#   geom_text(data = filter(del2, q< 0.05), aes(col = "red", label = paste(feature, sep ="")))

ggplot(del2, aes(x = reorder(feature, -del.mh), y = del.mh))+
  geom_histogram(stat = "identity")+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme(
    axis.text.x = element_text(angle = 45, hjust =1, vjust = 1),
    panel.grid = element_blank(),
    plot.margin = margin(10,10,10,55)
  )+
  ggtitle("Correlations with del.mh")+
  xlab("")+
  ylab("Pearson correlation")


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
