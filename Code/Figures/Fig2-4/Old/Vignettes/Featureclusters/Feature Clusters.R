#Cluster features
d = fread("Data/COMBINED_featureCounts.tsv")

#Add ct


library(M3C)
tsne(d$data,labels=as.factor(colnames(d)))