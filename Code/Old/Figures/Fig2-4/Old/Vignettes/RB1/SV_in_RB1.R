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

gene = "RB1"; cohort = c("Urinary tract", "Uterus"); Cutoff = 25

#Filter data set
d = filter(data, GENE %in% gene)
d2 = anti_join(all_samples, d)%>%
  mutate(LOF = 0)

d = bind_rows(d, d2)%>%
  filter( primaryTumorLocation %in% cohort)%>%
  dplyr::select(-cancertype)

d$y[d$LOF >= 6] = 1
d$y[d$LOF < 2] = 0
table(d$y)
d = dplyr::select(d, Sample_ID, primaryTumorLocation, y)

#Features
sv = grep("cluster", colnames(features), value = T)
feats = dplyr::select(features, Sample_ID, sv)
feats = data.frame(Sample_ID = feats$Sample_ID, SV = rowSums(feats[,2:ncol(feats)]))

d2 = left_join(d, feats)
d2$y[is.na(d2$y)] = "VUS"
ggplot(d2, aes(x = y, y = SV, fill = primaryTumorLocation))+
  geom_boxplot()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))

#Fold change
dt = d2%>%
  group_by(primaryTumorLocation, y)%>%
  mutate(med = median(SV),
         y = ifelse(y == 1, "LOF", y),
         y = ifelse(y == 0, "WT", y))%>%
  distinct(y, med)%>%
  filter(y != "VUS")%>%
  pivot_wider(names_from = y, values_from = med)%>%
  mutate(
    fold_change = LOF/WT
  )

#Per feature
feats = dplyr::select(features, Sample_ID, all_of(sv))
feats = dplyr::distinct(feats, Sample_ID, .keep_all = T)
d3 = left_join(d, feats)
d3 = pivot_longer(d3, sv)

ggplot(d3, aes(x = as.factor(y), y = value, fill = as.factor(y)))+
  geom_boxplot()+
  facet_wrap(~name, scales = "free")

#Fold change
d3$y[is.na(d3$y)] = "VUS"
dt = d3%>%
  group_by(name, y)%>%
  mutate(med = median(value),
         y = ifelse(y == 1, "LOF", y),
         y = ifelse(y == 0, "WT", y))%>%
  distinct(name, y, med)%>%
  filter(y != "VUS")%>%
  pivot_wider(names_from = y, values_from = med)%>%
  mutate(
    fold_change = LOF/WT
  )


###SBS3

#Per feature
feats = dplyr::select(features, Sample_ID, non.clustered_inv_100Kb.1Mb)
feats = dplyr::distinct(feats, Sample_ID, .keep_all = T)
d4 = left_join(d, feats)

ggplot(d4, aes(x = as.factor(y), y = non.clustered_inv_100Kb.1Mb, fill = primaryTumorLocation ))+
  geom_boxplot()+
  geom_jitter()+
  ggpubr::stat_compare_means(comparisons = list(c(1,2)))

t = filter(d4, primaryTumorLocation == "Uterus")
wilcox.test(x = t$y, y = t$non.clustered_inv_100Kb.1Mb)

#Fold change
dt = d2%>%
  group_by(y)%>%
  mutate(med = median(SV),
         y = ifelse(y == 1, "LOF", y),
         y = ifelse(y == 0, "WT", y))%>%
  distinct(y, med)%>%
  filter(y != "VUS")%>%
  pivot_wider(names_from = y, values_from = med)%>%
  mutate(
    fold_change = LOF/WT
  )

#Fold change
d4$y[is.na(d4$y)] = "VUS"
d4$y[d4$y == 1] = "LOF"
d4$y[d4$y == 0] = "WT"
dt = d4%>%
  group_by(y)%>%
  mutate(med = median(non.clustered_inv_100Kb.1Mb))%>%
  distinct(y, med)%>%
  filter(y != "VUS")%>%
  pivot_wider(names_from = y, values_from = med)%>%
  mutate(
    fold_change = LOF/WT
  )
