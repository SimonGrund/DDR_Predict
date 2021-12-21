library(data.table)
library(tidyverse)
library(caret)
library(cutpointr)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021//")

####
# Test models on other cohorts, genesets ..
###
d_lof = fread("Data/LOF_all.tsv")
#ct = distinct(d_lof, primaryTumorLocation)$primaryTumorLocation
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/COMBINED_featureCounts.tsv")
colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
d_lof = filter(d_lof, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)

# numerics = features%>%
#   dplyr::select(-Sample_ID)
# numerics <- sapply( numerics, as.numeric )
# numerics[is.na(numerics)] = 0
# 
# numerics = scale(log(numerics+0.1)) #Scale it allll together :D

#Result models
Res = fread("Figures/Figure_2/Models_booty.tsv")
Res = filter(Res, significant == T)
  
models = list.files("Results/Models//", full.names = T)

df = data.frame(models)%>%
  separate(models, into = c(NA, "Gene", "Cohort", "Cutoff", NA), sep = "\\|", remove = F)%>%
  rowwise()%>%
  mutate(ID = paste(Gene, Cohort, sep ="|"))%>%
#  filter(ID %in% Res$ID)%>%
  mutate(Cutoff = as.numeric(Cutoff))

df = left_join(df, dplyr::select(Res, ID, ROC, ROC_test))

####
# Loop across all models all cohorts
####
outer_first = T
for(row in 1:nrow(df)){
gene = df$Gene[row]
cohort = df$Cohort[row]
cat(gene, cohort)
c = readRDS(df$models[row])  

#Filter by the gene
tmp = filter(d_lof, GENE == gene)
tmp2 = anti_join(all_samples, tmp)%>%
  mutate(LOF = 0, GENE = gene, Germline = F, Double = F)
tmp = bind_rows(tmp, tmp2) #Should be 5991 rows now
tmp = tmp[tmp$MAX_CADD>25 | tmp$MAX_CADD<5 |LOF == 7, ] #Filter away VUS
tmp = mutate(tmp, y = ifelse(MAX_CADD >= 25|LOF == 7, 1,0))
tmp$y[tmp$y == 1 & tmp$LOF == 1] = 0
table(tmp$y)
#tmp = mutate(tmp, y = ifelse(LOF >= 6, 1,0))

tt =tmp%>%group_by(primaryTumorLocation)%>%
  mutate(n = sum(y == 1),
         n_bg = sum(y == 0)
         )%>%
  distinct(primaryTumorLocation, n, n_bg)%>%filter(n>=5 & n_bg>=5)

if(nrow(tt) == 0){next} 

#Loop through ct
first = T
for(row2 in 1:nrow(tt)){
  
  tryCatch({ #Need a try catch in case of missing signatures/all == zero
  cat("\nTest: " ,tt$primaryTumorLocation[row2])
  tmp_CT = filter(tmp, primaryTumorLocation == tt$primaryTumorLocation[row2])
  numerics = filter(features, Sample_ID %in% tmp_CT$Sample_ID)
  numerics2 = dplyr::select(numerics, -Sample_ID)
  numerics2 = scale(log(numerics2+0.1))
  numerics2[is.na(numerics2)] = 0
  numerics = data.frame(Sample_ID = numerics$Sample_ID, numerics2)
  tmp_CT = left_join(tmp_CT, numerics, by = "Sample_ID")
#  tmp_CT$"Signature.33" = 0
  ####Predict score on samples
  p = caret::predict.train(c, newdata = tmp_CT, type = "prob",
                           na.action = na.pass)

  p = p%>%
    mutate(
      Sample_ID = tmp_CT$Sample_ID,
      obs = tmp_CT$y,
      obs = ifelse(obs == 1, "LOF", "WT")
    )
  p$obs = factor(p$obs, levels = c("WT", "LOF"))

  cut = cutpointr::cutpointr(x = p$LOF, pos_class = "LOF", direction = ">=", class = p$obs, 
                             method = maximize_metric, metric = sum_sens_spec) #Could be changed for adaptive cutoff...
                            #method = oc_manual, cutpoint = df$Cutoff[row])
  
  out_tmp = data.frame(Gene = gene, Cohort = cohort, Test = tt$primaryTumorLocation[row2],
                       AUC = cut$AUC, nfg = sum(p$obs == "LOF"), nbg = sum(p$obs != "LOF"),
                       AUC_discovery = df$ROC[row], AUC_validation = df$ROC_test[row])
  if(first == T){
    inner_cycle = out_tmp
    first = F
  }else{
    inner_cycle = bind_rows(inner_cycle, out_tmp)
  }
  
  }, error=function(cond){print("Missing signature!")})
}

#Output
if(outer_first == T){
  result = inner_cycle
  outer_first = F
}else{
  result= bind_rows(result, inner_cycle)
}

}

result = filter(result, Test != Cohort)
result = distinct(result)

write.table(result, "Results/Models_performance_allCohorts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(result, "Figures/Figure_3//Models_performance_allCohorts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

