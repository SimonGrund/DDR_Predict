library(data.table)
library(tidyverse)
library(caret)
library(cutpointr)
library(scales)
library(patchwork)
library(extrafont)
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

#Result models — from looping the lasso all Data
Res = fread("Results/Full_valid_models.tsv")
#Res = filter(Res, significant == T)
  
models = list.files("Results/Valid_Models/", full.names = T)

df = data.frame(models)%>%
  separate(models, into = c(NA, "Gene", "Cohort", "Cutoff", NA), sep = "\\|", remove = F)%>%
  rowwise()%>%
  mutate(ID = paste(Gene, Cohort, sep ="|"))%>%
#  filter(ID %in% Res$ID)%>%
  mutate(Cutoff = as.numeric(Cutoff))

df = left_join(df, Res)

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
tmp2 = anti_join(all_samples, tmp, by = "Sample_ID")%>%
  mutate(LOF = 0, GENE = gene, Germline = F, Double = F, MAX_CADD = 0)
tmp = bind_rows(tmp, tmp2) #Should be 6072 rows now
tmp = tmp[tmp$MAX_CADD>Cutoff | tmp$MAX_CADD<5 |LOF == 7|LOF == 1 |LOF == 0, ]
tmp = mutate(tmp, y = ifelse(MAX_CADD >= Cutoff|LOF == 7, 1,0))
tmp$y[tmp$LOF < 2] = 0

# tmp = tmp[tmp$MAX_CADD>25 | tmp$MAX_CADD<5 |LOF == 7, ] #Filter away VUS
# tmp = mutate(tmp, y = ifelse(MAX_CADD >= 25|LOF == 7, 1,0))
# tmp$y[tmp$y == 1 & tmp$LOF == 1] = 0
# table(tmp$y)
#tmp = mutate(tmp, y = ifelse(LOF >= 6, 1,0))

tt =tmp%>%group_by(primaryTumorLocation)%>%
  mutate(n = sum(y == 1),
         n_bg = sum(y == 0)
         )%>%
  distinct(primaryTumorLocation, n, n_bg)%>%filter(n>=5 & n_bg>=5)

if(nrow(tt) == 0){next} 

#Make feature plot
t = as.data.frame(summary(c)$coefficients)
Features = data.frame(Gene = gene, Cohort = cohort, 
                      Feature = rownames(t), 
                      Estimate = t$Estimate,
                      stderror = t$`Std. Error`,
                      t_val = t$`t value`, 
                      p_val = t$`Pr(>|t|)`)
Features$Feature = str_replace_all(Features$Feature, "\\__._1", ">1")
Features$Feature = str_replace_all(Features$Feature, "b\\.", "b-")
Features$Feature = str_replace_all(Features$Feature, "1.1", "1-1")
Features$Feature = str_replace_all(Features$Feature, "\\.", " ")
Features$Feature = str_replace_all(Features$Feature, "del", "del.")
Features$Feature = str_replace_all(Features$Feature, "rep", "rep.")
Features$Feature = str_replace_all(Features$Feature, "mh", "Microhomology")
Features$Feature = str_replace_all(Features$Feature, "\\_", " ")

tmpo = filter(tmp, primaryTumorLocation == df$Cohort[row])
b1 = table(tmpo$y)

p1 = ggplot(Features, aes(x= Estimate, y = Feature))+
  geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
  theme_bw(base_size = 9, base_family = "Arial")+
  theme(
    # legend.position = c(0.75, 0.1),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = NA)
  )+
  geom_vline(xintercept = 0, lty = 2, col = "grey")+
  ylab("")
#p1

#Loop through ct
first = T
for(row2 in 1:nrow(tt)){
  
  tryCatch({ #Need a try catch in case of missing signatures/all == zero
  cat("\nTest: " ,tt$primaryTumorLocation[row2])
  tmp_CT = filter(tmp, primaryTumorLocation == tt$primaryTumorLocation[row2])
  numerics = filter(features, Sample_ID %in% tmp_CT$Sample_ID)
  numerics2 = dplyr::select(numerics, -Sample_ID, -Donor_ID)
  numerics2 = scale(log(numerics2+0.1))
  
  for(col in colnames(numerics2)){ 
    tmp_col = numerics2[,col]
    tmp_col = rescale(tmp_col, to = c(0,1)) #Rescale to get values pushed to the positive interval 0—1
    numerics2[,col] = tmp_col
  }
  
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
  
  AUC_tmp = as.data.frame(cut$roc_curve)
  AUC_tmp$label = paste(tt$primaryTumorLocation[[row2]], ", AUC=", round(cut$AUC, 2), sep ="")
  AUC_tmp$AUC = cut$AUC
  out_tmp = data.frame(Gene = gene, Cohort = cohort, Test = tt$primaryTumorLocation[row2],
                       AUC = cut$AUC, nfg = sum(p$obs == "LOF"), nbg = sum(p$obs != "LOF"),
                       AUC_discovery = df$ROC[row], n_discovery = df$Train_fg[row])
  if(first == T){
    AUC = AUC_tmp
    inner_cycle <<- out_tmp
    first = F
  }else{
    AUC = bind_rows(AUC, AUC_tmp)
    inner_cycle = bind_rows(inner_cycle, out_tmp)
  }
  
  }, error=function(cond){print("Missing signature!")})
}

if(first == T){next} #If no result from model

colorBlindBlack8  <- c("#000000", "#E69F00", "#56B4E9", "#009E73", 
                       "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

#Make AUC plot
p2 = ggplot()+
  geom_line(data = filter(AUC, AUC>= 0.6492155), aes(x = 1-tnr, y = tpr, color = label))+ #Update cutoff
  geom_line(data = filter(AUC, AUC < 0.6492155), aes(x = 1-tnr, y = tpr, group = label), color = "grey")+ #Update cutoff
  theme_bw(base_size = 9, base_family = "Arial")+
  theme(
    # legend.position = c(0.8, 0.2),
    legend.position = "right",
    legend.direction = "vertical",
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.border = element_rect(color = "black", fill = NA),
    axis.line = element_line(color = NA)
  )+
  xlab("1-Specificity")+
  ylab("Sensitivity")+
  scale_color_manual(values = colorBlindBlack8)+
 # scale_color_manual(values = c("Darkblue", "brown"))+
  guides(color = guide_legend(title = "", nrow = 4, byrow = T))
p2


p3=p1+p2 + 
  plot_annotation(
    subtitle = paste( gene, ", ", cohort, " cancer: WT: ", b1[1], ", LOF: ", b1[2], sep ="")
  ) & plot_layout(nrow = 1) & theme(text = element_text("Arial"))
p3

ggsave(plot = p3, paste("Figures/Figure_4_new/SubPlots/",gene, cohort,".pdf", sep = ""), device = "pdf", width = 7.2, height = 1.3)

if(outer_first == T){
  AUC_result = AUC
  result = inner_cycle
  outer_first = F
}else{
  AUC_result = bind_rows(AUC_result, AUC)
  result= bind_rows(result, inner_cycle)
}

}

#result = filter(result, Test != Cohort)
result = distinct(result)
n_distinct(result[c("Gene", "Cohort")])

write.table(result, "Results/Models_performance_allCohorts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)
write.table(result, "Figures/Figure_4_new//Models_performance_allCohorts.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

