library(data.table)
library(tidyverse)
library(glmnet)

source("Code/Generating_models//Lasso_for_evaluating_performance.R") #CV metrics
source("Code/Generating_models//Lasso_for_final_model.R") #final model
source("Code/Generating_models//featureFormatting.R")




data = fread("Data/LOF_25_all.tsv")
all_samples = fread("Data/all_donors.tsv")
features = fread("Data/Signatures/Scaled_features_27_05_21.tsv")


colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
data = filter(data, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)
table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

lasso = function(gene = "BRCA2", cohort = "Breast",  trainStudy = "HMF", Biallelic = T, nPermutations = 50){
#gene = "BRCA2"; cohort = "Breast"; trainStudy = "HMF"; Biallelic = T

  #Filter data set
  d = filter(data, Gene == gene)
  d2 = anti_join(all_samples, d)%>%
    mutate(Mutation_level = 0)
  d = bind_rows(d, d2)%>%
    filter( primaryTumorLocation == cohort)%>%
    dplyr::select(-cancertype)
  if(Biallelic == T){
    d = filter(d, !Mutation_level %in% c(1,2)) #Remove VUS
    d = mutate(d, y = ifelse(Mutation_level == 3, 1,0))
  }else{
    d = filter(d, !Mutation_level %in% c(1)) #Remove VUS
    d = mutate(d, y = ifelse(Mutation_level >1, 1,0))
  }
  
  table(d$y)
  
  #Add features
  #features_noSv = dplyr::select(features, Donor_ID, Sample_ID, all_of(grep("Signature", colnames(features), value = T)))
  #features_noSv = dplyr::select(features, Donor_ID, Sample_ID, all_of(grep("cluster", colnames(features), value = T, invert = T)))
  
  d = left_join(d, features)
  #Remove all NA columns
  d = d %>% dplyr::select_if(~all(!is.na(.)))
  
  #Split train test data
  train_data = dplyr::filter(d, Study == trainStudy)
  table(train_data$y)
  
  #Shuffle train_data mutation status
  for(j in 1:nPermutations){
  train_data$y = sample(train_data$y, size = nrow(train_data), replace = F)
  
  #Assign weights to the training
  w = table(train_data$y)
  w1 = w[1]/nrow(train_data)
  w0 = w[2]/nrow(train_data)

  train_data$weights[train_data$y == 0] = w0
  train_data$weights[train_data$y == 1] = w1
  
#  hist(train_data$weights, breaks = 30)
  
  cols_to_remove = c( "Donor_ID", "Sample_ID", "y", "weights" ,"Study", "primaryTumorLocation", "Mutation_level")
  
  ####
  # CV.LASSO
  ####
  
  #Generate the model metrics from CV
  ModelMetrics = lassoLoopNestedCv(train_data = train_data, cols_to_remove = cols_to_remove, gene = gene, cohort = cohort, study = trainStudy)
  
    if(j == 1){
      out_metrics =ModelMetrics
    }else{
      out_metrics = bind_rows(out_metrics, ModelMetrics)
    }
  
  }

  output = list(ModelMetrics)
  return(out_metrics)
}


