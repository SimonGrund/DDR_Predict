library(data.table)
library(tidyverse)
library(caret)
library(glmnet)
library(patchwork)
library(extrafont)
library(cutpointr)
library(scales)
source("Code/Generating_models//Lasso_for_evaluating_performance.R") #CV metrics
source("Code/Generating_models//Lasso_for_final_model.R") #final model
source("Code/Generating_models//featureFormatting.R")

data = fread("Data/LOF_25_all.tsv") #Loss-of-function annotations across all patients and all genes
all_samples = fread("Data/all_donors.tsv") #All 6,065 patients in the analyzis
features = fread("Data/Signatures/Scaled_features_27_05_21.tsv") #The scaled feature counts

colnames(features) = str_replace_all(colnames(features), "-",".")
colnames(features) = str_replace_all(colnames(features), ">","_._")
data = filter(data, Sample_ID %in% features$Sample_ID)
all_samples = filter(all_samples, Sample_ID %in% features$Sample_ID)
table(!duplicated(data$Sample_ID))
features$Signature.1 = NULL

lasso = function(gene = "BRCA2", cohort = "Breast",  trainStudy = "HMF", Biallelic = F){
#gene = "BRCA2"; cohort = "Breast"; trainStudy = "HMF"; Biallelic = T #Test case
  
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
    limit = 3
  }else{
    d = filter(d, !Mutation_level %in% c(1)) #Remove VUS
    d = mutate(d, y = ifelse(Mutation_level >1, 1,0))
    limit = 2
  }

  d = left_join(d, features)
  #Remove all NA columns
  d = d %>% dplyr::select_if(~all(!is.na(.)))
  
  #Split train test data
  train_data = dplyr::filter(d, Study == trainStudy)
  table(train_data$y)
  
  test_data = dplyr::filter(d, Study != trainStudy)
  
  #Assign weights to the training
  w = table(train_data$y)
  w1 = w[1]/nrow(train_data)
  w0 = w[2]/nrow(train_data)

  train_data$weights[train_data$y == 0] = w0
  train_data$weights[train_data$y == 1] = w1
  
  #hist(train_data$weights, breaks = 30)
  
  cols_to_remove = c( "Donor_ID", "Sample_ID", "y", "weights" ,"Study", "primaryTumorLocation", "Mutation_level")
  
  ####
  # CV.LASSO
  ####
  
  #Generate the model metrics from CV
  ModelMetrics = lassoLoopNestedCv(train_data = train_data, cols_to_remove = cols_to_remove, gene = gene, cohort = cohort, study = trainStudy)
  
  #Generate final model
  Test = F
  if(sum(test_data$y == 1) >= 5){Test = T}
  Model = lassoLoop(train_data = train_data, nIterations = 1, nfold = 5, Null = F,
                    test_data = test_data, cols_to_remove = cols_to_remove, test = Test,
                    gene = gene, cohort = cohort, study = trainStudy)
  if(nrow(Model) > 1){ #Don't do anything if only intercept
    Model = formatFeatures(Model = Model, train_data = train_data, test_data = test_data, test = Test)
  }else{
    Model$Feature_increase = NA; Model$Feature_increase_test = NA; Model$Feature_increase_p = NA
  }
  
  #Format the model
  Model$name = str_replace_all(Model$name, "\\__._1", ">1")
  Model$name = str_replace_all(Model$name, "b\\.", "b-")
  Model$name = str_replace_all(Model$name, "1.1", "1-1")
  Model$name = str_replace_all(Model$name, "\\.", " ")
  Model$name = str_replace_all(Model$name, "del", "del.")
  Model$name = str_replace_all(Model$name, "rep", "rep.")
  Model$name = str_replace_all(Model$name, "mh", "Microhomology")
  Model$name = str_replace_all(Model$name, "\\_", " ")
  
  Features = dplyr::select(Model, Gene, Cohort, Study, name, coefficient, Feature_increase, Feature_increase_test, Feature_increase_p)
  Model = dplyr::select(Model, Gene, Cohort, Study, AUC, AUC_test,
                        F1, F1_test,
                        diff_means,
                        n_train, n_test)%>%
    distinct(Gene, Cohort, Study, .keep_all = T)
  # 
  # 
  #Add event types
  eventTypes_train = as.data.frame(table(
    filter(data, Gene == gene, primaryTumorLocation == cohort, Study == trainStudy, Mutation_level >= limit)$State
  ))%>%
    mutate(
      Gene = gene, Cohort = cohort, Study = trainStudy
    )%>%
    dplyr::select(Gene, Cohort, Study, Event = Var1, Freq)
  # 
  if(Test == T){
  eventTypes_test = as.data.frame(
    table(
      filter(data, Gene == gene, primaryTumorLocation == cohort, Study != trainStudy, Mutation_level >= limit)$State
  ))%>%
  mutate(
    Gene = gene,
    Cohort = cohort,
    Study = ifelse(trainStudy == "HMF", "PCAWG", "HMF")
  )%>%
  dplyr::select(Gene, Cohort, Study, Event = Var1, Freq)
  }else{
    eventTypes_test = data.frame(Gene = gene, Cohort = cohort, Study = ifelse(trainStudy == "HMF", "PCAWG", "HMF"))
  }
  # 
  events = bind_rows(eventTypes_train, eventTypes_test)%>%
    pivot_wider(names_from = c(Event, Study), values_from = Freq)
  Model = left_join(Model, events)
  print("Done Model")
  # 
  # #generate full model
  # 
  # 
   output = list(ModelMetrics, Model, Features)
  return(output)
}


