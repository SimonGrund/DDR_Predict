library(data.table)
library(tidyverse)
library(caret)
setwd("/Volumes/Macintosh HD/Users/au460892/Documents/GitHub/DDR_Predict/")

#Load genes/cohorts to build with 
m= fread("Output/Step1_models/Models.tsv")
m = filter(m, ROC > 0.72, ROC_test > 0.65) #These cutoffs we learned from botstrapping all generated models in the paper
# We now have one high-confidence model to retrain!

source("Code/Generating_models/Step2_retrain_highconfidence_models/Retrain_models_full_data.R")

model_LOOP = "Empty"
for( row_number in 1:nrow(m)){
  cat(row_number, "/", nrow(m))
  row = m[row_number,]
  cat("\n\nGene=", row$Gene, "\nCohort=", paste(row$Cohort, collapse = ", ") , "\n")
  tryCatch({
    e_message = "Not known"
    gene = row$Gene
    cohort = row$Cohort
#    study = row$Study
    set.seed(seed = 1)
  # out = lasso(gene = "PTEN", cohort = "Lung")  #Run CV on samples
  #  out = lasso(gene = "BRCA2", cohort = "Pancreas")  #Run CV on samples
    out = lasso(gene = gene, cohort = cohort, Cutoff = 25)  #Run CV on samples
    if(model_LOOP == "Empty"){
      Model = out[[1]]
      Features = out[[2]]
      model_LOOP = "Going"
    }else{
      Model = bind_rows(Model, out[[1]])
      Features = bind_rows(Features, out[[2]])
    }
  }, error=function(e){print(e_message)})
}

Model

write.table(Model, file = "Output/Step2_retrained_models/retrained_models.tsv", sep = "\t", col.names = T, row.names = F)
write.table(Features, file = "Output/Step2_retrained_models/Features_retrained_models.tsv", sep = "\t", col.names = T, row.names = F)

