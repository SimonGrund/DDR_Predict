library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
#library(caret)
#setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/")
#51101994
#Load all mutation data
data = fread("Data/LOF_25_all.tsv")

genes = data%>%
  #  filter(Mutation_level > 1)%>%
  filter(Mutation_level >= 2)%>%
  group_by(Study, Gene, primaryTumorLocation)%>%
  mutate(n = n())%>%
  ungroup()%>%
  distinct(Gene, primaryTumorLocation, Study, n)

genes2 = filter(genes, n > 10)

source("Code/Generating_models/Generating_null_models/Null_model_lasso_regression.R")
Model = "Empty"

#Run in parraellel
cores=detectCores()
cl <- makeCluster(cores) #Save five cores if on own machine; must be manually set to the number of cores selected on the cluster
registerDoParallel(cl)

print("starting_loop")
strt<-Sys.time()
finalMatrix <- foreach(i=1:nrow(genes2), .combine=rbind, .packages = c("tidyverse", "glmnet", "data.table", "PRROC")) %dopar% {
  #do other things if you want
 # source("Feedback_from_reviewers/Split_by_dataset/Null_Model_LPC_16_06_21.R")
  row = genes2[i,]
  cat("\n\nGene=", row$Gene, "\nCohort=", row$primaryTumorLocation, "\nn=", row$n)
  gene = row$Gene
  cohort = row$primaryTumorLocation
  study = row$Study
  #Only calculate if not already calculated
  # tryCatch({
  #   f = fread(paste("Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_PERCASE_Monoallelic/", gene, "_",cohort,"_", study,".tsv", sep = ""))
  # }, error=function(e){
  # 
    tryCatch({
      e_message = "Not known"
     # set.seed(seed = 1)
      tempMatrix = lasso(gene = gene, cohort = cohort, trainStudy = study,  Biallelic = F, nPermutations = 5e3)  #Run CV on samples
      write.table(x = tempMatrix, 
                  paste("Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_PERCASE_Monoallelic/", gene, "_",cohort,"_", study,"3.tsv", sep = ""),
                  sep = "\t", 
                  col.names = T, 
                  row.names = F)
      tempMatrix #Return the tmpMatrix, will be rbind with the final matrxi
    }, error=function(e){print("Err")})
    
 # })
}
print(Sys.time()-strt)
#stop cluster
stopCluster(cl)


write.table(finalMatrix, file = "Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_Monoallelic_21_06.tsv", sep = "\t", col.names = T, row.names = F)
print("Alll dooooone")
