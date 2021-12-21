library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)


#Load all mutation data
data = fread("Data/LOF_25_all.tsv")

#Calculated
genes = data%>%
#  filter(Mutation_level > 1)%>%
  filter(Mutation_level == 3)%>%
  group_by(Study, Gene, primaryTumorLocation)%>%
  mutate(n = n())%>%
  ungroup()%>%
  distinct(Gene, primaryTumorLocation, Study, n)

genes$nAll = genes$max_n + genes$min_n
genes2 = filter(genes, n > 5)

#genes2 = filter(genes, n >= 12)
table(genes2$Study)
#genes2 = filter(genes, nAll >= 12)
source("Code/Generating_models/Generating_null_models/Null_model_lasso_regression.R")
  
#genes2 = genes2[93,]
Model = "Empty"

#Run in parraellel
cores=detectCores()
cl <- makeCluster(55) #Save five cores if on own machine; must be manually set to the number of cores selected on the cluster
registerDoParallel(cl)

print("starting_loop")
strt<-Sys.time()
finalMatrix <- foreach(i=1:nrow(genes2), .combine=rbind, .packages = c("tidyverse", "glmnet", "data.table", "PRROC")) %dopar% {
  #do other things if you want
 # source("Feedback_from_reviewers/Split_by_dataset/Null_Model_LPC_16_06_21.R")
  row = genes2[i,]
  gene = row$Gene
  cohort = row$primaryTumorLocation
  study = row$Study
  cat("\n\nGene=", row$Gene, "\nCohort=", row$primaryTumorLocation, "\nn=", row$n)
  
  #Only calculate if not already calculated
  # tryCatch({
  #   f = fread(paste("Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_PERCASE/", gene, "_",cohort,"_", study,".tsv", sep = ""))
  # }, error=function(e){  
  # 
    tryCatch({
      e_message = "Not known"
     # set.seed(seed = 1)
      tempMatrix = lasso(gene = gene, cohort = cohort, trainStudy = study,  Biallelic = T, nPermutations = 250)  #Run CV on samples
      write.table(x = tempMatrix, 
                  paste("Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_PERCASE/", gene, "_",cohort,"_", study,"2.tsv", sep = ""),
                  sep = "\t", 
                  col.names = T, 
                  row.names = F)
      tempMatrix #Return the tmpMatrix, will be rbind with the final matrxi
    }, error=function(e){})
 # })
}
print(Sys.time()-strt)
#stop cluster
stopCluster(cl)

write.table(finalMatrix, file = "Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_Biallelic_05_07.tsv", sep = "\t", col.names = T, row.names = F)
print("Alll dooooone")
