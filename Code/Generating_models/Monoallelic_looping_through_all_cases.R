library(data.table)
library(tidyverse)
library(caret)
library(foreach)
library(doParallel)

#45887813 allll
data = fread("Data/LOF_25_all.tsv")

genes = data%>%
#  filter(Mutation_level > 1)%>%
  filter(Mutation_level > 1)%>%
  group_by(Study, Gene, primaryTumorLocation)%>%
  mutate(n = n())%>%
  ungroup()%>%
  group_by(Gene, primaryTumorLocation)%>%
  mutate(max_n= max(n, na.rm = T),
         min_n= min(n, na.rm = T)
         )%>%
  distinct(Gene, primaryTumorLocation, Study, n, max_n, min_n)

genes$nAll = genes$max_n + genes$min_n
genes2 = filter(genes, n > 10)

#genes2 = filter(genes, n >= 12)
table(genes2$Study)
source("Code/Generating_models/Lasso_regression.R")
#genes2 = genes2[1:5,]
Model = "Empty"

#Run in parraellel

output_models <- function(ModelMetrics=NULL,Model= NULL, Features = NULL)
{
  me <- list(
    ModelMetrics = ModelMetrics,
    Model = Model,
    Features = Features
  )
  
  # Set the name for the class
  class(me) <- append(class(me),"Models")
  return(me)
}

#cores=detectCores()
cl <- makeCluster(20) #Save five cores if on own machine; must be manually set to the number of cores selected on the cluster
registerDoParallel(cl)

print("starting_loop")
strt<-Sys.time()
finalMatrix <- foreach(i=1:nrow(genes2), .packages = c("tidyverse", "glmnet", "data.table", "PRROC", "cutpointr")) %dopar% {
  #do other things if you want
  # source("Feedback_from_reviewers/Split_by_dataset/Null_Model_LPC_16_06_21.R")
  row = genes2[i,]
  cat("\n\nGene=", row$Gene, "\nCohort=", row$primaryTumorLocation, "\nn=", row$n)
  tryCatch({
    e_message = "Not known"
    gene = row$Gene
    cohort = row$primaryTumorLocation
    study = row$Study
    set.seed(seed = 1)
    tempMatrix = lasso(gene = gene, cohort = cohort, trainStudy = study,  Biallelic = F)  #Run CV on samples
    out = output_models()
    out$ModelMetrics = tempMatrix[[1]]
    out$Model = tempMatrix[[2]]
    out$Features = tempMatrix[[3]]
    return(out)
  }, error=function(e){cat("Error in", gene, ", ", cohort, ", ", study)})
}

print(Sys.time()-strt)
#stop cluster
stopCluster(cl)

for(j in 1:length(finalMatrix)){
  if(j == 1){
    ModelMetrics = finalMatrix[[j]]$ModelMetrics
    Model = finalMatrix[[j]]$Model
    Features = finalMatrix[[j]]$Features
  }else{
    ModelMetrics = bind_rows(ModelMetrics, finalMatrix[[j]]$ModelMetrics)
    Model = bind_rows(Model, finalMatrix[[j]]$Model)
    Features = bind_rows(Features, finalMatrix[[j]]$Features)
  }
}

write.table(ModelMetrics, file = "Output/ModelMetrics_Monoallelic_21_06.tsv", sep = "\t", col.names = T, row.names = F)
write.table(Model, file = "Output/Models_Monoallelic_21_06.tsv", sep = "\t", col.names = T, row.names = F)
write.table(Features, file = "Output/Features_Monoallelic_21_06.tsv", sep = "\t", col.names = T, row.names = F)

print("Alll dooooone")
