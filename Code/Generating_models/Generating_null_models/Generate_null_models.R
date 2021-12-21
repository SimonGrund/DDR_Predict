suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(tidyverse)))
suppressWarnings(suppressMessages(library(foreach)))
suppressWarnings(suppressMessages(library(doParallel)))
#library(caret)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/Feedback_from_reviewers/Split_by_dataset/Models/Null_models_10k/")
#Get command Args
args = commandArgs(trailingOnly=T)
args = as.integer(as.character(args))
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

task_nr = args[1]

cases = fread("Code/Generating_models/Generating_null_models/cases.tsv") #A file with one line per case you want to null-model

case = cases[task_nr,]
print(case)
print("Done")

#Load the functions needed
source("Code/Generating_models/Generating_null_models/Null_model_lasso_regression.R")


#Run in paraellel
cores=detectCores()
cl <- makeCluster(cores) #Save five cores if on own machine; must be manually set to the number of cores selected on the cluster
registerDoParallel(cl)

#Start the loop — loop through the same case 1k times for each core = 10k results
strt<-Sys.time()
finalMatrix <- foreach(i=1:cores, .combine=rbind, .packages = c("tidyverse", "glmnet", "data.table", "PRROC")) %dopar% {
  #do other things if you want
 # source("Feedback_from_reviewers/Split_by_dataset/Null_Model_LPC_16_06_21.R")
  #row = genes2[i,]
  row = case
  gene = row$Gene
  cohort = row$primaryTumorLocation
  study = row$Study
  allellic = ifelse(case$allelic == "Bi", T, F)
  cat("\n\nGene=", row$Gene, "\nCohort=", row$primaryTumorLocation, "\nn=", row$n)

  #Only calculate if not already calculated
  # tryCatch({
  #   f = fread(paste("Feedback_from_reviewers/Split_by_dataset/Models/Null_Models_PERCASE/", gene, "_",cohort,"_", study,".tsv", sep = ""))
  # }, error=function(e){
  #
    tryCatch({
      e_message = "Not known"
     # set.seed(seed = 1)
      tempMatrix = lasso(gene = gene, cohort = cohort, trainStudy = study,  Biallelic = allellic, nPermutations = 100)  #Run CV on samples
      write.table(tempMatrix, paste0("Output/null_models/",
                                     case$Gene[1],"_",
                                     case$Study[1],"_",
                                     str_replace_all(case$primaryTumorLocation[1], "\\/|\\ ","_") ,"_",
                                     case$allelic, "_",
                                     task_nr,"_",
                                     str_replace_all(str_replace_all(as.character(Sys.time()), ":|-", "")," ", "_"), #fix time
                                     ".tsv"), 
                  sep = "\t", col.names = T, row.names = F)
      
      tempMatrix #Return the tmpMatrix, will be rbind with the final matrxi
    }, error=function(e){})
 # })
}
print(Sys.time()-strt)
#stop cluster
stopCluster(cl)

write.table(finalMatrix, paste0("Output/null_models/Collected_null_models/",
                                case$Gene[1],"_",
                                case$Study[1],"_",
                                str_replace_all(case$primaryTumorLocation[1], "\\/|\\ ","_") ,"_",
                                case$allelic, "_",
                                task_nr,"_",
                                str_replace_all(str_replace_all(as.character(Sys.time()), ":|-", "")," ", "_"), #fix time
                               ".tsv"), 
            sep = "\t", col.names = T, row.names = F)

