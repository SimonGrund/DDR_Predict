library(data.table)
library(tidyverse)
library(caret)

#Load the mutational status of all samples in all DDR genes, where any mutations were present
data = fread("Data/Generated dataframes/LOF_all.tsv")

#Assemble list of genes with >=12 mutated samples
genes = data%>%
  filter(LOF >=6)%>%
  group_by(GENE, primaryTumorLocation)%>%
  mutate(n = n())%>%
  ungroup()%>%
  distinct(GENE, primaryTumorLocation, n)
genes2 = filter(genes, n >= 12) #n = 11

source("Code/Analysis/Generating_models_discover_validation/Discover_and_validate_models.R")

#Loop through all the cases in which we have sufficient data
Model = "Empty"
for( row_number in 1:nrow(genes2)){
  cat(row_number, "/", nrow(genes2))
  row = genes2[row_number,]
  cat("\n\nGene=", row$GENE, "\nCohort=", row$primaryTumorLocation, "\nn=", row$n)
  tryCatch({
    e_message = "Not known"
    gene = row$GENE
    cohort = row$primaryTumorLocation
    set.seed(seed = 1)
  # out = lasso(gene = "PTEN", cohort = "Lung")  #Run CV on samples
  #  out = lasso(gene = "BRCA2", cohort = "Breast", Cutoff = 25)  #Run CV on samples
    out = lasso(gene = gene, cohort = cohort, Cutoff = 25)  #Run CV on samples
    if(Model == "Empty"){
      Model = out[[1]]
      Features = out[[2]]
    }else{
      Model = bind_rows(Model, out[[1]])
      Features = bind_rows(Features, out[[2]])
    }
  }, error=function(e){print(e_message)})
}

#See which models were made — We expect just one model, of SMAD4 in pancreas
Model

write.table(Model, file = "Results/Models.tsv", sep = "\t", col.names = T, row.names = F)
write.table(Features, file = "Results/Features.tsv", sep = "\t", col.names = T, row.names = F)

