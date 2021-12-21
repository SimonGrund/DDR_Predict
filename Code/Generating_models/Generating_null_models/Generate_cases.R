library(data.table)
library(tidyverse)
library(foreach)
library(doParallel)
#library(caret)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/GenomeDK/PCAWG/simong_SV/DDR_project_May_2021/")
# 
#Load all mutation data
data = fread("Data/LOF_25_all.tsv")

#Find biallelic cases
genes = data%>%
  filter(Mutation_level == 3)%>%
  group_by(Study, Gene, primaryTumorLocation)%>%
  mutate(n = n())%>%
  distinct(Gene, primaryTumorLocation, Study, n)%>%
  filter(n>5)

biallelic = genes
biallelic$"allelic" = "Bi"

#Fnd Monoallelic cases
genes = data%>%
  filter(Mutation_level > 1)%>%
  group_by(Study, Gene, primaryTumorLocation)%>%
  mutate(n = n())%>%
  distinct(Gene, primaryTumorLocation, Study, n)%>%
  filter(n>10)

monoallelic = genes
monoallelic$"allelic" = "Mono"

#All cases
cases = bind_rows(
  biallelic, monoallelic
  )%>%
  dplyr::select(Gene, Study, primaryTumorLocation, allelic,n )

# 
write.table(cases, "Code/Generating_models/Generating_null_models/cases.tsv", sep ="\t",
            col.names = T, row.names = F)


## New cases where no higher Nulls were found — in case you wanna upadte the list and run more cases on sum... 

# d = fread("Output(Model_performances.tsv") #Twelve fewer at 08
# cases = fread("Code/Generating_models/Generating_null_models/cases.tsv")
# 
# cases3 = filter(d, n_Null < 1e4) #Filter if any got less than 10k null models to give them another run
# cases = filter(cases, ID %in% cases3$ID)
# 
# write.table(cases, "Code/Generating_models/Generating_null_models/cases.tsv", sep ="\t",
#             col.names = T, row.names = F)

