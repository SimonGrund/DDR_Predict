library(data.table)
library(tidyverse)
#library(maftools)
#library(extrafont)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021//")

#Patients with CDK12 in prostate
d = fread("Data/Pathogenic_variants.tsv")
#Add CT
ct = fread("Data/Cancertypes_jan2021.tsv")
ct = distinct(ct)
d = left_join(d, ct, by = c("Donor_ID"))

#4 doubles, total 14 patients as expected. 12 are HMF so no expression data :o 
d = filter(d, Gene == "CDK12", primaryTumorLocation %in% c("Prostate"))
samples = distinct(d, Sample_ID)$Sample_ID

#

d = fread("Data/Expression/joint_fpkm.tsv", nrows = 10)
