library(tidyverse)
library(data.table)
#setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/Github/")

####
# This script formats and annotates the public PCAWG variant files .maf so that it may be used as input
# for an example of the analysis.
# It is the most demanding script of the project in terms of ressources, and you will need ~20gb RAM 
# and 5 min runtime. 
#
# Simon Grund Sørensen, 21-03-2021
#
###


### Load somatic samples extracted for all 735 DDR genes, across the publicly available PCAWG samples
### First we just load the SNV_MNV_INDELS from PCAWG, which may be obtained here: https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz
PCAWG = fread("Data/download?fn=%2FPCAWG%2Fconsensus_snv_indel%2Ffinal_consensus_passonly.snv_mnv_indel.icgc.public.maf")
PCAWG = dplyr::select(PCAWG,
        CHROM = Chromosome, POS=Start_position, REF = Reference_Allele, ALT = Tumor_Seq_Allele2,
        Variant_Type, VAF = i_VAF, Donor_ID, Gene = Hugo_Symbol
                      )%>%
  filter(VAF > 0.2)

#Filter to the 735 genes of interest
GOI = readxl::read_excel("Data/Supp_Tables/Sup_tables.xlsx", sheet = 3)
PCAWG = filter(PCAWG, Gene %in% GOI$GENE)

#Filter to the samples we have included in our study after whitelisting. Only one sample per donor ID, so we can
#simply join the sample ID to the donor ID
SOI = readxl::read_excel("Data/Supp_Tables/Sup_tables.xlsx", sheet = 2)
PCAWG = filter(PCAWG, Donor_ID %in% SOI$Donor_ID)
PCAWG = left_join(PCAWG, SOI)
PCAWG = dplyr::select(PCAWG, Donor_ID, Sample_ID, primaryTumorLocation, Gene, everything())

#Write a table with all samples present in the data set
all_samples = PCAWG%>%
  distinct(Donor_ID, Sample_ID, primaryTumorLocation)

write.table(all_samples, "Data/Generated_dataframes/all_samples.tsv", sep = "\t", col.names = T, row.names = F)

#Annotate gnomAD V2.1.1 (https://cadd.gs.washington.edu/download)
gnomAD = fread("Data/gnomAD_DDR_COI.tsv")
colnames(gnomAD)= c("CHROM", "POS", "REF", "ALT", "Gnomad_AF")
gnomAD = separate(gnomAD, col = Gnomad_AF, into = c(NA, "Gnomad_AF"), sep = "=")
gnomAD$CHROM = as.character(gnomAD$CHROM)
PCAWG = left_join(PCAWG, gnomAD)

#Annotate CADD v1.6 
CADD = fread("Data/CADD16_DDR.tsv")
colnames(CADD) = c("CHROM", "POS", "REF", "ALT", "CADD_raw", "CADD_phred")
CADD$CADD_raw = NULL
PCAWG = left_join(PCAWG, CADD)
 
#Variants with CADD < 10 are omitted to make the file smaller. Let's set all of these to CADD = 5
PCAWG = mutate(PCAWG,
               CADD_phred = ifelse(is.na(CADD_phred), 5, CADD_phred))
cat("\n#NA VAF: ", sum(is.na(PCAWG$VAF)))

#Annotate ClinVar
clinvar = fread("Data/ClinVar_colsOfInterest.tsv")
PCAWG = left_join(PCAWG, clinvar)


### Let's start structuring the data.
### We want to have one row per patient per gene, giving info on mutational status (if any)

d = PCAWG
cat("Distinct PCAWG samples with variants:", n_distinct(PCAWG$Sample_ID) )
cat("\nAll variants: ", nrow(d))

#Filter away common variants (gnomAD > 0.5%, and low VAF (<0.2), and variants occuring >50 times)
d = d%>%
  filter(Gnomad_AF < 0.005|is.na(Gnomad_AF))%>%
  filter(VAF > 0.20 | is.na(VAF))%>%
  group_by(CHROM, POS, REF, ALT)%>%
  mutate(n = n())%>%
  ungroup()%>%
  filter(n < 50)%>%
  dplyr::select(-n)

#
d2 = dplyr::select(d, Donor_ID, Sample_ID, primaryTumorLocation, CHROM, POS, REF, ALT, Gene,Variant_Type, Gnomad_AF, VAF, CADD_phred, ClinVar)

#Annotate the infered state as a string, and the mutation level as a number from 1-7. 0 Will indicate no mutational change
d2=d2%>%
  mutate(
    infered_state = ifelse(ClinVar %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic, risk factor"), "Clinvar Pathogenic", "Uncertain significance"),
    infered_state = ifelse(ClinVar %in% c("Benign", "Benign/Likely benign", "Likely benign"), "Clinvar Benign", infered_state),
    infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred>=25, "Pathogenic VUS (CADD > 25)", infered_state),
    infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 20, "VUS (CADD > 20)", infered_state),
    infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 15, "VUS (CADD > 15)", infered_state),
    infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred >= 10, "VUS (CADD > 10)", infered_state),
    infered_state = ifelse(infered_state == "Uncertain significance" & CADD_phred <= 10, "Benign VUS (CADD < 10)", infered_state)
  )%>%
  mutate(Mutation_level = ifelse(infered_state == "Clinvar Pathogenic", 7, NA),
         Mutation_level = ifelse(infered_state == "Pathogenic VUS (CADD > 25)", 6, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 20)",  5, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 15)",  4, Mutation_level),
         Mutation_level = ifelse(infered_state == "VUS (CADD > 10)",  3, Mutation_level),
         Mutation_level = ifelse(infered_state == "Benign VUS (CADD < 10)",  2, Mutation_level),
         Mutation_level = ifelse(infered_state == "Clinvar Benign",  1, Mutation_level)

  )

#Remove any samples with more than 25 pathogenic evnets in DDR genes
d2 = group_by(d2, Sample_ID)%>%
  mutate(n = sum(Mutation_level >= 6))%>%
  filter(n < 25)%>%
  dplyr::select(-n)

#Save a separate sheet with pathogenic variants
patho = filter(d2, Mutation_level >= 6) #nrow = 1921
write.table(patho, "Data/Generated_dataframes/Pathogenic_variants.tsv", sep ="\t", col.names = T, row.names = F)

#Save a separate sheet with all variants
write.table(d2, "Data/Generated_dataframes/All_filtered_variants.tsv", sep ="\t", col.names = T, row.names = F) #nrow = 3859

#Some sanity checks
cat("\nPathogenic variants: ", sum(d2$Mutation_level>= 6))
cat("\nNo impact variants: ", sum(d2$Mutation_level > 3))
cat("\nClinVar annotated variants: ", sum(d2$ClinVar %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic, risk factor", "Benign", "Benign/Likely benign", "Likely benign")))
cat("\nClinVar NOT annotated variants: ", sum(!d2$ClinVar %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic, risk factor", "Benign", "Benign/Likely benign", "Likely benign")))
cat("\nClinVar benign: ", sum(d2$ClinVar %in% c("Benign", "Benign/Likely benign", "Likely benign")))
cat("\nClinVar Pathogenic: ", sum(d2$ClinVar %in% c("Pathogenic", "Pathogenic/Likely pathogenic", "Pathogenic, risk factor")))
cat("\nCADD Pathogenic variants: ", sum(d2$Mutation_level == 6))
cat("\nCADD Benign variants: ", sum(d2$Mutation_level %in% c(2,3)))
cat("\nVUS: ", sum(d2$Mutation_level %in% c(3,4,5)))

#Format the data to one line per patient per gene
d3= d2%>%
  group_by(Sample_ID, Gene)%>%
  mutate(
    Double = ifelse(sum(Mutation_level>=6)>1, T, F),
    Germline = ifelse(max(Mutation_level) >=6 && grepl("Germ", Variant_Type), T, F),
    LOF =  max(Mutation_level),
    MAX_CADD = max(CADD_phred)
  )%>%
  ungroup()%>%
  dplyr::select(Donor_ID, Sample_ID, primaryTumorLocation, Gene, Germline, LOF, Double, MAX_CADD)%>%
  distinct()%>% #Get just one row per sample_ID / gene
  group_by(Donor_ID, Sample_ID, Gene)%>%
  mutate(n= n()) #Check that all n == 1

#Get max VAF, among pathogenic somatic variants — used to calculate weights in the regression models
VAF = d2%>%
  filter(Mutation_level >= 6)%>%
  group_by(Sample_ID, Gene)%>%
  mutate(VAF = max(VAF))%>%
  dplyr::select(Sample_ID, Gene, VAF)%>%
  distinct()

d3 = left_join(d3, VAF)

#LOF indicates the mutational status. Let's see how it looks at a per-patient level
table(d3$LOF)
d3$n = NULL
d3 = rename(d3, GENE = Gene)
write.table(d3, "Data/Generated_dataframes/LOF_all.tsv", sep = "\t", col.names = T, row.names = F, quote = F)

