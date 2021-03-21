library(data.table)
library(tidyverse)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/")

var = fread("Data/Pathogenic_variants.tsv")
#Add CT
ct = fread("Data/Cancertypes_jan2021.tsv")
ct = distinct(ct)
var = left_join(var, ct, by = c("Donor_ID"))
#Add signatures
features = fread("Data/COMBINED_featureCounts.tsv")
features = left_join(features, ct)


Variant_extract = function(gene = "BRCA2", cohort = "Breast"){
  d = filter(var, Gene %in% gene, primaryTumorLocation %in% cohort)%>%
    dplyr::select(Sample_ID,primaryTumorLocation, CADD_phred, CHROM, POS, REF, ALT, EFFECT, Variant_Type, HGVS_P, ClinVar)
  
  feats = filter(features, Sample_ID %in% d$Sample_ID)%>%distinct()
  d = left_join(d, feats)

d

}
t = Variant_extract(gene = c("MTOR"), cohort = c(unique(var$primaryTumorLocation)))
# 
# t = Variant_extract(gene = c("BRCA2", "BRCA1"), cohort = c("Bone/Soft tissue", "CNS", "Breast", "Ovary"))
# 
# #EGFR, MDM2, CDK4, NF1, ERBB2, TP53, PIK3R1, TERT 
# t = Variant_extract(gene = "TERT", cohort = "CNS")
# 
# t = Variant_extract(gene = "RAD50", cohort = "Colorectal")
# t = Variant_extract(gene = "RAD50", cohort = c("Colorectal", "Uterus", "Prostate"))
# 
# t = Variant_extract(gene = "CDK12", cohort = c("Breast", "Ovary"))
# 
# t = Variant_extract(gene = "CDK12", cohort = c("Breast", "Ovary", "Prostate", "Skin"))
# t2 = dplyr::select(t, primaryTumorLocation, CHROM, POS, REF, ALT, HGVS_P, ClinVar, Variant_Type, EFFECT, `non-clustered_inv_100Kb-1Mb`)
# 
# t = Variant_extract(gene = "ARID1A", cohort = c("Biliary","Breast", "Kidney", "Stomach", 
#                                                 "Bone/Soft tissue", "Esophagus", "NET", "Unknown"))
# t2 = dplyr::select(t, Sample_ID, primaryTumorLocation, CHROM, POS, REF, ALT, HGVS_P, ClinVar, Variant_Type, EFFECT, Signature.2, Signature.13)
# table(t2$primaryTumorLocation)
# table(t2$Variant_Type)
# d = t2
# n_distinct(d$Sample_ID) #122 variants across 106 patients
# 
# table(d$Variant_Type)
# 
# d2 = d%>%
#   #filter(Variant_Type == "Somatic SNV")%>%
#   rowwise()%>%
#   mutate(VT = paste(REF, "—>", ALT, collapse = ""))
# 
# d2 = d2[grep("Somatic SNV", d2$Variant_Type),]
# t = as.data.frame(table(d2$VT))
# t
# 
# ggplot(filter(t, Freq>1), aes(x = reorder(Var1, -Freq), y = Freq))+
#   geom_histogram(stat = "identity")+
#   theme(
#     axis.text.x = element_text(angle = 45)
#   )
# 

