d = fread("Data/Pathogenic_variants.tsv")
#Add CT
ct = fread("Data/Cancertypes_jan2021.tsv")
ct = distinct(ct)
d = left_join(d, ct, by = c("Donor_ID"))

#Add sig13
feats = fread("Data/COMBINED_featureCounts.tsv")
feats = left_join(feats, ct)
#feats = filter(feats, primaryTumorLocation == "Breast")
#feats = dplyr::select(feats, Sample_ID, Signature.13)

gene = "ARID1A"; cohort = c("Billiary", "Breast", "Kidney", "Stomach", "Bone/Soft tissue", 
                            "Esophagus", "NET", "Unknown"); Cutoff = 25

d = filter(d, Gene == gene, primaryTumorLocation %in% cohort)
n_distinct(d$Sample_ID) #122 variants across 106 patients

table(d$Variant_Type)

d2 = d%>%
  #filter(Variant_Type == "Somatic SNV")%>%
  rowwise()%>%
  mutate(VT = paste(REF, "—>", ALT, collapse = ""))

d2 = d2[grep("Somatic", d2$Variant_Type),]
t = as.data.frame(table(d2$VT))
t

ggplot(filter(t, Freq>1), aes(x = reorder(Var1, -Freq), y = Freq))+
  geom_histogram(stat = "identity")+
  theme(
    axis.text.x = element_text(angle = 45)
  )



