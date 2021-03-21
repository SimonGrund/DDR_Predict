d = fread("Data/Pathogenic_variants.tsv")
#Add CT
ct = fread("Data/Cancertypes_jan2021.tsv")
ct = distinct(ct)
d = left_join(d, ct, by = c("Donor_ID"))

#Add sig13
feats = fread("Data/COMBINED_featureCounts.tsv")
feats = left_join(feats, ct)
feats = filter(feats, primaryTumorLocation == "Breast")
feats = dplyr::select(feats, Sample_ID, Signature.13)

d = filter(d, Gene == "CDK12", primaryTumorLocation %in% c("Breast"))
d = left_join(d, feats)


#3 patients have C -> CT : Are they APOBEC?

p = c("CPCT02220029T", "DRUP01220001T", "DRUP01220003sT")
p = d$Sample_ID

all = fread("Data/all_donors.tsv")
t = filter(all, Sample_ID %in% p)
t$Study = NULL; t$primaryTumorLocation = NULL; t$Sig = T

feats = fread("Data/COMBINED_featureCounts.tsv")
feats = left_join(feats, ct)

t = left_join(feats, t)

ggplot(t, aes(x = "Sig. 13", y = Signature.13))+
  geom_boxplot()+
  geom_point(aes(col = Sig))
