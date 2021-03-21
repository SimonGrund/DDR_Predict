library(tidyverse)
library(data.table)

d = fread("combined.ann.filtered.vcf", skip = 0, sep = "\t", header = T)
d = rename(d, HGVS_P = 'ANN[0].HGVS_P', GENE = 'ANN[0].GENE', EFFECT = "ANN[0].EFFECT", AA_POS = "ANN[0].AA_POS", CDS_POS = "ANN[0].CDS_POS")
d = separate(d, col = ID, into = c("ID","Tumor_sample_barcode", "INFO"), sep = "\\|") #Double check this fits to the labelling of your files
print("Loaded combined.ann.filtered.vcf")

KVsep <- fixed(";")  #key-value separator
Vsep <- fixed("=")     #value separator

d2 <-  d %>%
  mutate(KVpairs = str_split(INFO, KVsep)) %>%
  unnest(KVpairs) %>%
  separate(KVpairs, into = c("key", "value"), Vsep) %>% #Split INFO
  spread(key, value) %>%
  dplyr::select(-INFO)%>%
  select_if(~sum(!is.na(.)) > 0)

write.table(d2, "DDR_variants.table",sep = "\t", col.names = T, row.names = F, quote = F)
