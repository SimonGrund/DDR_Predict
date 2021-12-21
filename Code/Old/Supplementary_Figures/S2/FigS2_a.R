library(tidyverse)
library(data.table)
library(extrafont)
setwd("/Volumes/Macintosh HD/Users/au460892/Documents/GitHub/DDR_Predict/")

d = fread("Data/Generated_dataframes/All_filtered_variants.tsv")

#Split variants by the annotation
d2 = d%>%
  mutate(
    subgroup = "VUS",
    subgroup = ifelse(Mutation_level < 3, "Benign", subgroup),
    subgroup = ifelse(Mutation_level > 5, "Pathogenic", subgroup)
  )
d2$subgroup = factor(d2$subgroup, levels = c("Benign", "VUS", "Pathogenic"))
d2$state = "Somatic"

d3 = d2%>%
  group_by(Variant_Type, subgroup, state)%>%
  mutate(
    n = n()
  )%>%
  distinct(Variant_Type, subgroup, state, n)

#Add scientific labelling to plot
scientific <- function(x) {
  parse(text=gsub("1e\\+*", "10^", scales::scientific_format()(x)))
}

all = fread("Data/Generated_dataframes/all_samples.tsv")

gB <<- ggplot(d3, aes(x = subgroup, y = n,fill = state))+
  geom_histogram(stat = "identity", position = "dodge", width = 0.6, color = "black", size = 0.2)+
  geom_text(aes(label = paste("n=",n,"  ",sep =""), angle = 90), 
            position = position_dodge(width = 0.6), hjust = 1, size = 2.5, family = "Arial")+
  theme_bw(base_size = 9, base_family = "Arial") +
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      #  legend.position = c(0.92, 0.9),
        panel.border = element_rect(color = "black", fill = NA),
        legend.background = element_rect(fill =NA),
        axis.line = element_line(color = NA))+
  guides(fill = guide_legend(title = ""))+
  ggh4x::facet_nested(~Variant_Type, scales = "free_x")+
  ylab("#Variants")+
  #ggtitle("Variants across 735 DDR genes in 5,990 patients")+
  xlab("")+
  scale_y_log10(label=scientific)+
  scale_fill_manual(values= c("Goldenrod2", "darkgreen"))+
  labs(color  = "")

gB

####Part B: This part is not generated as we have not included SNPeff or any other effect predictor in the data, 
#### and it isn't used in the analysis.
#
# db = d2 %>%
#   mutate(ClinVar = ifelse(ClinVar %in% c("Benign", "Benign/Likely benign", "Likely benign"), "Benign", ClinVar),
#          ClinVar = ifelse(ClinVar %in% c("Pathogenic, risk factor ","Pathogenic","Likely pathogenic", "Pathogenic/Likely pathogenic", "Likely pathogenic, drug response"), "Pathogenic", ClinVar),
#          ClinVar = ifelse(ClinVar %in% c("Benign", "Pathogenic"), ClinVar,  "Uncertain")
#   )
# 
# ggplot(db, aes(x = reorder(EFFECT, med), y = CADD_phred, fill = state))+
#   geom_boxplot(outlier.shape =  NA, size = 0.2)+
#   geom_point(size = 0.3, position = position_jitterdodge(jitter.width = 0.1), aes(col = ClinVar, group = state))+
#   geom_hline(yintercept = 25, lty = 2, col = "grey")+
#   geom_hline(yintercept = 10, lty = 2, col = "grey")+
#   theme_bw(base_size = 9, base_family = "Arial") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
#         #        legend.position = c(0.92, 0.9),
#         panel.border = element_rect(color = "black", fill = NA),
#         legend.background = element_rect(fill =NA),
#         axis.title = element_text(size = 7),
#         axis.line = element_line(color = NA))+
#   scale_y_continuous(breaks = c(0, 10,  25,  40), 
#                      labels = c(0, 10,  25, 40))+
#   xlab("snpEff predicted variant effect")+
#   ylab("CADD-phred score")+
#   scale_fill_manual(values= c("Goldenrod2", "darkgreen"))+
#   scale_color_manual(values = c("grey", "cornflowerblue", "brown"))
# 
# 
# 
# 


ggsave(plot = gB, device = "pdf", "Output/Supp_Figures/S2A.df", width = 5.5, height = 3)

