library(tidyverse)
library(data.table)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021/")
source("Code/R_analysis/Lasso_regression/Helper_functions/sumlog.R")

Features = fread("Results/Features_30Fold.tsv")
Features = Features%>%
  rowwise()%>%
  mutate(
    ID =  paste(Gene,"|",Cohort, sep =""),
    direction = ifelse(t_val > 0,"Positive", "Negative"),
  )%>%
  group_by(ID)%>%
  mutate(max_std = max(stderror))%>%
  ungroup()%>%
  filter(max_std < 20)
  
Models = fread("Results/Models_30Folds.tsv") 
Models <<- Models%>%
  rowwise()%>%
  mutate(
    ID =  paste(Gene,"|",Cohort, sep =""),
    p_value_combined = sumlog(c(p_val_study, p_val_test))$p,
    q_val = p.adjust(p_value_combined, method = "fdr", n = 265)
  )%>%
  filter(ID %in% Features$ID)%>%
  distinct(ID, .keep_all = T)%>%
  mutate(
    col = ifelse(ROC_test>=0.7 & ROC >= 0.7, "ROC > 0.7", "ROC < 0.7"),
    label = ifelse(ROC_test >= 0.7 & ROC >= 0.7 , paste(ID,"\np=", round(p_value_combined, 5),sep = ""), NA)
  )

Features = filter(Features, ID %in% filter(Models, !is.na(label))$ID)

write.table(Models, file = "Figures/Figure_3_new//models.tsv", sep ="\t", col.names = T, row.names = F)
write.table(Features, file = "Figures/Figure_3_new//features.tsv", sep ="\t", col.names = T, row.names = F)


  

# ggplot(filter(Features, ID %in% filter(Models, !is.na(label))$ID), aes(x= Estimate, y = Feature))+
#   geom_errorbar(aes(xmin = Estimate - stderror, xmax = Estimate + stderror))+
#   facet_wrap(~Gene+Cohort, scales = "free")
# 
# ggplot(Models, aes(x = ROC_test, y = -log(p_value_combined)))+
#   geom_point(aes(color = col))+
#   scale_color_manual(values = c("grey", "darkblue", "brown"))+
#   geom_text(aes(label = label), hjust = 1, vjust = 1)+
#   geom_hline(yintercept = -log(0.05), lty = 2)+
#   geom_vline(xintercept = 0.7, lty = 2, col = " grey")+
#   theme_bw(base_size = 8, base_family = "Arial") +
#   theme(panel.grid.major = element_blank(),
#         panel.grid.minor = element_blank(),
#         legend.position = c(0.9, 0.1),
#         panel.border = element_rect(color = "black", fill = NA),
#         axis.line = element_line(color = NA))+
#   guides(color = guide_legend(title =""))

#write.table(Features, "Results/Features.tsv", sep = "\t", col.names = T, row.names = F)

#Models = filter(Models, ID %in% Features$ID)
#length(unique(Models$Gene))
