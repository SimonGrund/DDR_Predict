

resultPlot = function(gene = "BRCA2", cohort = "Breast", study = "HMF"){
#setwd("Article/Figures/Fig4/")


#gene = "BRCA2"; cohort = "Breast"; study = "HMF"; Study = "HMF"
#gene = "RB1"; cohort = "Lung"; study = "HMF"; Study = "HMF"
#  gene = "ARID1A"; cohort = "Esophagus"; study = "HMF"; Study = "HMF"
####
# Lollipop
####
source("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/lollipop.R")
l = lollipoppy(gene = gene, cohort = cohort, Study = study, labels = F)
l
l = l + labs(title = "",
         subtitle = paste(gene, "mutations in", cohort, "cancer patients in the", study, "study", "( N =",n_FG,")"))+
  theme(
  legend.justification = "bottom"
)
l
####
# Features
####
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/")
source("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/feature_plot.R")
feature_plot(gene = gene, cohort = cohort, Study = study)

ftop = g1
ftop = ftop+labs(title = "",
                 ssubtitle = "All features")
                # subtitle = paste("Features counts of", cohort, "cancer patients in the", study, "study"))
fbot = g2+
  theme(
    legend.justification = "bottom"
  )



####
# Results Models
####
#gene = c("BRCA1", "CDK12", "BAP1")

setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/")
source("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/results.R")

#set.seed(2)
lasso(gene = gene, cohort = cohort, study = study)
Posterior_probabilities
#dfeat_out
#ROC_curve
####
# Results Features
####
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig4/")
source("SimulateFeatureInfluence.R")

featureImportance = FeatureImportance(gene = gene, cohort = cohort, study = study)
featureImportance = featureImportance+ggtitle(label = "", subtitle = "Features from modelling")+
  theme(
    legend.justification = "bottom"
  )
###Combine figure

###
# Features and mutations
###

library(patchwork)
#ROC = ROC + theme(legend.position="none")

# space = ggplot()+theme_void() + labs(title = "C",
#                       subtitle = paste("Model performance and features from model trained and LOOCV validation on", study, " data"))
# 

#bottom = space/(Posterior_probabilities + ROC_curve + dfeat_out + featureImportance+ patchwork::plot_layout(widths =c(2, 1, 2, 1))) + patchwork::plot_layout(heights = c(0.2,1))
bottom =  (dfeat_out + featureImportance) / (Posterior_probabilities + ROC_curve) + 
  patchwork::plot_layout(widths =c(3, 1, 3, 1))

bottom = bottom & theme(
  legend.position = "bottom"
)
#bottom
# 
# bottom = space/(Posterior_probabilities + ROC_curve) / (dfeat_out + featureImportance) + patchwork::plot_layout(widths =c(3, 1, 3, 1)) + patchwork::plot_layout(heights = c(0.2,1,1))
# bottom &theme(
#   legend.position = "bottom"
# )& ggtitle(paste(
#   paste(gene, collapse = "|"), 
#   ", ", 
#   paste(cohort, collapse = "|"))
# )



Input = l / ftop/fbot / bottom & 
  theme(legend.position = "right") & 
  guides(color = guide_legend(nrow = 6, title.position="top", title.hjust = 0.5),
         fill = guide_legend(nrow = 6, title.position="top", title.hjust = 0.5))
  
out = Input +   plot_annotation(tag_levels = 'A')  + patchwork::plot_layout(heights = c(1,0.7,0.3,1,1))
out

layout <- "
AAAAAA
AAAAAA
BBBBBB
BBBBBB
CCCCCC
DDDDGG
DDDDGG
EEFFGG
EEFFGG
EEFFGG
"
out = l + ftop + fbot + dfeat_out + Posterior_probabilities + ROC_curve + featureImportance + 
  patchwork::plot_layout(design = layout, guides = "collect") +
  plot_annotation(tag_levels = 'A', theme = theme(legend.position = 'bottom') )&
  theme(
    axis.text = element_text(size = 12),
    axis.title = element_text(size = 14),
    plot.subtitle = element_text(size = 14),
    plot.title = element_text(size = 18),
    strip.text = element_text(size = 14),
    legend.text = element_text(size = 12),
    legend.title = element_text(size = 14)
  )


#out = Input + patchwork::plot_layout(heights = c(1,0.2,1,0.3,3))
#out
                                     #guides = "collect")

#ggsave(plot = out, "/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures_new/Figures/Fig3.pdf",
#       device = "pdf", units = "in", width = 18, height = 18)



ggsave(plot = out,
       paste("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Plots/ResultFigures/",gene, cohort, study, "_Fig4.pdf", sep = ""),
       device = "pdf", units = "in", width = 14, height =20)

}

resultPlot("BRCA2","Breast","HMF")
#resultPlot()
# 
#resultPlot(gene = "RB1", "Lung", "HMF")
# #
# 
# resultPlot(gene = "RB1", "Lung", "HMF")
# 
# #setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/")
# #source("Article/Figures/Fig3/Filter_results.R")
# source("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig3/Filter_results.R")
# #Models = filter(Models, ID %in% Features$ID)   #Only look at models which have any significant features
# #Models = Models[grep("\\|", Models$Gene, invert = T),] #Fix it so i can look on groups of genes as well
# #Models = Models[1,]
source("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_June_2020/Article/Figures/Fig3/Filter_results.R")

for(i in 1:nrow(Models)){
  resultPlot(gene = Models$Gene[i], cohort = Models$Cohort[i], study = Models$Study[i])
}

