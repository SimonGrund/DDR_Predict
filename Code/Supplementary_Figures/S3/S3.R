library(tidyverse)
library(data.table)
setwd("/Volumes/GenomeDK/PCAWG/simong_SV/DDR_project_Jan_2021/")

all = fread("Data/all_donors.tsv")
  
  ###Add features
  d = fread("Data/COMBINED_featureCounts.tsv")
  d = filter(d, Sample_ID %in% all$Sample_ID)
  d$Donor_ID = NULL;
  
  d = pivot_longer(d, cols = all_of(colnames(d)[2:ncol(d)]))

  d$"Type" = "Indels"
  d$"Type"[grep(pattern = "*Sig*", d$name)] = "SBS Signatures"
  d$"Type"[grep(pattern = "*clustered*", d$name)] = "Structural variants"
  #Subtypes
  d$"Subtype" = "Reference Signatures, Alexandrov et al.[24]"
  #d$"Subtype"[grep(pattern = "*Signature.N*", d$name)] = "Novel"
  #d$"Subtype"[grep(pattern = "*.N[[:digit:]]*", d$name, ignore.case = F)] = "Novel"
  d$"Subtype"[grep(pattern = "*.N[[:digit:]]*", d$name, ignore.case = F)] = "Degaspari et al.[14]"
  d$"Subtype"[grep(pattern = "*.MMR*", d$name)] = "Degaspari et al.[14]"
  d$"Subtype"[grep(pattern = "*PLATINUM*", d$name)] = "Degaspari et al.[14]"
  d$"Subtype"[grep(pattern = "*del.*", d$name)] = ""
  d$"Subtype"[grep(pattern = "*ins*", d$name)] = ""
  d$"Subtype"[grep(pattern = "*clustered*", d$name)] = "Clustered"
  d$"Subtype"[grep(pattern = "*non.clustered*", d$name)] = "Non-clustered"
  
  d$name =str_replace(d$name, "Signature.", "")
  d$name =str_replace(d$name, "clustered_", "")
  d$name =str_replace(d$name, "non-", "")
  
  d$col[d$Subtype == "Reference Signatures, Alexandrov et al.[24]"] = 1
  d$col[d$Subtype != "Reference Signatures, Alexandrov et al.[24]"] = 2
  
  scientific_10 <- function(x) {
    parse(text=gsub("1e\\+*", "10^", scales::scientific_format()(x)))
  }
  
  
  #Zeroes
  d2 = d%>%
    group_by(name, Type, Subtype)%>%
    mutate(
      n_zero = sum(value == 0),
      n = n(),
      prop_zero = n_zero/n()
    )%>%
    distinct(name, .keep_all = T)%>%
    dplyr::select(-c(value, Sample_ID))
  
  d3 = left_join(d, d2, by = c("name", "Type","Subtype"))
  
  d3$name = str_replace(d3$name, "_", " ")
  
  order_SNV=as.character(c(1:30, 33,36,38,51,52, "N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9", "N10", "N11", "N12", "MMR1", "MMR2", "PLATINUM"))
  order_INDELS = c("del.mh", "del.rep", "del.none","ins")
  order_SV = c("del 1-10Kb","del 10-100Kb","del 100Kb-1Mb", "del 1Mb-10Mb", "del >10Mb",
               "inv 1-10Kb","inv 10-100Kb","inv 100Kb-1Mb", "inv 1Mb-10Mb", "inv >10Mb",
               "tds 1-10Kb","tds 10-100Kb","tds 100Kb-1Mb", "tds 1Mb-10Mb", "tds >10Mb",
               "trans")
  all = c(order_SNV, order_INDELS, order_SV)
  
  d3$name = factor(d3$name, levels = all)
  
d3$Type = factor(d3$Type, levels = c("SBS Signatures", "Indels", "Structural variants"  ))

d3$col = 1
d3$col[d3$Type != "SBS Signatures"] = 2
#dtop = filter(d3, col == 1)
  #g2 <<- 
  gt1 <<-  ggplot(filter(d3, col == 1), aes(x = name, y = value))+
    geom_boxplot(color = "black", show.legend = T, outlier.shape = NA, size = 0.2, width = 0.6)+
    geom_rug(data = distinct(d3, name, .keep_all = T)%>%filter(col == 1), 
             aes(color = 1-prop_zero), sides = "b", size = 2, show.legend = F)+
    geom_jitter(width = 0.01, alpha = 0.1, size = 0.1)+
    ggh4x::facet_nested(.~Type+Subtype, scales = "free", space = "free")+
    scale_y_log10(label=scientific_10)+
    theme_bw(base_size = 9, base_family = "Arial") +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
        #  legend.position = c(0.9, 0.8),
          panel.border = element_rect(color = "black", fill = NA),
          axis.line = element_line(color = NA),
        axis.title = element_text(size = 7)
          # strip.background = element_blank(),
          # strip.placement = "outside",
        #  strip.text = element_text(color = "white"),
          # axis.text.x = element_blank(),
          # axis.ticks.x = element_blank()
    )+
    ylab("#Mutations\nper patient")+
    # xlab("Feature")+
    xlab("")+
    #geom_vline(aes(xintercept = name), col = "darkslategrey", alpha = 0.2)+
    #ggtitle("Features")+
    #scale_color_manual(values = c("black", "black"))+
    scale_fill_manual("Gene status",
                      values = c("lightgrey", "white"),
                      labels=c("Mutated", "Not mutated"))+
    scale_color_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0.5, limits = c(0,1))+
    guides(color = guide_colorbar(title = "Prop. of samples w/ 0"))
  
 #gt1
 
 gt2 <<-  ggplot(filter(d3, col == 2), aes(x = name, y = value))+
   geom_boxplot(color = "black", show.legend = T, outlier.shape = NA, size = 0.2)+
   geom_jitter(width = 0.01, alpha = 0.1, size = 0.1)+
   geom_rug(data = distinct(d3, name, Subtype, .keep_all = T)%>%filter(col == 2), 
            aes(color = 1-prop_zero), sides = "b", size = 2)+
   #   geom_jitter(width = 0.01, alpha = 0.1, size = 0.2)+
   ggh4x::facet_nested(.~Type+Subtype, scales = "free", space = "free")+
   scale_y_log10(label=scientific_10)+
   theme_bw(base_size = 9, base_family = "Arial") +
   theme(panel.grid.major = element_blank(),
         panel.grid.minor = element_blank(),
         axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
         legend.position = "bottom",
         legend.direction = "horizontal",
         legend.key.width = unit(1, "in"),
         legend.key.height = unit(0.1, "in"),
         panel.border = element_rect(color = "black", fill = NA),
         axis.line = element_line(color = NA),
         axis.title = element_text(size = 7)
         # strip.background = element_blank(),
         # strip.placement = "outside",
         #  strip.text = element_text(color = "white"),
         # axis.text.x = element_blank(),
         # axis.ticks.x = element_blank()
   )+
   ylab("#Mutations\nper patient")+
   # xlab("Feature")+
   xlab("")+
   #geom_vline(aes(xintercept = name), col = "darkslategrey", alpha = 0.2)+
   #ggtitle("Features")+
   scale_color_gradient2(low = "red", mid = "grey", high = "blue", midpoint = 0.5, limits = c(0,1))+
  # scale_color_manual(values = c("black", "black"))+
   scale_fill_manual("Gene status",
                     values = c("lightgrey", "white"),
                     labels=c("Mutated", "Not mutated"))+
   guides(color = guide_colorbar(title = "Prop. of samples with mutations", title.position="top", title.hjust = 0.5))
 
 #gt2
 
 library(patchwork)
 layout = "
 AA
 BB
"
 out = gt1 + gt2 +
   plot_annotation(
     title = "Supp. Figure 3",
     subtitle = 'a',
     caption = ''
   ) +
   patchwork::plot_layout(design = layout) & theme(text = element_text("Arial", size = 9),
                                                   title =  element_text("Arial", size = 9))
 #out¨
 
 ggsave(plot = out, "Supp_figures/S3_mutational_patterns/S3a.jpeg", height = 5.5, width = 7, device = "jpeg", dpi = 500)
 
# ggsave(plot = out, "Supp_figures/S3_mutational_patterns/S3a.pdf", height = 5.5, width = 7, device = "pdf")
 
 
 
  

