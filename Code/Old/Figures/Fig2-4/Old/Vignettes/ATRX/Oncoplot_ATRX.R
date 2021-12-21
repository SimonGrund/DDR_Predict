######
# Figure 1 — Summary figure 
# - 1. Overview of samples per study, cancertype and metastasis
# - 2. Overview of DDR deficiency across samples
# - 3. Overview of types of mutational signals
######
library(data.table)
library(tidyverse)
library(maftools)
library(extrafont)
setwd("/Volumes/Macintosh HD/Users/au460892/Desktop/R-projects/DDR_project_Jan_2021//")
source("Code/Plotting_functions//custom_onco.R")
source("Code/Plotting_functions/maftoolsHelperFuncs.R")

#####
# Make maftools plot
#####

d_all = fread("Data/Pathogenic_variants.tsv")

###Load genes of interest

d = d_all%>%
  separate(col = EFFECT, sep = "&", into = c("EFFECT", NA))

d = filter(d, VAF>0.20| is.na(VAF))
d = filter(d, Gnomad_AF < 0.005|is.na(Gnomad_AF))

#Add CT
ct = fread("Data/Cancertypes_jan2021.tsv")
ct = distinct(ct)
d = left_join(d, ct, by = c("Donor_ID"))

#d$EFFECT[d$EFFECT %in% c("protein_protein_contact", "structural_interaction_variant", "sequence_feature","synonymous_variant", "intron_variant", "intragenic_variant")] = "Other"
d$EFFECT[d$EFFECT %in% c("splice_acceptor_variant", "splice_donor_variant", "splice_region_variant")] = "splice site variant"
d$EFFECT[d$EFFECT %in% c("frameshift_variant", "disruptive_inframe_insertion")] = "Frameshift variant"
d$EFFECT[d$EFFECT %in% c("start_lost", "stop_lost", "stop_gained")] = "Start/stop variant"
d$EFFECT[d$EFFECT %in% c("3_prime_UTR_variant", "5_prime_UTR_variant")] = "3'/5' UTR variant"

table(d$EFFECT)

d$EFFECT = str_replace_all(d$EFFECT, "_", " ")
table(d$EFFECT)

d = separate(d, Variant_Type, into = c("state", "Variant_Type"), sep = " ")

#########
d = dplyr::select(d, Hugo_Symbol = Gene, Chromosome = CHROM, Start_Position = POS, End_Position = POS, Reference_Allele = REF, 
                  Alternative_Allele = ALT, Tumor_Sample_Barcode = Sample_ID, Variant_Classification = EFFECT, 
                  Variant_Type = Variant_Type, CADD_phred, cancertype = primaryTumorLocation, Study, state,VAF, HGVS_P)

d$Tumor_Seq_Allele1 = d$Reference_Allele
d$Tumor_Seq_Allele2 = d$Alternative_Allele

d3 = d
#d$Variant_Classification = as.factor(d$Variant_Classification)
levels(as.factor(d3$Variant_Classification))

n_distinct(d3$Tumor_Sample_Barcode)


#Add missing donors
msi = fread("Data/all_donors.tsv")
msi = filter(msi, !Sample_ID %in% d3$Tumor_Sample_Barcode)%>%dplyr::select(Tumor_Sample_Barcode = Sample_ID, everything())
d4 = bind_rows(d3, msi)
n_distinct(d4$Tumor_Sample_Barcode)

# #d4  = d3
# genes = fread("Figures/Figure_3_new/Models_booty.tsv")
# g = distinct(genes, Gene)$Gene
# 
# d5 = dplyr::filter(d4, Hugo_Symbol %in% g)
d5 = filter(d4, cancertype %in% c("Breast", "Bone/Soft tissue", "Ovary", "CNS"))
laml.maf <<- read.maf(
  maf = d5,
  clinicalData = d4,
  vc_nonSyn = levels(as.factor(d4$Variant_Classification)),
) 
d5 = filter(d4, cancertype %in% c("CNS"))
laml.maf2 <<- read.maf(
  maf = d5,
  clinicalData = d4,
  vc_nonSyn = levels(as.factor(d4$Variant_Classification)),
) 

#pathways = filter(pathways, Pathway %in% c("Mismatch Repair (MMR)", "Damage Sensor etc.",))

aml_genes_vaf = subsetMaf(maf = laml.maf, genes = distinct(d4, Hugo_Symbol)$Hugo_Symbol, fields = "VAF", mafObj = FALSE)[,mean(VAF, na.rm = TRUE), Hugo_Symbol]
#genes = fread("Figures/Figure_4_new/Models_performance_allCohorts.tsv")


dev.off()
#grDevices::cairo_pdf(filename = "Article/Figures_new/Figure_2/Fig2_D_oncoplot.pdf", width = 12, height = 6,  family = "Liberation Sans")
#png(file = "Article/Figures/Fig1/Oncoplot.png",width = 210, height = 120, units = "mm", res = 180)
pdf("Figures/Figure_4_new/Vignettes/ATRX/oncoplot_ATRX.pdf", 
    width = 7.5, height = 3.5,  
    useDingbats=FALSE, family = "Arial", pointsize = 9)

 oncoplot(laml.maf, top = 25,
          genes = c("ATRX", "IDH1", "BRCA2", "BRCA1"),
          #sepwd_samples = -1,
        #  clinicalFeatures = c('cancertype', 'Study'), 
          sortByAnnotation = F, sortByMutation = T,
          gene_mar = 20, colors = yarrr::piratepal(palette = "xmen"),
          bgCol = "grey90",
        #  pathways = 'auto',
          #pathways = pathways, ,sampleOrder = sOrder$Tumor_Sample_Barcode,
          additionalFeature = c("state", "Germline")
       #selectedPathways = c("MMR","HR")
     ##  leftBarData = aml_genes_vaf,
       #leftBarLims = c(0, 1)
          )
    #      bgCol = NA, 
     #     borderCol = "grey60", 
      #    gene_mar = 6, 
       #   sepwd_genes = 1 
      #    annotationColor = colors_study, 
         # colors = colors_features
 dev.off()

#Lollipoppy
source("Code/Plotting_functions/customLolly.R")
#lollipopPlot(laml.maf, gene = "ATRX", AACol = "HGVS_P") 
 
###Somatic interactions
somaticInteractions(maf = laml.maf2, 
                    #top = 10,
                    genes = c("BRCA2", "BRCA1", "ATRX", "IDH1", "TP53", "PBRM1"),
                  #  top = 10, 
                     pvalue = c(0.05, 0.1), 
                    showCounts = T, 
                    showSigSymbols = F)


#lollipopPlot(m = laml.maf, gene = "ATRX", AACol = HGVS_P)


