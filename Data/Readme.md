# Data for running the analysis pipeline

We have made available a series of annotation data sets from public sources, to make it easy to assemble a data set for testing the code used in the project. For the actual variants you need to download the public somatic PCAWG variants (https://dcc.icgc.org/api/v1/download?fn=/PCAWG/consensus_snv_indel/final_consensus_passonly.snv_mnv_indel.icgc.public.maf.gz)

We have assembled gnomAD and clinvar annotations fitting with this data set, and if you run the Code/Assemble_Data/generate_dataset.R it will generate the apropriate data set for the rest of the analysis.

CADD scores (filtered to fit the regions of the DDR genes) were too big to include in the project, and we have made them available for download at 

