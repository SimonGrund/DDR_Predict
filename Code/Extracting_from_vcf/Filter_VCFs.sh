#!/bin/bash
#SBATCH --partition normal
#SBATCH --mem-per-cpu 40G
#SBATCH -c 1
#SBATCH -t 40:00:00

#################################################################
###                                                           ###
### Filter, combine and snpEff-annotate a folder of VCF files ###
### Requires a .bed file of regions of interest, vcftools     ###
### as well as snpEff and snpSift jar-files                   ###
### and an environment with R3.5 installed                    ###
### June 2020, Simon Grund Sørensen                           ###
###                                                           ###
#################################################################

################################
###Filter all VCF files in folder and export the filtered VCFs into a folder called snpSift_filtered
###############################

#Run script in folder of .vcf files.
#Have an folder called "software" with the snpEff files from https://pcingola.github.io/SnpEff/download/
#Have a folder called "regions" with a bed file of regions of interest
mkdir filtered_VCF

#Filter to the set of variants using snpSift. Alternatively use tabix.
for F in *.vcf.gz; do \
    s=$F; \

zcat $F | \
java -jar software/snpEff_latest_core/snpEff/SnpSift.jar \
intervals regions/intervalList_allChrom.bed \
>  filtered_VCF/$s;
done

cd filtered_VCF/


###################
### Add filename and INFO column to the ID column (sep = |), 
### and merge all files into one combined table 
###################
mkdir combined_table

header="CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO"
echo -e "$header" > barcoded/combined.vcf; \
for F in *vcf; do \
    s=$F; \
    echo  "$header" > combined_table/combined.vcf; \
    grep -v "#" ${F} | awk -F "\t" -v var="$s" '{print $1"\t"$2"\t"$3"|"var"|"$8"\t"$4"\t"$5"\t"$6"\t"$7"\t" }' >> combined_table/combined.vcf; \
done

cd combined_table/

###########################
### snpEff on the combined file to get snpEff annotation. Not strictly necessary
###########################
java -jar /home/simong/PCAWG/simong_SV/software/snpEff_latest_core/snpEff/snpEff.jar \
-v GRCh37.75 \
-noStats \
-canon  \
-onlyProtein  \
-t \
-hgvs -lof -no-downstream -no-upstream -no-intergenic \
combined.vcf \
> combined.ann.vcf

###################################
### Extract fields of interest: genename, HGVS_p, "ANN[*].EFFECT"  "ANN[0].AA_POS" "ANN[0].CDS_POS"
###################################
 java -jar /home/simong/PCAWG/simong_SV/software/snpEff_latest_core/snpEff/SnpSift.jar \
    extractFields -e "NA" \
    combined_filtered.vcf \
    CHROM POS REF ALT ID QUAL FILTER "ANN[0].GENE" "ANN[0].HGVS_P" "ANN[0].EFFECT" "ANN[0].AA_POS" "ANN[0].CDS_POS"\
    > combined.ann.filtered.vcf

###################################
#### Use the attached R-code to split ID into ID, filename and INFO
#### INFO is then split into key-value pairs.
#### Also annotates CADD score and GNOMAD score
###################################

source activate R35
Rscript software/split_ID_column.R

#####################
### The output table with all variants in the regions of interest and the sample barcode is in 
### filtered_VCF/combined_table/DDR_variants.table
#####################
rm combined* #Clean up


