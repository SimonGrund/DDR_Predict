# Code

### Obtaining data sets
The analysis was based on two pre-existing data sets, that of the Pan-Cancer analysis of Whole Genomes and that of the Hartwig Medical Foundation. The public parts of the PCAWG data set may be accessed at https://dcc.icgc.org/releases/PCAWG, whereas controlled files may be accessed through gbGaP and daco, as instructed on this site https://docs.icgc.org/pcawg/data/. 
The ICGC study ID of the project is EGAS00001001692. 

The HMF data used in this project may be found by accession code DR-044, and can be obtained upon request at the Hartwig Medical Foundation (https://www.hartwigmedicalfoundation.nl/en).

### Pre-processing original .vcf files for annotation of DDR gene deficiencies
In Extracting_from_vcf/ is a script for extracting DDR gene variants from the original .vcf files

### Pathogenicity annotations for analysis

- CADD Scores: CADD scoring (v1.6) may be downloaded and installed from https://cadd.gs.washington.edu/download 
- ClinVar annotation may be obtained from NCBI: https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/
- Gnomadd allele frequencies (to filter away common variants with population freq. > 0.5%) may be obtained from Broad: https://gnomad.broadinstitute.org/downloads 

### Summarizing whole-genome mutational statistics

We have provided th mutational summaries in *Supplementary Table 3 & 4* but they can be re-generated using the sig.tools.lib which can be installed from https://github.com/Nik-Zainal-Group/signature.tools.lib [1]

### Generating predictive models

The code we used to generate predictive models of DDR deficiencies may be found in /Generating_models/ . There is a seperate Readme file in this folder.

[1] Degasperi A, Amarante TD, Czarnecki J, Shooter S, Zou X, Glodzik D, Morganella S, Nanda AS, Badja C, Koh G, Momen SE, Georgakopoulos-Soares I, Dias JML, Young J, Memari Y, Davies H, Nik-Zainal S. A practical framework and online tool for mutational signature analyses show inter-tissue variation and driver dependencies. Nat Cancer. 2020 Feb;1(2):249-263. doi: 10.1038/s43018-020-0027-5. Epub 2020 Feb 17. PMID: 32118208; PMCID: PMC7048622.
### Generating data permutation models (null-models)
may be found in /Generating_models/Generating_null_models . There is a seperate Readme file in this folder.