# Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns

This is the GitHub repository for the paper "Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns"
 
# Data (*To be updated*)
It is possible to run the analysis and all figures of the paper based on somatic variants of the publicly available PCAWG samples (https://dcc.icgc.org/releases/PCAWG/). We have necessary annotations included CADD (v1.6), ClinVar annotations and gnomAD (v2.1.1) populatioon frequencies in the data folder.

## Mutational patterns
Signature contributions of single base substitutions are obtained from the supplementary information of the article "A practical framework and online tool for mutational signature analyses show inter-tissue variation and driver dependencies" by Degaspari et al. [1]. The number of indels and SVs are counted using a local installation of the Sig.tols.lib (https://github.com/Nik-Zainal-Group/signature.tools.lib). The combined set of features is made available in the Data folder (*To be uploaded*) and in th supplementary tables S3 and S4 of the paper.

# Code
We will here provide all code needed for running an example of the analysis, shown in the paper "Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns". The code is designed to run in R4.0.2. Package versions are specified in the scripts.

[1] Degasperi A, Amarante TD, Czarnecki J, Shooter S, Zou X, Glodzik D, Morganella S, Nanda AS, Badja C, Koh G, Momen SE, Georgakopoulos-Soares I, Dias JML, Young J, Memari Y, Davies H, Nik-Zainal S. A practical framework and online tool for mutational signature analyses show inter-tissue variation and driver dependencies. Nat Cancer. 2020 Feb;1(2):249-263. doi: 10.1038/s43018-020-0027-5. Epub 2020 Feb 17. PMID: 32118208; PMCID: PMC7048622.




