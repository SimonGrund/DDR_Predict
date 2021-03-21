# Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns

This is the GitHub repository for the paper "Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns"

# Data
It is possible to run the analysis and all figures of the paper based on somatic variants of the publicly available PCAWG samples (https://dcc.icgc.org/releases/PCAWG/). We have included ClinVar annotations and gnomAD (v2.1.1) poopulatioon frequencies in the data folder. The CADD values however, is too large a file and we suggest that users obtain them for their set of variants thorugh the scoring system at https://cadd.gs.washington.edu/

## Mutational patterns
Signature contributions of single base substitutions are obtained from the supplementary information of the article "A practical framework and online tool for mutational signature analyses show inter-tissue variation and driver dependencies" by Degaspari et al. [1]. The number of indels and SVs are counted using a local installation of the Sig.tols.lib (https://github.com/Nik-Zainal-Group/signature.tools.lib). The combined set of features is made available in the Data folder.

# Code
We here provide all code needed for running an example of the analysis, figures and supplementary figures  shown in the paper "Pan-cancer association of DNA repair deficiencies with whole genome mutational patterns". The code is designed to run in R4.0.2.

[1] Degasperi A, Amarante TD, Czarnecki J, Shooter S, Zou X, Glodzik D, Morganella S, Nanda AS, Badja C, Koh G, Momen SE, Georgakopoulos-Soares I, Dias JML, Young J, Memari Y, Davies H, Nik-Zainal S. A practical framework and online tool for mutational signature analyses show inter-tissue variation and driver dependencies. Nat Cancer. 2020 Feb;1(2):249-263. doi: 10.1038/s43018-020-0027-5. Epub 2020 Feb 17. PMID: 32118208; PMCID: PMC7048622.




