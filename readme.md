## Utilization and release of fossil carbon from eroding permafrost coastlines

Bioinformaice code for bacterial community analyses in "Utilization and release of fossil carbon from eroding permafrost coastlines" by Manuel Ruben et al. (doi: XXX). Metabarcoding reads were obtained using 16S rRNA primers, followed by Illumina MiSeq sequencing and generation of amplicon sequence variants (ASVs) using DADA2 following the developer's guidelines (https://benjjneb.github.io/dada2/tutorial_1_8.html). The repo includes the full pipleline (raw fastq reads to ASV taxonomy) and all subsequent steps that result in the paper.

### Organization of directories and files

yedomaBacteria.Rmd: Rmarkdown describing primer clipping and DADA2 pipeline 

Resulting (ASV table), taxonony table, and ASV sequences

ENA accession numbers of all raw fastq files

metadata: physicochemical measurements and general sample information, needed for community analyses

Scripts for analyses/figures in the paper: loading/formatting data, rarefaction, community, and correlations

