## Release and microbial utilization of fossil carbon from eroding permafrost coastlines

Bioinformatic code for bacterial community analyses in "Release and microbial utilization of fossil carbon from eroding permafrost coastlines" by Ruben and colleagues (doi: XXX). Metabarcoding reads were obtained using 16S rRNA primers, followed by Illumina MiSeq sequencing and generation of amplicon sequence variants (ASVs) using [DADA2](https://benjjneb.github.io/dada2/tutorial_1_8.html). The repo includes the full pipeline (raw fastq reads to ASV taxonomy) and all subsequent analyses.

### Overview of content

[Rmarkdown](yedomaBacteria.Rmd) describing primer clipping and DADA2 pipeline 

Resulting [ASV table](seqtab.txt), [taxonomy table](tax.txt), and [ASV sequences](ASV.fasta)

[ENA accession numbers](ENA_accessions.txt) of all raw fastq files

[Physicochemical data and sample info](metadata.txt), needed for community analyses

Scripts for analyses/figures, to be run in this order: [load/format data](DataLoad.R), [calculate alpha-diversity](RarefacDiversity.R), [community analyses](Res_Communities.R), [environmental correlations](Res_Correlations.R)

