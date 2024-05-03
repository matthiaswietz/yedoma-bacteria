
############################################################################################
###  YEDOMA BACTERIA - AMPLICON ANALYSIS  ###
############################################################################################

# Contribution to paper "Microbial communities degrade ancient permafrost-derived organic matter in Arctic seawater"
# This script: Loading and formatting data

setwd("/AWI_MPI/collaborations/Permafrost/Rstats")

library(dplyr)
library(ggplot2)
library(tibble)
library(tidyr)
library(gtools)
library(stringr)
library(ampvis2)
library(vegan)
library(psych)
library(fishualize)


############################################################################################
   ###  LOAD AMPLICONS ###
############################################################################################

# Import ASV + taxonomy tables
ASV <- read.table(
  "seqtab.txt",
  h = T,
  sep = "\t",
  row.names = 1, 
  check.names = F)

TAX <- read.table(
  "tax.txt",
  h = T, 
  sep = "\t", 
  stringsAsFactors = F, 
  row.names = 1)

# Rename NAs with last known taxrank + "uc"
k <- ncol(TAX)-1
for (i in 2:k) {
  if (sum(is.na(TAX[, i])) >1) {
    temp <- TAX[is.na(TAX[, i]), ]
    for (j in 1:nrow(temp)) {
      if (sum(is.na(
        temp[j, i:(k+1)])) == length(temp[j, i:(k+1)])) {
        temp[j, i] <- paste(temp[j, (i-1)], "_uc", sep = "")
        temp[j, (i+1):(k+1)] <- temp[j, i]
      }
    }
    TAX[is.na(TAX[, i]), ] <- temp}
  if (sum(is.na(TAX[, i]))==1) {
    temp <- TAX[is.na(TAX[, i]), ]
    if (sum(is.na(temp[i:(k+1)])) == length(temp[i:(k+1)])) {
      temp[i] <- paste(temp[(i-1)], "_uc", sep="")
      temp[(i+1):(k+1)] <- temp[i]
    }
    TAX[is.na(TAX[, i]),] <- temp
  }
}
TAX[is.na(TAX[, (k+1)]), (k+1)] <- paste(
  TAX[is.na(TAX[, (k+1)]), k], "_uc", sep="")

# Shorten/modify names
TAX <- TAX %>%
  mutate(across(everything(),~gsub("Clade","SAR11_Clade", .))) %>%
  mutate(across(everything(),~gsub("Candidatus","Cand", .))) %>%
  mutate(across(everything(),~gsub("Roseobacter_clade_NAC11-7_lineage","Roseobacter_NAC11-7", .))) %>%
  mutate(across(everything(),~gsub("_marine_group","", .))) %>%
  mutate(across(everything(),~gsub("_terrestrial_group","", .))) %>%
  mutate(across(everything(),~gsub("_CC9902","", .))) %>%
  mutate(across(everything(),~gsub("(Marine_group_B)","", ., fixed=T))) %>%
  mutate(across(everything(),~gsub("(SAR406_clade)","SAR406", ., fixed=T)))

# Calculate mean of negative controls; round
ASV$negative <- rowMeans(ASV[,c(
  "w-negativ-25c_clip",
  "w-negativ-30c_clip")])
ASV$negative <- round(ASV$negative, 0)

# Subtract negative counts from each column
ASV = ASV - ASV$negative

# Set neg. values to zero
ASV[ASV < 0] <- 0

# Remove Mitochondria and Chloroplasts
# Remove more potential contaminants
TAX <- TAX[-grep(
  'Mitochondria|Chloroplast', TAX$Family),] %>%
  filter(!Family %in% c(
    "Corynebacteriaceae","Bacillaceae",
    "Weeksellaceae","Enterococcaceae",
    "Carnobacteriaceae","Streptococcaceae",
    "Propionibacteriaceae","Xanthobacteraceae",
    "Burkholderiaceae"))

# Match TAX after contaminant removal
ASV <- ASV[row.names(TAX),]

# Remove NegCtr columns
ASV <- ASV[, !grepl(
  'negativ', names(ASV))]


############################################################################################
   ###  LOAD METADATA  ###
############################################################################################

ENV <- read.csv(
  "metadata.txt", h=T, sep="\t", 
  stringsAsFactors=F, check.names=F) %>%
  filter(type !="negCtr")

# Sort data 
ASV <- ASV[,mixedsort(names(ASV))]
ENV <- ENV[mixedorder(ENV$clip_id),]

# Rename
names(ASV) = ENV$sample_title
rownames(ENV) = ENV$sample_title

# Sort factor levels
ENV$type <- factor(ENV$type, levels=c(
  "orgYE","orgSW","SW-postT",
  "SW+YE_s","fSW+YE_s","SW+YE_w",
  "fSW+YE_w","fSW","SW+YE_w_mon",
  "fSW+YE_w_mon","fSW_mon",
  "SW+YE_s_mon","fSW+YE_s_mon"))
ENV$days <- factor(ENV$days, levels=c(
  "orgYE","orgSW","0","4","8","12","17",
  "23","29","43","57","71","85"))


############################################################################################
   ###  FORMATTING ###
############################################################################################

# Subset to ASVs >3 counts in >3 samples
ASV.abs <- ASV %>% 
  filter(rowSums(.> 3) > 3)

# Relative abundances
ASV.rel = as.data.frame(
  apply(ASV.abs, 2, function(x) x / sum(x) * 100)) 

# Hellinger-transform
ASV.hel = as.data.frame(
  apply(ASV.rel, 2, function(x) sqrt(x / sum(x))))

# Filter TAX based on refined ASV table 
TAX <- TAX[row.names(ASV.abs),]

# Ampvis load
amp.rel <- data.frame(
  OTU = rownames(ASV.abs),
  ASV.abs, TAX, check.names=F)
amp.rel <- amp_load(amp.rel, ENV)

amp.hel <- data.frame(
  OTU = rownames(ASV.hel),
  ASV.hel, TAX, check.names=F)
amp.hel <- amp_load(amp.hel, ENV)

# Load AlphaDiv indices; merge with ENV
AlphaDiv <- read.table(
  "AlphaDiv.txt",
  h = T,  sep ="\t") %>%
  left_join(ENV)

# Define plotting colors
col = c(
  "orgYE"="gray8",
  "orgSW"="mediumorchid1", #deepskyblue2
  "SW-postT"="mediumorchid3",
  "SW+YE_s"="lightgoldenrod2", #tan3
  "fSW+YE_s"="goldenrod3", #darkgoldenrod1
  "SW+YE_w"="lightseagreen", #cyan3
  "fSW+YE_w"="cadetblue2", #aquamarine cyan1
  "fSW"="azure3",
  "fSW_mon"="azure4",
  "SW+YE_w_mon"="azure4",
  "fSW+YE_w_mon"="azure4",
  "SW+YE_s_mon"="gray44",
  "fSW+YE_s_mon"="gray44")

# Remove temporary data
rm(temp,i,j,k)

# Save workspace
save.image("Yedoma.Rdata")
