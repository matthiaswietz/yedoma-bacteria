
############################################################################################
###  YEDOMA BACTERIA - AMPLICON ANALYSIS  ###
############################################################################################

# This script: Rarefaction and alpha-diversity  

setwd("/AWI_MPI/collaborations/Permafrost/Rstats")

# Load packages
library(iNEXT)
library(olsrr)
library(ggplot2)
library(cowplot)
library(dplyr)


############################################################################################
   ###  RAREFACTION AND COVERAGE  ###
############################################################################################

iNEXT <- iNEXT(
  ASV, q=c(0),
  datatype="abundance", 
  conf=0.95, nboot=50)

###################################

rarefac <- fortify(iNEXT, type=1) %>%
  left_join(ENV, by=c("Assemblage"="sample_title")) %>%
  filter(!type.convert() %in% c(
    "fSW+YE_w_mon","SW+YE_w_mon",
    "fSW+YE_s_mon","SW+YE_s_mon"))

rarefac.point <- rarefac[which(
  rarefac$Method == "Observed"),]
rarefac.line <- rarefac[which(
  rarefac$Method != "Observed"),]
rarefac.line$Method <- factor(rarefac.line$Method,
   # c("interpolated","extrapolated"),
    c("Rarefaction","Extrapolation"))

rarefaction <- ggplot(rarefac, 
  aes(x=x, y=y, colour=type)) +
  geom_line(aes(linetype = Method), 
            linewidth = 0.5, data = rarefac.line) +
  # scale_colour_discrete(guide="none") +
  # scale_colour_manual(values = rev(
  #  pnw_palette("Sunset", 34))) +
  # scale_color_manual(
  #   values = scico(77, palette='lajolla')) +
  scale_x_continuous(limits = c(0,1e+5)) +
  labs(x="Sample size", y="Species richness") +
  facet_grid(~type) +
  theme_bw(base_size=12) + 
  theme(legend.position="none",
        axis.ticks = element_blank())
  
###################################

cover <- fortify(iNEXT, type=2) %>%
  left_join(ENV, by=c("Assemblage"="sample_title")) 

cover.point <- cover [which(
  cover$Method == "Observed"),]
cover.line <- cover [which(
  cover$Method != "Observed"),]
cover.line$Method <- factor(cover.line$Method,
    #c("interpolated","extrapolated"),
    c("Rarefaction","Extrapolation"))

coverage <- ggplot(cover, 
  aes(x=x, y=y, colour=type))+ 
geom_line(aes(linetype = Method), 
  linewidth = 0.5, data = cover.line) +
scale_colour_discrete(guide="none") +
  scale_color_fish_d(option="Chaetodon_ephippium") +
scale_x_continuous(
  limits = c(0,1e+5)) +
scale_y_continuous(
  breaks = seq(0.9,1,0.05), 
  limits = c(0.9,1)) +
labs(x="Sample size", y="Sample coverage") +
  facet_grid(~type) +
theme_bw(base_size=12) + 
theme(legend.position="none",
      axis.ticks = element_blank())

# Plot curves
plot_grid(
  rarefaction, 
  coverage,
  ncol=1,
  nrow=2,
  align="tblr")

###################################

## CREATE SUMMARIES ##

richness <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Species richness",] %>%
  arrange(Assemblage) 
simpson <- iNEXT$AsyEst[
  iNEXT$AsyEst$Diversity=="Simpson diversity",] %>%
  arrange(Assemblage) 

AlphaDiv <- data.frame(
  sample_title = ENV$sample_title,
  richness = richness$Observed,
  simpson = simpson$Observed) %>%
  left_join(ENV)

###################################

# remove temp-data
rm(richness, simpson, shannon,
   rarefac, rarefac.point, rarefac.line,
   cover, cover.point, cover.line)

