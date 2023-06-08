############################################################################################
###  YEDOMA BACTERIA - AMPLICON ANALYSIS  ###
############################################################################################

## This script: community ordination + composition
# Identify origin per Genus (max-abundance in SW or YE)

tax <- TAX %>%
  rownames_to_column("asv") %>%
  distinct(Genus, .keep_all=T)

oriGen <- ASV.rel %>% 
  rownames_to_column("asv")   %>% 
  reshape2::melt() %>% 
  left_join(ENV, by=c("variable"="sample_title")) %>%
  left_join(tax)  %>%
  drop_na(Genus) %>%
  filter(type %in% c("orgYE","orgSW","SW-postT")) %>%
  group_by(Genus, type) %>%
  summarize(maxAb=mean(value)) %>%
  filter(maxAb == max(maxAb)) %>%
  ungroup() %>%
  mutate(type=case_when(
    grepl("SW", type)~"seawater", 
    TRUE~"permafrost")) 

# Identify origin per ASV (max-abundance in SW or YE)
oriASV <- ASV.rel %>% 
  rownames_to_column("asv")   %>% 
  reshape2::melt() %>% 
  left_join(ENV, by=c("variable"="sample_title")) %>%
  left_join(tax)  %>%
  drop_na(Genus) %>%
  filter(type %in% c("orgYE","orgSW","SW-postT")) %>%
  group_by(asv, type) %>%
  summarize(maxAb=mean(value)) %>%
  filter(maxAb == max(maxAb)) %>%
  ungroup() %>%
  mutate(type=case_when(
    grepl("SW", type)~"seawater", 
    TRUE~"permafrost"))  %>%
  dplyr::select(c(asv, type, maxAb)) 

###################################################

## Broad overview -- NMDS
# only days sampled in all treatments
# export size 4x6
amp_ordinate(amp_subset_samples(
  amp.rel, type %in% c(
    "orgYE","orgSW","SW-postT",
    "SW+YE_s","fSW+YE_s",
    "SW+YE_w","fSW+YE_w") & !days %in% c(
    "4","12","29","43","71")),
  type = "PCoA",
  distmeasure = "bray",
  transform = "hellinger", 
  sample_color_by = "type", 
  sample_shape_by = "period",
  sample_point_size = 4) +
  scale_color_manual(values=col) +
  theme_classic() + 
  theme(#axis.text = element_blank(),
        #axis.title = element_blank(),
        axis.ticks = element_blank(),
        #panel.border = element_blank(),
        panel.grid.minor = element_blank(),
        legend.position = "right")

###################################################

## Broad overview -- Order level
# export size 4x5
amp_heatmap(amp_subset_samples(
  amp.hel, type %in% c(
    "orgYE","orgSW","SW-postT","SW+YE_s",
    "fSW+YE_s","SW+YE_w","fSW+YE_w") &
    !days %in% c("4","12","29","43","71")),
  tax_aggregate = "Order",
  order_x_by=c(
    "orgSW","SW-postT","orgYE","SW+YE_s",
    "fSW+YE_s","SW+YE_w","fSW+YE_w"),
  order_y_by=c(
    "SAR11 clade","Pseudomonadales",
    "Flavobacteriales","Enterobacterales",
    "Rhodobacterales","Pirellulales",
    "Burkholderiales","Microtrichales",
    "Micrococcales","Thermomicrobiales",
    "KD4-96_uc","Gitt-GS-136_uc",
    "Gaiellales","Rhizobiales",
    "Solirubrobacterales"),
  tax_show = 15,
  normalise = F,
  plot_values = F,
  group_by = "type")  +
  scale_fill_fish(
    option="Ostorhinchus_angustatus",
    begin = 0, end = 0.9) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle=45, vjust=1),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10))

###################################################

## Genus details
# export size 11x6
amp_heatmap(amp_subset_samples(
  amp.hel, !type %in% c(
    "fSW+YE_w_mon","SW+YE_w_mon",
    "fSW+YE_s_mon","SW+YE_s_mon") &
    !days %in% c("4","12","29","43","71")),
  tax_aggregate = "Genus",
  tax_show = c(
    "Psychrobacter","Colwellia",
    "Polaribacter","Moritella",
    "Flavobacterium","Maribacter",
    "Saccharospirillaceae_uc",
    "Marmoricola","KD4-96_uc",
    #"Pseudoalteromonas",
    #"Kordiimonadales_uc","Peredibacter",
    #"Croceibacter","Pedobacter",
    #"RS62 marine group","Pelagicoccus",
    #"Cryomorphaceae_uc","Roseibacillus",
    #"Reinekea",#"Yoonia-Loktanella",
    # "Zhongshania","Crocinitomix",
    "SAR11_Clade Ia","NS9 marine group_uc",
    "Magnetospiraceae_uc",
    "Marinimicrobia (SAR406 clade)_uc",
    "Methyloceanibacter","Nocardioides",
    "Oryzihumus","Gitt-GS-136_uc",
    "Aurantivirga","Pseudohongiella",
    "Gaiella","Methyloligellaceae_uc",
    "Alphaproteobacteria_uc",
    "Solirubrobacter","Iamia",
    "Gemmatimonadaceae_uc",
    "Flavobacteriaceae_uc",
    "67-14_uc","MB-A2-108_uc",
    "Vicinamibacterales_uc"),
  order_y_by = c(
    #"Peredibacter","Croceibacter",
    #"Cryomorphaceae_uc","Roseibacillus",
    #"Kordiimonadales_uc","Pedobacter",
    #"Pelagicoccus","RS62 marine group",
    #"Zhongshania","Crocinitomix",
    #"Pseudoalteromonas",
    "Pseudohongiella","Aurantivirga",
    "Maribacter","Flavobacteriaceae_uc",
    "Moritella","Colwellia",    
    "Saccharospirillaceae_uc",
    "Psychrobacter","Polaribacter",
    "Alphaproteobacteria_uc",
    "Magnetospiraceae_uc",
    "Marinimicrobia (SAR406 clade)_uc",
    "NS9 marine group_uc",
    "SAR11_Clade Ia","Marmoricola",
    "WCHB1-32","Methyloceanibacter",
    "Nocardioides","Flavobacterium",
    "Methyloligellaceae_uc","Iamia",
    "Oryzihumus","Gemmatimonadaceae_uc",
    "Vicinamibacterales_uc",
    "Solirubrobacter","MB-A2-108_uc",
    "Gaiella","Gitt-GS-136_uc",
    "KD4-96_uc","67-14_uc"),
  normalise = F,
  plot_values = F,
  facet_by = "type",
  group_by = "days") +
  scale_fill_fish(
   option="Exallias_brevis", #
   direction = -1,
   begin = 0.05, end = 0.99) +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(angle=0, vjust=1, hjust=0.5),
        legend.text = element_text(size=10),
        legend.title = element_text(size=10))

# Assign origin YE/SW
oriGen %>%
  filter(Genus %in% c("Psychrobacter","Colwellia",
"Polaribacter","Moritella",
"Flavobacterium","Maribacter",
"Saccharospirillaceae_uc",
"Marmoricola","KD4-96_uc",
"SAR11_Clade Ia","NS9 marine group_uc",
"Magnetospiraceae_uc",
"Marinimicrobia (SAR406 clade)_uc",
"Methyloceanibacter","Nocardioides",
"Oryzihumus","Gitt-GS-136_uc",
"Aurantivirga",
"Pseudohongiella","Gaiella",
"Methyloligellaceae_uc",
"Alphaproteobacteria_uc",
"Solirubrobacter","Iamia",
"Gemmatimonadaceae_uc",
"Flavobacteriaceae_uc",
"67-14_uc","MB-A2-108_uc",
"Vicinamibacterales_uc")) %>%
  mutate(Genus=factor(Genus, levels=c(
    #"Peredibacter","Croceibacter",
    #"Cryomorphaceae_uc","Roseibacillus",
    #"Kordiimonadales_uc","Pedobacter",
    #"Pelagicoccus","RS62 marine group",
    #"Zhongshania","Crocinitomix",
    #"Pseudoalteromonas",
    "Pseudohongiella","Aurantivirga",
    "Maribacter","Flavobacteriaceae_uc",
    "Moritella","Colwellia",    
    "Saccharospirillaceae_uc",
    "Psychrobacter","Polaribacter",
    "Alphaproteobacteria_uc",
    "Magnetospiraceae_uc",
    "Marinimicrobia (SAR406 clade)_uc",
    "NS9 marine group_uc",
    "SAR11_Clade Ia","Marmoricola",
    "WCHB1-32","Methyloceanibacter",
    "Nocardioides","Flavobacterium",
    "Methyloligellaceae_uc","Iamia",
    "Oryzihumus","Gemmatimonadaceae_uc",
    "Vicinamibacterales_uc",
    "Solirubrobacter","MB-A2-108_uc",
    "Gaiella","Gitt-GS-136_uc",
    "KD4-96_uc","67-14_uc"))) %>%
  mutate(data="data") %>%  # dummy for plotting
  mutate(type=str_to_upper(str_sub(type, 1, 1))) %>%
  ggplot(aes(x=data,y=Genus,fill=type)) +
  geom_tile(color="black") +
  geom_text(
    aes(label=type), 
    color="white", size=6) +
   scale_fill_manual(values=c(
     "S"="darkseagreen3",
     "P"="darkorchid4")) +
  theme_bw() +
  theme(axis.ticks = element_blank(),
        axis.text = element_blank(),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        legend.position = "none",
        axis.title = element_blank())
  
# Check original abundances of pot. responders
amp_heatmap(amp_subset_samples(
  amp.rel, type %in% c(
    "orgYE","orgSW","SW-postT")),
  tax_aggregate = "Genus",
  tax_show = c(
    "Psychrobacter","Colwellia",
    "Polaribacter","Moritella",
    "Aurantivirga","Maribacter",
    "Pseudohongiella",
    "Saccharospirillaceae_uc"),
  normalise = T,
  plot_values = T,
  plot_colorscale = "sqrt",
  round = 1,
  group_by = "type") +
  scale_fill_fish(
    option="Ostorhinchus_angustatus",
    begin = 0, end = 0.9)

###################################################

### DIVERSITY

AlphaDiv %>% filter(!type %in% c(
  "fSW+YE_w_mon","SW+YE_w_mon",
  "fSW+YE_s_mon","SW+YE_s_mon") &
    !days %in% c("4","12","29","43","71")) %>%
  reshape2::melt() %>%
  filter(variable %in% c("simpson","richness")) %>%
ggplot() +
  aes(x=type, y=value, fill=type) +
  geom_boxplot(outlier.shape = NA) + 
  scale_fill_manual(values=col) +
  facet_grid(variable~ ., scales="free") +
  #ylab("invSimpson index") +
  #ggtitle("AlphaDiversity") +
  theme_classic() +
  theme(
    #panel.grid.minor = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(
      angle = 45, vjust=1, hjust=1),
    legend.position = "none",
    axis.title.x = element_blank())

save.image("Permafrost.Rdata")
