############################################################################################
###  Correlations -- Fig. S9  ###
############################################################################################

# filter, only 1 ASV per Genus
asv <- ASV.rel %>%
  filter(rowSums(.>= 0.5) >= 3) 

tax <- TAX %>%
  rownames_to_column("asv") %>%
  distinct(Genus, .keep_all=T)

asv = as.data.frame(
  apply(asv, 2, function(x) sqrt(x / sum(x)))) %>%
  rownames_to_column("asv") %>%
  filter(asv %in% tax$asv) %>%
column_to_rownames("asv")

# ENV data
env <- ENV %>%
  dplyr::select(where(is.numeric)) %>%
  dplyr::select(-c("sampleVolume_mL"))

#######################

# Correlations
cor <- corr.test(
    t(asv), env, 
    use = "pairwise",
    method="spearman",
    adjust = "holm",
    alpha=.05, ci=T,
    minlength=5, normal=T) 

# extract, subset most significant; 
r <- cor$r 
r[abs(r) <= -0.35 | abs(r) <= 0.35 ] = NA 

# reformat; remove all-NA taxa
r <- as.data.frame(r) %>%
  rownames_to_column("asv") %>%
  filter_at(vars(2:7), any_vars(!is.na(.))) %>%
  #filter(asv!="sample volume (mL)" & is.na(`sample volume (mL)`)) %>%
  left_join(tax) %>% 
  distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, round, 2) 

# reformat p-values
p <- as.data.frame(cor$p) %>%
  rownames_to_column("asv") %>%
# filter(asv!="sample volume (mL)" & is.na(`sample volume (mL)`)) %>%
    distinct(asv, .keep_all=T) %>%
  reshape2::melt() %>%
  mutate_if(is.numeric, ~case_when(
    . < 0.05 & . > 0.01 ~ "*",
    . < 0.01 & . > 0.001 ~ "**",
    . < 0.001 ~ "***"),
    TRUE ~ NA) %>%
  dplyr::rename(pvalue=value)

##############################################

# Combine all info
# including ASV origin
cor <- right_join(r, p) %>%
  unite(sign, value, pvalue, sep="", remove=F) %>%
  mutate_at(c("sign"), ~gsub("NA", NA, .)) %>%
  mutate(across(value, as.numeric)) %>%
  filter(!Genus %in% c("Oleispira","Subgroup 7_uc") & 
           !Phylum %in% c("Firmicutes")) %>%  
  mutate(Class = case_when(
    !grepl("Alpha|Gamma|eroidia", Class)~"other", 
    TRUE ~ Class)) %>%
  mutate_at(c("Class"), ~gsub(
    "bacteria", "", .)) %>%
  mutate(Genus = str_replace(
    Genus,"_uc|uc|_clade|_cluster|_Clade_", "")) %>%
  mutate(id = paste(Genus,asv)) %>%
  mutate(id = str_replace(
    id,"_", " ")) %>%
  drop_na(sign) %>%
  left_join(oriASV) %>%
  group_by(asv) %>%
  mutate(origin = case_when(
    n_distinct(type) > 1~"unknown",
    TRUE ~ type)) %>%
  filter(origin!="unknown") %>%
  mutate(id=factor(id, levels=c(
    "Peredibacter asv136","Bacteriovoracaceae asv62",
    "AKYG1722 asv129","Marine Group II asv177",
    "Gemmatimonadaceae asv152","Oryzihumus asv30",
    "Marmoricola asv123","Gitt-GS-136 asv24",
    "Haliangium asv99","SAR11 Clade Ia asv28",
    "SAR11 clade asv160","SAR11 Clade II asv125",
    "Crocinitomix asv14","Croceibacter asv48",
    "Luteibaculum asv61","NS5 marine group asv78",
    "Methylophagaceae asv96","SAR92 clade asv118",
    "Rhodobacteraceae asv174","Halocynthiibacter asv54",
    "Sphingorhabdus asv66","Porticoccus asv94",
    "Pseudohongiella asv52","Marinoscillum asv85",
    "Saprospiraceae asv82","Flavobacteriaceae asv21",
    "Psychromonas asv72","Maribacter asv11",
    "Aquibacter asv25","Moritella asv6","Colwellia asv3",
    "Saccharospirillaceae asv26","Psychrobacter asv1"))) %>%
  mutate(variable=factor(variable, levels=c(
    "CO2_prod","NH4","NO2","PO4","SO4"))) %>%
  mutate(Class=factor(Class, levels=c(
    "Gammaproteo","Bacteroidia","Alphaproteo","other"))) 

# export size 7x5
ggplot(cor) + aes(
  x=id, 
  # x=paste(asv, Genus, sep = "_"),
  y=variable, 
  fill=value, label=pvalue) +
  geom_tile(color="gray84") +
  geom_text(color="white", size=3.8) +
  scale_fill_fish(
    alpha = NULL,
    begin = 0.1,
    end = 0.98,
    # direction = -1,
    na.value = "aliceblue",
    option = "Chaetodon_larvatus",  #Epinephelus_lanceolatus Chaetodon_ephippium
    limits = c(-0.65, 0.65),
    breaks = c(-0.6,-0.3,0,0.3,0.6)) +
  # scale_x_discrete(labels = cor$id) +
  facet_grid(
    Class~origin,  
    scales="free", 
    space="free") +
  coord_flip() +
  theme_classic() +
  theme(axis.ticks = element_blank(),
        axis.text.x = element_text(
          angle=45, hjust=1, vjust=1),
        panel.grid = element_blank(),
        panel.border = element_blank(),
        # strip.background = element_blank(),
        #legend.position = "none",
        axis.title = element_blank())

