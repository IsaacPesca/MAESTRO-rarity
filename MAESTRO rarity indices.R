
# Code Used to calculate the species rarity indices using the MAESTRO list of species 
# Isaac Trindade Santos
# MAESTRO workgroup
library(plotly)
library(htmlwidgets)
library(ggnewscale)
library(ape)
library(fishtree)
library(stringr)
library(picante)
library(pals)
library(ggpubr)
library(ggtree)
library(phyloregion)
library(caper)
library(bench)
library(phytools)

library(psych)
library(tidyverse)
library(gt)
library(glue)
library(paleotree)
library(dispRity)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree) # use the lines above if this package ggtree below is not installed on your library
library(castor)
library(rotl)
library(gridExtra)
library(rfishbase)
library(dplyr)
library(tidyr)
library(missForest)
library(RRphylo)


######## B - calculate evolutionary distinctiveness ########
## phylogenetic dis

subtree_MAESTRO666 <- readRDS("subtree_MAESTRO6.RDS")
# novo aqui
phy_tree  <- subtree_MAESTRO666

phy_tree$tip.label <- str_replace_all(phy_tree$tip.label, c("_" = " "))

# evolutionary distinctiveness
phy_tree.cm <- clade.matrix(phy_tree)
Evol_Di <- ed.calc(phy_tree.cm)
head(Evol_Di$spp)
most_distinct_species<-Evol_Di$spp[order(Evol_Di$spp$ED,decreasing=T),];head(most_distinct_species)
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
most_distinct_species$dist_scaled <- scale_values(most_distinct_species$ED)

MAESTRO_most_distinct_species<- most_distinct_species

#write.csv(most_distinct_species, "MAESTRO_most_distinct_speciesA.csv")

###### C - prepare trait matrix ######
# recorte comeca aqui

# load MAESTRO file
biomass_maestro_old <- read.csv("datras_medits_biomass_withtraits.csv", dec = ",", sep = ";")

biomass_maestro <- readRDS("MAESTRO_biomass_assemblage666.RDS")

new_list_of_species1 <- unique(biomass_maestro$genus_sp)

length(unique(biomass_maestro$genus_sp))
#  [1] 555
# MUDEI AQUI 27 -07 - 2025
biomass_maestro <- biomass_maestro_old

length(unique(biomass_maestro$genus_sp))
# 620 species
# add a column with a new ID to use at the ArcGIS
biomass_maestro$ID_m<-row.names(biomass_maestro)
species_maestro<- unique(biomass_maestro$genus_sp)
species_maestro<- as.character(species_maestro)

traits_MAESTRO<- biomass_maestro %>% dplyr::select(genus_sp, length.maturity, 
                                                   age.maturity, growth.coefficient, 
                                                   tl, fecundity, offspring.size, 
                                                   spawning.type, habitat, feeding.mode)
traits_MAESTRO<- unique(traits_MAESTRO)
## we have species with more than one trait value
## we need to select the maximum value available

traits_MAESTRO1<- traits_MAESTRO %>% 
  group_by(genus_sp) %>% 
  dplyr::summarise("length.maturity"  = max(length.maturity),
                   "age.maturity"  = max(age.maturity),
                   "growth.coefficient"  = max(growth.coefficient),
                   "tl"  = max(tl),
                   "fecundity"  = max(fecundity),
                   "offspring.size"  = max(offspring.size),
                   "spawning.type"  = unique(spawning.type),
                   "habitat"  = unique(habitat),
                   "feeding.mode"  = unique(feeding.mode))
## some species has multiple categorical trait values
## we need to select only one to proceed
traits_MAESTRO1[duplicated(traits_MAESTRO1$genus_sp),]


## Syngnathus_typhle
traits_MAESTRO1<- subset(traits_MAESTRO1[-575,]) 

## Hippocampus_guttulatus
traits_MAESTRO1<-traits_MAESTRO1[-253,]

## Dipturus_oxyrinchus
traits_MAESTRO1<-traits_MAESTRO1[-181,]

## Dipturus_nidarosiensis 
traits_MAESTRO1<-traits_MAESTRO1[-179,]

## now subset trait matrix using the new list of species from MAESTRO
## this new species comes from the SAC process
## and from the subset of 0.25 degree grid cells

traits_MAESTRO123 <- traits_MAESTRO1 %>% 
  dplyr::filter(genus_sp %in% new_list_of_species1)


write.csv(traits_MAESTRO123, "traits_MAESTRO1.csv")
row.names(traits_MAESTRO123)<- traits_MAESTRO123$genus_sp

# load all traits from fishbase with imputation
all_species_traits_fishbase_imputed <- read.csv("traits_FishBase_imputed.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(all_species_traits_fishbase_imputed)[colnames(all_species_traits_fishbase_imputed) == 'X'] <- 'genus_sp'
all_species_traits_fishbase_imputed1<- all_species_traits_fishbase_imputed %>% dplyr::select(genus_sp, SwimMode_fam, BodyShapeI_sp)

# corrigi aqui tb
traits_MAESTRO2 <- merge(traits_MAESTRO123, all_species_traits_fishbase_imputed1, by = "genus_sp", all.x = TRUE)

# lets combine with the fishbase table
# keep the MAESTRO traits
# use traits from FishBase only for the AquaMaps species that are missing from MAESTRO

MAESTRO_only1<- read.csv( "MAESTRO_only.txt")
MAESTRO_only1<- MAESTRO_only1$x
AQUAMAPS_only1<- read.csv("AQUAMAPS_only.txt", header=TRUE, stringsAsFactors=FALSE)
AQUAMAPS_only1<-AQUAMAPS_only1$x
species_in_common1<- read.csv("species_in_common.txt", header=TRUE, stringsAsFactors=FALSE )
species_in_common1<-species_in_common1$x

rownames(all_species_traits_fishbase_imputed)<- all_species_traits_fishbase_imputed$genus_sp
AQUAMAPS_only11<- all_species_traits_fishbase_imputed[rownames(all_species_traits_fishbase_imputed) %in% AQUAMAPS_only1, ]

# make the column names the same in both files
# traits MAESTRO
# traits AQUAMAPS
colnames(traits_MAESTRO2)

colnames_maestro1<- c("genus_sp", "lm", "tm", "k", 
                      "tl", "fecundity", "offspring.size", "repr.guild",
                      "habitat", "feeding.mode", "swim.mode", "body.shape"
)

colnames(traits_MAESTRO2)<- colnames_maestro1

colnames(AQUAMAPS_only11)

traits_AquaMaps_subset<- AQUAMAPS_only11 %>% 
  dplyr::select(genus_sp, Lm_ma,tm_ma,K_es,
                Troph_es,                 ReprGuild_fam,
                DemersPelag_sp,
                SwimMode_fam,BodyShapeI_sp
  )

colnames_AquaMaps1<- c("genus_sp", "lm", "tm", "k", 
                       "tl",  "repr.guild",
                       "habitat", 
                       "swim.mode", "body.shape")

colnames(traits_AquaMaps_subset)<- colnames_AquaMaps1
rownames(traits_AquaMaps_subset)<- NULL
traits_AquaMaps_subset$fecundity<- NA
traits_AquaMaps_subset$offspring.size<- NA
traits_AquaMaps_subset$feeding.mode<- NA

MAESTRO_AquaMaps_traits<- rbind(traits_MAESTRO2,traits_AquaMaps_subset)
rownames(MAESTRO_AquaMaps_traits)<- MAESTRO_AquaMaps_traits$genus_sp
# a few species has more than one entry in the matrix, keep only one
# Acipenser_sturio
MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits[-656,]
rownames(MAESTRO_AquaMaps_traits)<- MAESTRO_AquaMaps_traits$genus_sp

# Gobius_cobitis
# MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits[-1151,]

# make sure that the traits have the same names (MAESTRO trait names and FishBase trait names)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("nonguarders", "non-guarder", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("bearers", "bearer", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("guarders", "guarder", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("bearers_and_guarders", "mixed", MAESTRO_AquaMaps_traits$repr.guild)

# save rownames for later
row_names_MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits$genus_sp
row_names_MAESTRO_AquaMaps_traits<- as.matrix(row_names_MAESTRO_AquaMaps_traits)
row_names_MAESTRO_AquaMaps_traits<- data.frame(row_names_MAESTRO_AquaMaps_traits)

# corrigi aqui tb
MAESTRO_AquaMaps_traits123 <- MAESTRO_AquaMaps_traits %>% 
  dplyr::filter(genus_sp %in% new_list_of_species1)


# remove species names from matrix, the imputation does not require that
MAESTRO_AquaMaps_traits2<- MAESTRO_AquaMaps_traits123 %>% dplyr::select(-one_of('genus_sp' )) 
MAESTRO_AquaMaps_traits2 <- as.data.frame(unclass(MAESTRO_AquaMaps_traits2),stringsAsFactors=TRUE)

# make sure that columns are numeric or factors 
varClass(MAESTRO_AquaMaps_traits2) 



# this line below will run the Random Forest to fill gaps in trait matrix
MAESTRO_AquaMaps_traits3<- missForest(MAESTRO_AquaMaps_traits2, maxiter = 10, ntree = 100, variablewise = TRUE,
                                      decreasing = FALSE, verbose = TRUE,
                                      mtry = floor(sqrt(ncol(MAESTRO_AquaMaps_traits2))), replace = TRUE,
                                      classwt = NULL, cutoff = NULL, strata = NULL,
                                      sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                                      xtrue = NA, parallelize = c('no', 'variables', 'forests'))


traits_MAESTRO666_error_estimates<- MAESTRO_AquaMaps_traits3$OOBerror
MAESTRO_AquaMaps_traits_imputed<- MAESTRO_AquaMaps_traits3$ximp

# corrigi aqui tb
MAESTRO_AquaMaps_traits_imputed$genus_sp<- MAESTRO_AquaMaps_traits123$genus_sp

#traits_MAESTRO666<- traits_MAESTRO3$ximp
row.names(MAESTRO_AquaMaps_traits_imputed)<- MAESTRO_AquaMaps_traits_imputed$genus_sp
write.csv(MAESTRO_AquaMaps_traits_imputed, "MAESTRO_AquaMaps_traits_imputed.csv")

###### D - calculate functional distinctiveness ######
library(funrar)

# novo aqui
# importante 
# subset rownames with character vecotr containing the AquaMaps species names
MAESTRO_most_distinct_species$species<- gsub(" ", "_", MAESTRO_most_distinct_species$species)

rownames(MAESTRO_most_distinct_species) <- MAESTRO_most_distinct_species$species
rownames(MAESTRO_AquaMaps_traits_imputed) <- MAESTRO_AquaMaps_traits_imputed$genus_sp

traits_MAESTRO666_subset <- MAESTRO_AquaMaps_traits_imputed[rownames(MAESTRO_AquaMaps_traits_imputed) %in% MAESTRO_most_distinct_species$species, ] 
traits_MAESTRO666_subset<- traits_MAESTRO666_subset %>% dplyr::select(-one_of('genus_sp' )) 
# novo aqui
traitT<-compute_dist_matrix(traits_MAESTRO666_subset, metric="gower", center=FALSE, scale=FALSE)
# novo aqui
saveRDS(traitT, "traitT_distance_matrix_MAESTRO.RDS")
functional_distinctiveness<-distinctiveness_global(traitT, di_name="global_di")
colnames(functional_distinctiveness)[colnames(functional_distinctiveness) == 'species'] <- 'genus_sp'
MAESTRO_species_indices<- most_distinct_species
colnames(MAESTRO_species_indices)[colnames(MAESTRO_species_indices) == 'dist_scaled'] <- 'evol_dist_scaled'
colnames(MAESTRO_species_indices)[colnames(MAESTRO_species_indices) == 'species'] <- 'genus_sp'

MAESTRO_species_indices$genus_sp <- chartr(" ", "_", MAESTRO_species_indices$genus_sp )
MAESTRO_species_indices <- merge(MAESTRO_species_indices, functional_distinctiveness, by = "genus_sp", all.x = TRUE)

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_species_indices$funct_dist_scaled <- scale_values(MAESTRO_species_indices$global_di)

# recorte termina aqui
###### E - prepare biomass matrix #########

# load biomass novo aqui
MAESTRO_biomass_assemblage<- read.csv("MAESTRO_biomass_assemblage.csv", h=T)
# # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # 
# # # # # # # # # # # # # # # # # # 

length(unique(MAESTRO_biomass_assemblage$genus_sp))

biomass_matrix<- MAESTRO_biomass_assemblage %>% dplyr::select(genus_sp, Year, biomass)

biomass_matrix1<- biomass_matrix %>% 
  group_by(genus_sp, Year) %>% 
  dplyr::summarise( "Year" = unique(Year),
                    "biomass"  = sum(biomass))

biomass_matrix2<- biomass_matrix1 %>% pivot_wider(names_from = genus_sp, values_from = biomass, values_fill = 0)
biomass_matrix2<-as.matrix(biomass_matrix2)

biomass_matrix2<-as.data.frame(biomass_matrix2)
row.names(biomass_matrix2)<- biomass_matrix2$Year
biomass_matrix3<- biomass_matrix2[,-1]
biomass_matrix3<-as.matrix(biomass_matrix3)
biomass_matrix4<-make_relative(biomass_matrix3)
biomass_matrix5<-biomass_matrix4[,MAESTRO_species_indices$genus_sp]

######### F - calculate taxonomic scarcity ######
# novo aqui
#most_distinct_species$species <- chartr(" ", "_", most_distinct_species$species )
row_names1<- MAESTRO_species_indices$genus_sp
# MAESTRO_AquaMaps_traits_imputed
# novo aqui
traits_MAESTRO6666 <- MAESTRO_AquaMaps_traits_imputed[rownames(MAESTRO_AquaMaps_traits_imputed) %in% row_names1, ]
traits_MAESTRO6666<- traits_MAESTRO6666 %>% dplyr::select(-one_of('genus_sp' )) 
traits_distance_matrix2<-compute_dist_matrix(traits_MAESTRO6666, metric="gower", center=FALSE, scale=FALSE)

funrar_MAESTRO<- funrar(biomass_matrix5, traits_distance_matrix2, rel_abund = TRUE)
# this is scarcity across years
scarcity_MAESTRO<- funrar_MAESTRO$Si

library(tidyr)

scarcity_MAESTRO1<- colMeans(scarcity_MAESTRO, na.rm=TRUE)
scarcity_MAESTRO1<- as.matrix(scarcity_MAESTRO1)
scarcity_MAESTRO1<- as.data.frame(scarcity_MAESTRO1)
colnames(scarcity_MAESTRO1)[colnames(scarcity_MAESTRO1) == 'V1'] <- 'scarcity_average'

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_species_indices$scarcity_average <- scarcity_MAESTRO1$scarcity_average
MAESTRO_species_indices$scarcity_scaled <- scale_values(MAESTRO_species_indices$scarcity_average)

biomass_matrix_log<- MAESTRO_biomass_assemblage %>% 
  group_by(genus_sp) %>% 
  dplyr::summarise( "biomass"  = sum(biomass))
biomass_matrix_log$biomass_log<- log10(biomass_matrix_log$biomass)
biomass_matrix_log$biomass_log_scaled <- scale_values(biomass_matrix_log$biomass_log)
MAESTRO_species_indices <- merge(MAESTRO_species_indices, biomass_matrix_log, by = "genus_sp", all.x = TRUE)
# this line below will add a new column that represents scarcity
# scarcity is the opposite of abundant, then we can use 1 - ABUNDANCE/BIOMASS
# that is why I included 1 - below
MAESTRO_species_indices$scarcity_logbio_scaled<- 1 - MAESTRO_species_indices$biomass_log_scaled

write.csv(MAESTRO_species_indices,"MAESTRO_species_indices1.csv")



###### E - prepare presence absence matrix #########

MAESTRO_pres_absc_long11<- MAESTRO_pres_absc_long1 %>% dplyr::select(genus_sp, OID_First)
MAESTRO_pres_absc_long11<- unique(MAESTRO_pres_absc_long11)
MAESTRO_pres_absc_long11$OID_First<- as.factor(MAESTRO_pres_absc_long11$OID_First)

library(reshape2)

# transform the long format to a wide format
# this is the presence absence matrix from AquaMap
# correct column name
# cname<- c("OID_First", "genus_sp")
# colnames(MAESTRO_pres_absc_long11)<- cname

# AquaMaps_species_list<- AquaMaps_species_indices$genus_sp
# keep only species with phylogeny information
# create object with list of species without phylogeny information
# AQUAMAPS_only <- setdiff(AquaMaps_pres_absc_long1$genus_sp, AquaMaps_species_indices$genus_sp )
# use object to remove from the list
# AquaMaps_pres_absc_long16 <- AquaMaps_pres_absc_long1[!AquaMaps_pres_absc_long1$genus_sp %in% AQUAMAPS_only, ]

# importante
# traits_AquaMaps666_subset <- AquaMaps_pres_absc_long1[rownames(AquaMaps_pres_absc_long1) %in% AquaMaps_most_distinct_species$species, ] 

MAESTRO_pres_absc_wide<- dcast(MAESTRO_pres_absc_long11, OID_First~genus_sp, length)

# make sure that it is a matrix
typeof(MAESTRO_pres_absc_wide)

rownames(MAESTRO_pres_absc_wide)<-MAESTRO_pres_absc_wide$OID_First

MAESTRO_pres_absc_wide1<-MAESTRO_pres_absc_wide[,-1]

MAESTRO_pres_absc_wide2<- data.frame(MAESTRO_pres_absc_wide1)
MAESTRO_pres_absc_wide2<- as.matrix(MAESTRO_pres_absc_wide2)
# importante


MAESTRO_sp_restrictedness<- restrictedness(MAESTRO_pres_absc_wide2, relative = FALSE)

rest_colnames<- c("genus_sp", "Ri")
colnames(MAESTRO_sp_restrictedness)<- rest_colnames
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_sp_restrictedness$restrictedness_scaled <- scale_values(as.numeric(MAESTRO_sp_restrictedness$Ri))

MAESTRO_species_indices123<- merge(MAESTRO_species_indices,MAESTRO_sp_restrictedness, by = "genus_sp", all.x = TRUE)
str(MAESTRO_species_indices123)

MAESTRO_species_indices1234<- MAESTRO_species_indices123 %>% dplyr::select(genus_sp,evol_dist_scaled, funct_dist_scaled, scarcity_logbio_scaled, restrictedness_scaled )
row.names(MAESTRO_species_indices1234)<- MAESTRO_species_indices1234$genus_sp
MAESTRO_species_indices1234<- MAESTRO_species_indices1234[,-1]

# now that we have all indices needed for MAESTRO we can check the correlations between the indices
# low correlations indicates complementarity between the indices

#create pairs plot
pairs.panels(MAESTRO_species_indices1234, main= "MAESTRO - Biodiversity Facet`s Correlations", lm=TRUE)

write.csv(MAESTRO_species_indices123, "MAESTRO_species_indices123.csv")

