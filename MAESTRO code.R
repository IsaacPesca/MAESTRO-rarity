
#title: "MAESTRO"
#author: "Isaac Trindade Santos"
#date: "2024-06-01"

### A - build a tree for MAESTRO species
### B - calculate evolutionary distinctiveness
### C - prepare trait matrix
### D - calculate functional distinctiveness
### E - prepare biomass matrix
### F - calculate taxonomic scarcity 
### G - build a 3D plot with taxonomic, functional and phylogenetic indices
### H - bring the MAESTRO species`s vulnerability indices
### I - use the vulnerability indices to color code the species inside the 3D plot

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

######## Table from FishBase with all species (rows) and traits (columns) available 
fishbase_all_sp6 <- read.csv("fishbase_all_sp555.csv", header=T, na.strings=c("","NA"))
fishbase_all_sp7 <- fishbase_all_sp6[-c(1:2), ]
fishbase_all_sp7 <- fishbase_all_sp7[,-c(1)]
fishbase_all_sp7 <- data.frame(fishbase_all_sp7)
species_families <- frame()
species_families <- fishbase_all_sp7 %>% select(Sp, Family, Class, Order) 
write.csv(species_families , "species_families.csv")

######## A - build a tree for MAESTRO species  #########


#___________________FIRST PART___________________________________
####### load JAWLESS fish phylogeny

jawless_fish1<-read.tree("jawlessFishTree.tre")
class(jawless_fish1)<- "phylo"

jawless_fish1<-fix.poly(jawless_fish1,type="resolve")
#jawless_fish1$Nnode
#[1] 16
jawless_fish1$node.label<- c(1:16)

joutgroup<-c("Myxine_glutinosa", "Petromyzon_marinus","Geotria_australis", "Lampetra_fluviatilis", "Ichthyomyzon_greeleyi", "Ichthyomyzon_unicuspis" )

rooted_tree = root(jawless_fish1, outgroup = joutgroup, resolve.root = TRUE)
is.rooted(jawless_fish1)

class(jawless_fish1)

#ggtree(jawless_fish1, layout='circular') + theme_tree2() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

jawless_fish1<- multi2di(jawless_fish1)
max(branching.times(jawless_fish1))
#plot(jawless_fish1)
is.ultrametric(jawless_fish1)


#___________________SECOND PART___________________________________
# load the shark consensus tree
shark_consensus_MAESTRO<- readRDS("shark_consensus_updated.RDS")
shark_consensus_MAESTRO$Nnode
#[1] 1153
shark_consensus_MAESTRO$node.label<- c(1:1153)

shark_consensus_MAESTRO<-fix.poly(shark_consensus_MAESTRO,type="resolve")

shark_consensus_MAESTRO<- multi2di(shark_consensus_MAESTRO)
is.rooted(shark_consensus_MAESTRO)

#ggtree(shark_consensus_MAESTRO, layout='circular') + theme_tree2() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

# is it binary?
is.binary(shark_consensus_MAESTRO)
# [1] TRUE
# is it ultrametric?
is.ultrametric(shark_consensus_MAESTRO)
# [1] TRUE
max(branching.times(shark_consensus_MAESTRO))
# [1] 382.128
list_of_sharks<- shark_consensus_MAESTRO$tip.label
list_of_sharks<-as.data.frame(list_of_sharks)


#___________________THIRD PART___________________________________
# get the bony fish trees
fish_tree_file_MAESTRO32 <- readRDS("treefish32_ALL_updated_corrected.RDS") 
sampled_numbers1<-sample(100, size = 1)  
species32<- fish_tree_file_MAESTRO32[[sampled_numbers1]]

species32$node.label<- c(1:31080)
class(species32) <- "phylo"

is.rooted(species32)
#[1] TRUE
max(branching.times(species32))
# [1] 368.027

#___________________FOURTH PART___________________________________
# bind both trees here
sharks_bonyfish_MAESTRO<-bind.tree(shark_consensus_MAESTRO, species32, where = "root", position = 0, interactive = FALSE)

#ggtree(sharks_bonyfish_MAESTRO, layout='circular') + theme_tree2()


#sharks_bonyfish_MAESTRO <- phytools::force.ultrametric(sharks_bonyfish_MAESTRO)
#is.ultrametric(sharks_bonyfish_MAESTRO)

#sharks_bonyfish_MAESTRO<- multi2di(sharks_bonyfish_MAESTRO)
#is.ultrametric(sharks_bonyfish_MAESTRO)

species<- sharks_bonyfish_MAESTRO$tip.label
length(species)
#[1] 32235
max(branching.times(sharks_bonyfish_MAESTRO))
#[1] 382.128
max(sharks_bonyfish_MAESTRO$edge.length)
#[1] 382.128

# is it binary?
is.binary(sharks_bonyfish_MAESTRO)
# [1] TRUE
length(sharks_bonyfish_MAESTRO$tip.label)
# [1] 32235


ape::write.tree(sharks_bonyfish_MAESTRO, file='sharks_bonyfish_MAESTRO_tree.txt')



## bind jawless fish now
jawlessfish_sharks_bonyfish_MAESTRO<-bind.tree(sharks_bonyfish_MAESTRO, jawless_fish1, where = "root", position = 0, interactive = FALSE)
jawlessfish_sharks_bonyfish_MAESTRO<- multi2di(jawlessfish_sharks_bonyfish_MAESTRO)
length(species)
#[1] 32235
#write.csv(species, "fish_shark_species.csv")
max(branching.times(jawlessfish_sharks_bonyfish_MAESTRO))

max(jawlessfish_sharks_bonyfish_MAESTRO$edge.length)
#[1] 382.128

# is it binary?
is.binary(jawlessfish_sharks_bonyfish_MAESTRO)
# [1] TRUE
length(jawlessfish_sharks_bonyfish_MAESTRO$tip.label)
# [1] 32235
ape::write.tree(jawlessfish_sharks_bonyfish_MAESTRO, file='mega_tree.txt')
saveRDS(jawlessfish_sharks_bonyfish_MAESTRO, "jawlessfish_sharks_bonyfish_MAESTRO.RDS")
#mega_tree<- ape::read.tree("mega_tree.txt")

#mega_tree_circular <-ggtree(mega_tree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
#  theme_tree2() +
#  ggtitle(label = "Bony, Cartilaginous and Jawless Fish Tree",subtitle = " 32235 Species")




#___________________FIFTH PART___________________________________
####
## load species list from MAESTRO 
####

biomass_maestro <- read.csv("datras_medits_biomass_withtraits.csv", dec = ",", sep = ";")
length(unique(biomass_maestro$genus_sp))
# 620 species
# add a column with a new ID to use at the ArcGIS
biomass_maestro$ID<-row.names(biomass_maestro)
species_maestro<- unique(biomass_maestro$genus_sp)
species_maestro<- as.character(species_maestro)

# remove species from tree
species_maestro <-species_maestro[! species_maestro %in% c("Acipenser_sturio", "Cetorhinus_maximus")]



## subset our MAESTRO tree here

#subtree <- castor::get_subtree_with_tips(jawlessfish_sharks_bonyfish_MAESTRO, only_tips=species_maestro)$subtree
subtree_MAESTRO666<- readRDS("subtree_MAESTRO6.RDS")
#ggtree(subtree, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Jawless, Bony and Cartilaginous fish")

is.binary(subtree_MAESTRO666)
#[1] TRUE
subtree <- phytools::force.ultrametric(subtree)
is.ultrametric(subtree)
#[1] TRUE

#___________________PLOT TREES TOGETHER_______________________________
# shark tree
# length(shark_consensus_MAESTRO$tip.label)
#[1] 1154
shr2<-ggtree(shark_consensus_MAESTRO, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label="Sharks - consensus tree", subtitle = "1154 Species")

####
# jawless fish
# length(jawless_fish1$tip.label)
#[1] 1154
jaw2<-ggtree(jawless_fish1, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label="Jawless Fish", subtitle = "17 Species")
####
# bony fish tree
# length(species32$tip.label)
#[1] 31081
bf2<-ggtree(species32, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label = "Bony fish - fish tree of life", subtitle = "31081 Species")
####
# binded tree - bony + shark
# length(jawlessfish_sharks_bonyfish_MAESTRO$tip.label)
# [1] 32235
MAESTRO1<- ggtree(jawlessfish_sharks_bonyfish_MAESTRO, color="grey", size=0.1,) + 
  theme_tree2() + 
  ggtitle(label = "Jawless, Bony and Shark Fish Tree", subtitle = "32252 Species")
####
#MAESTRO species tree
#length(subtree$tip.label) 
# [1] 609
#MAESTRO all species tree
#length(unique(biomass_maestro$genus_sp))
#  [1] 620


MAESTRO2 <-ggtree(subtree, color="grey", size=0.001, ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "602 Species")


MAESTRO_circular <-ggtree(subtree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "602 Species")

####
grid.arrange(jaw2,bf2,shr2,MAESTRO1,MAESTRO2, ncol=5)
####
## save final tree for MAESTRO project
subtree$node.label<-c(1:604)
is.binary(subtree)
is.ultrametric(subtree)
subtree$Nnode
saveRDS(subtree,"MAESTRO_subtree_604_species_Jun2024.RDS")
# novo
length(subtree_MAESTRO666$tip.label)
# novo
MAESTRO_circular <-ggtree(subtree_MAESTRO666, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "589 Species")

MAESTRO_circular

MAESTRO_tree_species<- subtree_MAESTRO666$tip.label
ALL_species<- unique(biomass_maestro$genus_sp)

species_not_in_the_tree<- MAESTRO_tree_species[!ALL_species %in% MAESTRO_tree_species ]
species_not_in_the_tree
# list of species from MAESTRO that are missing in the tree
# novo
# [1] "Bathyraja_richardsoni"     "Scymnodon_ringens"         "Muraena_helena"            "Ariosoma_balearicum"      
# [5] "Nansenia_oblita"           "Manducus_maderensis"       "Synodus_saurus"            "Lestidiops_sphyrenoides"  
# [9] "Lampanyctus_intricarius"   "Notoscopelus_bolini"       "Symbolophorus_veranyi"     "Gaidropsarus_vulgaris"    
# [13] "Gaidropsarus_argentatus"   "Ciliata_mustela"           "Molva_dypterygia"          "Echiodon_drummondii"      
# [17] "Symphodus_bailloni"        "Ephippion_guttifer"        "Eutrigla_gurnardus"        "Leptoclinus_maculatus"    
# [21] "Lumpenus_lampretaeformis"  "Anarhichas_minor"          "Icelus_bicornis"           "Myoxocephalus_scorpioides"
# [25] "Pomatomus_saltatrix"       "Benthodesmus_simonyi"      "Lesueurigobius_suerii"     "Remora_brachyptera"       
# [29] "Citharus_linguatula"       "Scophthalmus_maximus"

######## B - calculate evolutionary distinctiveness ########
## phylogenetic dis

subtree_MAESTRO6 <- readRDS("subtree_MAESTRO6.RDS")
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

write.csv(most_distinct_species, "MAESTRO_most_distinct_speciesA.csv")

###### C - prepare trait matrix ######
# recorte comeca aqui

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

write.csv(traits_MAESTRO1, "traits_MAESTRO1.csv")
row.names(traits_MAESTRO1)<- traits_MAESTRO1$genus_sp

# load all traits from fishbase with imputation
all_species_traits_fishbase_imputed <- read.csv("traits_FishBase_imputed.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(all_species_traits_fishbase_imputed)[colnames(all_species_traits_fishbase_imputed) == 'X'] <- 'genus_sp'
all_species_traits_fishbase_imputed1<- all_species_traits_fishbase_imputed %>% dplyr::select(genus_sp, SwimMode_fam, BodyShapeI_sp)


traits_MAESTRO2 <- merge(traits_MAESTRO1, all_species_traits_fishbase_imputed1, by = "genus_sp", all.x = TRUE)

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

# a few species has more than one entry in the matrix, keep only one
# Gobius_cobitis
MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits[-1151,]
# Acipenser_sturio
MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits[-656,]

MAESTRO_AquaMaps_traits$repr.guild<- gsub("nonguarders", "non-guarder", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("bearers", "bearer", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("guarders", "guarder", MAESTRO_AquaMaps_traits$repr.guild)
MAESTRO_AquaMaps_traits$repr.guild<- gsub("bearers_and_guarders", "mixed", MAESTRO_AquaMaps_traits$repr.guild)

# save rownames for later
row_names_MAESTRO_AquaMaps_traits<- MAESTRO_AquaMaps_traits$genus_sp
row_names_MAESTRO_AquaMaps_traits<- as.matrix(row_names_MAESTRO_AquaMaps_traits)
row_names_MAESTRO_AquaMaps_traits<- data.frame(row_names_MAESTRO_AquaMaps_traits)

# remove species names from matrix, the imputation does not require that
MAESTRO_AquaMaps_traits2<- MAESTRO_AquaMaps_traits %>% dplyr::select(-one_of('genus_sp' )) 
MAESTRO_AquaMaps_traits2 <- as.data.frame(unclass(MAESTRO_AquaMaps_traits2),stringsAsFactors=TRUE)

# make sure that columns are numeric or factors 
varClass(MAESTRO_AquaMaps_traits2) 

MAESTRO_AquaMaps_traits3<- missForest(MAESTRO_AquaMaps_traits2, maxiter = 10, ntree = 100, variablewise = TRUE,
                                      decreasing = FALSE, verbose = TRUE,
                                      mtry = floor(sqrt(ncol(MAESTRO_AquaMaps_traits2))), replace = TRUE,
                                      classwt = NULL, cutoff = NULL, strata = NULL,
                                      sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                                      xtrue = NA, parallelize = c('no', 'variables', 'forests'))


traits_MAESTRO666_error_estimates<- MAESTRO_AquaMaps_traits3$OOBerror
MAESTRO_AquaMaps_traits_imputed<- MAESTRO_AquaMaps_traits3$ximp
MAESTRO_AquaMaps_traits_imputed$genus_sp<- row_names_MAESTRO_AquaMaps_traits$row_names_MAESTRO_AquaMaps_traits

#traits_MAESTRO666<- traits_MAESTRO3$ximp
row.names(MAESTRO_AquaMaps_traits_imputed)<- MAESTRO_AquaMaps_traits_imputed$genus_sp
write.csv(MAESTRO_AquaMaps_traits_imputed, "MAESTRO_AquaMaps_traits_imputed.csv")

###### D - calculate functional distinctiveness ######
library(funrar)

# novo aqui
# importante 
# subset rownames with character vecotr containing the AquaMaps species names
MAESTRO_most_distinct_species$species<- gsub(" ", "_", MAESTRO_most_distinct_species$species)
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
MAESTRO_species_indices$scarcity_logbio_scaled<- 1 - MAESTRO_species_indices$biomass_log_scaled

write.csv(MAESTRO_species_indices,"MAESTRO_species_indices1.csv")

MAESTRO_species_indices1<- MAESTRO_species_indices %>% dplyr::select(genus_sp,evol_dist_scaled, funct_dist_scaled, scarcity_logbio_scaled )
row.names(MAESTRO_species_indices1)<- MAESTRO_species_indices1$genus_sp
MAESTRO_species_indices1<- MAESTRO_species_indices1[,-1]
#create pairs plot
pairs.panels(MAESTRO_species_indices1, main= "MAESTRO - Biodiversity Facet`s Correlations", lm=TRUE)
library(ppcor)
pcor(MAESTRO_species_indices1)
library(dispRity)
tree.age(subtree)



## here starts plots

###### G - build a 3D plot with taxonomic, functional and phylogenetic indices ######
#saveRDS(subtree, "subtree_MAESTRO.rds")
#subtree_MAESTRO<- readRDS("subtree_MAESTRO.rds")

phy_tree <- subtree_MAESTRO #### mudei aqui
phy_tree$tip.label <- str_replace_all(phy_tree$tip.label, c(" " = "_"))

# plot phylogenetic tree
# Define colors
cols <- rev(jet(10))

# add several variables                

p <- ggtree(phy_tree, layout='circular', color="grey50", size=0.3) +
  annotate('text', x=0, y=40, label='', family='mono', size=16)# +
#geom_tiplab(cex=3,offset = 32,size = 3,geom = "text",linetype = "blank",align = TRUE)
#p2 <- gheatmap(p, EvolDi, width=0.2, hjust='left', colnames_angle=-10, font.size=1.5)  +
#      new_scale_fill()+
#      scale_fill_manual(values=c("#E41A1C","#377EB8","#FC8D59")) + theme_tree()

p2 <- gheatmap(p, MAESTRO_species_indices1[,c(1:3)], width=0.2, hjust='left', colnames_angle=-10, font.size=3)  +
  # new_scale_fill() +
  #scale_fill_gradient2 # if we have negative values in the variable, use scale_fill_gradient2 
  scale_fill_gradient2(low="blue", high="red", guide="colorbar") 
theme_tree()     

open_tree(p2, 20) %>% rotate_tree(70)                


library(rgl)
plot3d(MAESTRO_species_indices1$scarcity_logbio_scaled,MAESTRO_species_indices1$funct_dist_scaled,MAESTRO_species_indices1$evol_dist_scaled, xlab="Scarcity", ylab="Functional Distinc.", zlab="Evolutionary Distinct.", pch=5 )
#library(ggplot2)
#library(plotly)
#ggplot(MAESTRO_species_indices1, aes(x=scarcity_scaled, y=funct_dist_scaled, z=evol_dist_scaled)) +   #, color=Species
#  theme_void() 
# rgl::axes_3D() +
#  stat_3D()



library(plotly)

plot_3d <- MAESTRO_species_indices1 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                y=~funct_dist_scaled,
                                                z=~evol_dist_scaled, alpha = 0.5)


plot_3d




q.evol<- quantile(MAESTRO_species_indices1$evol_dist_scaled, c(0.25, 0.5, 0.75), type=1)
#25%        50%        75% 
#0.08712192 0.12865583 0.18725516 

q.funct<-quantile(MAESTRO_species_indices1$funct_dist_scaled, c(0.25, 0.5, 0.75), type=1)
#25%       50%       75% 
#0.1051994 0.1828711 0.2750250 

q.scarcity<-quantile(MAESTRO_species_indices1$scarcity_logbio_scaled, c(0.25, 0.5, 0.75), type=1)
#25%       50%       75% 
#0.9782451 0.9986295 0.9999192
#MAESTRO_species_indices2$abundant_scaled<- 1- MAESTRO_species_indices2$scarcity_scaled
MAESTRO_species_indices1$abundant_logbio_scaled<- 1- MAESTRO_species_indices1$scarcity_logbio_scaled


q.abundant<-quantile(MAESTRO_species_indices2$abundant_scaled, c(0.25, 0.5, 0.75), type=1)
#25%          50%          75% 
#0.0000778267 0.0013383711 0.0213682585 


MAESTRO_species_indices2<- MAESTRO_species_indices1 %>% 
  mutate(evol_dist_quartile = case_when(evol_dist_scaled <= q.evol[[1]] ~ "Quartile1",
                                      evol_dist_scaled >= q.evol[[1]]  & evol_dist_scaled <= q.evol[[2]] ~ "Quartile2",
                                      evol_dist_scaled >= q.evol[[2]]  & evol_dist_scaled <= q.evol[[3]] ~ "Quartile3",
                                      evol_dist_scaled >= q.evol[[3]] ~ "Quartile4")) 

MAESTRO_species_indices2<- MAESTRO_species_indices2 %>% 
  mutate(funct_dist_quartile = case_when(funct_dist_scaled <= q.funct[[1]] ~ "Quartile1",
                                         funct_dist_scaled >= q.funct[[1]]  & funct_dist_scaled <= q.funct[[2]] ~ "Quartile2",
                                         funct_dist_scaled >= q.funct[[2]]  & funct_dist_scaled <= q.funct[[3]] ~ "Quartile3",
                                         funct_dist_scaled >= q.funct[[3]] ~ "Quartile4"))

MAESTRO_species_indices2<- MAESTRO_species_indices2 %>% 
  mutate(scarcity_scaled_quartile = case_when(scarcity_logbio_scaled <= q.scarcity[[1]] ~ "Quartile1",
                                              scarcity_logbio_scaled >= q.scarcity[[1]]  & scarcity_logbio_scaled <= q.scarcity[[2]] ~ "Quartile2",
                                              scarcity_logbio_scaled >= q.scarcity[[2]]  & scarcity_logbio_scaled <= q.scarcity[[3]] ~ "Quartile3",
                                              scarcity_logbio_scaled >= q.scarcity[[3]] ~ "Quartile4")) 

MAESTRO_species_indices3<- MAESTRO_species_indices2 %>% 
  mutate(RARE_fun_evol = case_when(funct_dist_scaled <= q.funct[[1]] &  evol_dist_scaled <= q.evol[[1]]~ "Common2D",
                                   funct_dist_scaled >= q.funct[[3]] & evol_dist_scaled >= q.evol[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_fun_evol= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_tax_evol = case_when(scarcity_logbio_scaled <= q.scarcity[[1]] &  evol_dist_scaled <= q.evol[[1]]~ "Common2D",
                                   scarcity_logbio_scaled >= q.scarcity[[3]] & evol_dist_scaled >= q.evol[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_tax_evol= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_tax_funct = case_when(funct_dist_scaled <= q.funct[[1]]  &  scarcity_logbio_scaled <= q.scarcity[[1]]~ "Common2D",
                                    funct_dist_scaled >= q.funct[[3]] & scarcity_logbio_scaled >= q.scarcity[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_tax_funct= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_3D = case_when(funct_dist_scaled <= q.funct[[1]]  &  scarcity_logbio_scaled <= q.scarcity[[1]]&  evol_dist_scaled <= q.evol[[1]]~ "Common3D",
                                    funct_dist_scaled >= q.funct[[3]] & scarcity_logbio_scaled >= q.scarcity[[3]]&  evol_dist_scaled >= q.evol[[3]]~ "Rare3D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_3D= "Mid3D"))

ggplot(MAESTRO_species_indices3, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=RARE_fun_evol)) + theme_bw() +
  geom_point(aes(size=scarcity_logbio_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")

#ggplot(MAESTRO_species_indices3, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=RARE_3D)) + theme_bw() +
#  geom_point(aes(size=scarcity_logbio_scaled), alpha= 0.5)+
#  #  geom_smooth(method=lm)+
#  geom_vline(xintercept=q.funct[[3]])+
#  geom_hline(yintercept=q.evol[[3]])+ 
#  xlab("Functional Distinctiveness")+
#  ylab("Evolutinary Distinctiveness")


ggplot(MAESTRO_species_indices3, aes(x=scarcity_logbio_scaled, y=evol_dist_scaled, color=RARE_tax_evol)) + theme_bw() +
   geom_point(aes(size=funct_dist_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.scarcity[[3]], colour = "blue")+
  geom_vline(xintercept=q.scarcity[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Taxonomic Rarity")+
  ylab("Evolutinary Distinctiveness")

ggplot(MAESTRO_species_indices3, aes(x=scarcity_logbio_scaled, y=funct_dist_scaled, color=RARE_tax_funct)) +  theme_bw() +
  geom_point(aes(size=evol_dist_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.scarcity[[3]], colour = "blue")+
  geom_vline(xintercept=q.scarcity[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.funct[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.funct[[1]],linetype = "dashed", colour = "red")+ 
  xlab("Taxonomic Rarity")+
  ylab("Functional Distinctiveness")

###### H - bring the MAESTRO species`s vulnerability indices ######

# vdata - spatially-explicit species risk estimates

# number of species with VULNERABILITY INIDICES
vdata<- vdata
length(unique(vdata$latname))
# [1] 366

vdata_MAESTRO<- vdata %>% dplyr::select(latname, S.TSMr, 
                                         S.rlstatus, S.vind, 
                                         S.HII, AC.tvar, AC.hrange,
                                        AC.lmax, AC.hfrag, E.toe,E.vel,
                                        E.plost,E.nrchng)
vdata_MAESTRO$latname <- sub(" ", "_", vdata_MAESTRO$latname)
vdata_MAESTRO<- unique(vdata_MAESTRO)


# vspec - spatially-non-explicit species risk estimates
# number of species with VULNERABILITY INIDICES
length(unique(vspec$latname))
# [1] 366
vspec<- vspec

vspec_MAESTRO<- vspec %>% dplyr::select(latname,rcp, S.TSMr, 
                                        S.rlstatus, S.vind, 
                                        S.HII, AC.tvar, AC.hrange,
                                        AC.lmax, AC.hfrag, E.toe,E.vel,
                                        E.plost,E.nrchng,
                                        sens, adcap,
                                        expo, vuln,
                                        scat,acat,vcat,ecat
                                        )

vspec_MAESTRO$latname <- sub(" ", "_", vspec_MAESTRO$latname)
colnames(vspec_MAESTRO)[colnames(vspec_MAESTRO) == 'latname'] <- 'genus_sp'
vspec_MAESTRO26<- subset(vspec_MAESTRO, rcp== "2.6")
vspec_MAESTRO85<- subset(vspec_MAESTRO, rcp== "8.5")


MAESTRO_species_indices3$genus_sp <- rownames(MAESTRO_species_indices3)

MAESTRO_species_indices4<- MAESTRO_species_indices3 %>% dplyr::select(genus_sp,
  evol_dist_scaled,scarcity_logbio_scaled, funct_dist_scaled
                                              )
MAESTRO_species_indices26 <- merge(MAESTRO_species_indices4,vspec_MAESTRO26, by = "genus_sp" )
MAESTRO_species_indices85 <- merge(MAESTRO_species_indices4,vspec_MAESTRO85, by = "genus_sp" )

MAESTRO_species_indices2266 <- merge(MAESTRO_species_indices4,vspec_MAESTRO26, by = "genus_sp", all.x = TRUE )
MAESTRO_species_indices8855 <- merge(MAESTRO_species_indices4,vspec_MAESTRO85, by = "genus_sp" , all.x = TRUE)
MAESTRO_species_indices2266 <- MAESTRO_species_indices2266 %>% replace(is.na(.), 0)
MAESTRO_species_indices2266$ecat<- gsub("0", "No_Info", MAESTRO_species_indices2266$ecat)
MAESTRO_species_indices2266$vcat<- gsub("0", "No_Info", MAESTRO_species_indices2266$vcat)
MAESTRO_species_indices2266$acat<- gsub("0", "No_Info", MAESTRO_species_indices2266$acat)
MAESTRO_species_indices2266$scat<- gsub("0", "No_Info", MAESTRO_species_indices2266$scat)
MAESTRO_species_indices8855 <- MAESTRO_species_indices8855 %>% replace(is.na(.), 0)
MAESTRO_species_indices8855$ecat<- gsub("0", "No_Info", MAESTRO_species_indices8855$ecat)
MAESTRO_species_indices8855$vcat<- gsub("0", "No_Info", MAESTRO_species_indices8855$vcat)
MAESTRO_species_indices8855$acat<- gsub("0", "No_Info", MAESTRO_species_indices8855$acat)
MAESTRO_species_indices8855$scat<- gsub("0", "No_Info", MAESTRO_species_indices8855$scat)


MAESTRO_species_indices26 <- MAESTRO_species_indices26 %>% replace_na(list(RARE_tax_funct= "Mid2D"))

#MAESTRO_species_indices26 <- MAESTRO_species_indices26 %>% replace(is.na(.), 0)
#MAESTRO_species_indices85 <- MAESTRO_species_indices85 %>% replace(is.na(.), 0)

MAESTRO_species_indices266<- MAESTRO_species_indices26[,-1]
MAESTRO_species_indices266<- MAESTRO_species_indices266[,-4]
MAESTRO_species_indices266_S<-MAESTRO_species_indices266[,c(1:7)]
MAESTRO_species_indices266_A<-MAESTRO_species_indices266[,c(1:3,8:11)]
MAESTRO_species_indices266_E<-MAESTRO_species_indices266[,c(1:3,12:15)]
MAESTRO_species_indices266_SCORE<-MAESTRO_species_indices266[,c(1:3,16:19)]
MAESTRO_species_indices266_CAT<-MAESTRO_species_indices266[,c(1:3,20:23)]
#create pairs plot
pairs.panels(MAESTRO_species_indices266_S, main= "SSP1-2.6 Biodiversity and Vulnerability Sensitivity Correlations", stars = TRUE)
pairs.panels(MAESTRO_species_indices266_A, main= "SSP1-2.6 Biodiversity and Vulnerability Adaptivity Correlations", stars = TRUE)
pairs.panels(MAESTRO_species_indices266_E, main= "SSP1-2.6 Biodiversity and Vulnerability Exposure Correlations", stars = TRUE)
pairs.panels(MAESTRO_species_indices266_SCORE, main= "SSP1-2.6 Biodiversity and Vulnerability Scores Correlations", stars = TRUE)
pairs.panels(MAESTRO_species_indices266_CAT, main= "SSP1-2.6 Biodiversity and Vulnerability Categories Correlations", stars = TRUE)


###### I - use the vulnerability indices to color code the species inside the 3D plot ######


# 3D plot with colour by vulnerability
MAESTRO_species_indices226 <- merge(MAESTRO_species_indices4,vspec_MAESTRO26, by = "genus_sp",, all.x = TRUE  )
MAESTRO_species_indices226 <- MAESTRO_species_indices226 %>% replace(is.na(.), 0)



plot_3d_S.TSMr <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                y=~funct_dist_scaled,
                                                z=~evol_dist_scaled, alpha = 0.5, 
                                                color =~S.TSMr,marker_color=~S.TSMr, marker_colorscale='Viridis',
                                                size = ~S.TSMr,
                                                type = "scatter3d",
                                                mode= 'markers')


plot_3d_S.TSMr


#_____________
# SSP1-2.6: Climate Risk (circle colors) and Vulnerability Score (circle size)'
colors1 <- c('salmon', 'cyan', '#965F8A','grey' , 'red2')
MAESTRO_species_indices2266$vcat<-MAESTRO_species_indices2266$vcat

fig <- plot_ly(MAESTRO_species_indices2266, x=~scarcity_logbio_scaled,
               y=~funct_dist_scaled,
               z=~evol_dist_scaled, color = ~vcat, size = ~vuln, colors = colors1,
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 20),
               text = ~paste('Species:', genus_sp, '<br>Vulnerability score:', vuln, '<br>Climate risk score:', vcat))

fig <- fig %>% layout(title = 'SSP1-2.6: Climate Risk (circle colors) and Vulnerability Score (circle size)',
                      scene = list(xaxis = list(title = 'Taxonomic Rarity'
                                                ),
                                   yaxis = list(title = 'Functional Rarity'),
                                   zaxis = list(title = 'Evolutionary Rarity')))

fig
library(htmlwidgets)

saveWidget(fig, file = "myplot.html")
#_____________



#_____________
# SSP1-2.6: Scat: Sensitivity risk score
colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
MAESTRO_species_indices2266$scat<-MAESTRO_species_indices2266$scat

figA <- plot_ly(MAESTRO_species_indices2266, x=~scarcity_logbio_scaled,
               y=~funct_dist_scaled,
               z=~evol_dist_scaled, color = ~scat, size = ~sens, colors = colors1,
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 20),
               text = ~paste('Species:', genus_sp, '<br>Sensitivity Risk Score:', scat, '<br>Climate Sensitivity Score:', sens))

figA <- figA %>% layout(title = 'SSP1-2.6: Sensitivity Risk Score (circle colors) and Climate Sensitivity Score (circle size)',
                      scene = list(xaxis = list(title = 'Taxonomic Rarity'
                      ),
                      yaxis = list(title = 'Functional Rarity'),
                      zaxis = list(title = 'Evolutionary Rarity')))

figA


saveWidget(figA, file = "myplot_sensi.html")
#_____________

#_____________

colors1 <- c('salmon', 'cyan', 'grey' ,'red2' )
um<- ggplot(MAESTRO_species_indices2266, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=vcat)) + theme_bw() +
  geom_point(aes(size=vuln), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")+ ggtitle("SSP1-2.6: Climate Risk and Vulnerability Score") 


colors1 <- c('salmon', 'cyan', 'grey' ,'red2' )
dois<-ggplot(MAESTRO_species_indices8855, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=vcat)) + theme_bw() +
  geom_point(aes(size=vuln), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness") + ggtitle("SSP1-8.5: Climate Risk and Vulnerability Score") 

#_____________

colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
tres<- ggplot(MAESTRO_species_indices2266, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=scat)) + theme_bw() +
  geom_point(aes(size=sens), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")+ ggtitle("SSP1-2.6: Sensitivity Risk Score and Climate Sensitivity Score") 


colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
quatro<-ggplot(MAESTRO_species_indices8855, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=scat)) + theme_bw() +
  geom_point(aes(size=sens), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness") + ggtitle("SSP1-8.5: Sensitivity Risk Score and Climate Sensitivity Score") 

#_____________



#_____________

colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
cinco<- ggplot(MAESTRO_species_indices2266, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=ecat)) + theme_bw() +
  geom_point(aes(size=expo), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")+ ggtitle("SSP1-2.6: Exposure Risk Score and Climate Exposure Score") 


colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
seis<-ggplot(MAESTRO_species_indices8855, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=ecat)) + theme_bw() +
  geom_point(aes(size=expo), alpha= 0.5)+
  scale_color_manual(values = colors1 )+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]], colour = "blue")+
  geom_vline(xintercept=q.funct[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness") + ggtitle("SSP1-8.5: Exposure Risk Score and Climate Exposure Score") 

#_____________
ggarrange(um, dois, cinco,seis,
          labels = c("A", "B", "C", "D"),
          ncol = 2, nrow = 2)


write.csv(MAESTRO_species_indices8855, "MAESTRO_85_vulnerability_and_rarity.csv")
write.csv(MAESTRO_species_indices2266, "MAESTRO_26_vulnerability_and_rarity.csv")



colors1 <- c('salmon', 'cyan', 'grey',  'red2', '#965F8A' )
dez<-ggplot(MAESTRO_species_indices8855, aes(x=scarcity_logbio_scaled, y=evol_dist_scaled)) + theme_bw() +
  geom_point(alpha= 0.5)+
  scale_color_manual(values = colors1 )+
    geom_smooth(method=lm)+
  geom_vline(xintercept=q.scarcity[[3]], colour = "blue")+
  geom_vline(xintercept=q.scarcity[[1]],linetype = "dashed", colour = "red")+
  geom_hline(yintercept=q.evol[[3]], colour = "blue")+ 
  geom_hline(yintercept=q.evol[[1]],linetype = "dashed", colour = "red")+
  xlab("Scarcity")+
  ylab("Evolutinary Distinctiveness") + ggtitle("SSP1-8.5: Exposure Risk Score and Climate Exposure Score") 












#_____________
# S.TSMr: Thermal safety margin index SSP5-8.5
colors1 <- c('salmon', 'cyan', '#965F8A','grey' , 'red2')
MAESTRO_species_indices8855$vcat<-MAESTRO_species_indices8855$vcat

fig1 <- plot_ly(MAESTRO_species_indices8855, x=~scarcity_logbio_scaled,
               y=~funct_dist_scaled,
               z=~evol_dist_scaled, color = ~vcat, size = ~vuln, colors = colors1,
               marker = list(symbol = 'circle', sizemode = 'diameter'), sizes = c(5, 20),
               text = ~paste('Species:', genus_sp, '<br>Vulnerability score:', vuln, '<br>Climate risk score:', vcat))

fig1 <- fig1 %>% layout(title = 'SSP5-8.5: Climate Risk (circle colors) and Vulnerability Score (circle size)',
                      scene = list(xaxis = list(title = 'Taxonomic Rarity'
                      ),
                      yaxis = list(title = 'Functional Rarity'),
                      zaxis = list(title = 'Evolutionary Rarity')))

fig1

#_____________







#_____________
# S.TSMr: Thermal safety margin index SSP1-2.6
fig_AC.tvar <- plot_ly(MAESTRO_species_indices226, x=~scarcity_logbio_scaled,
               y=~funct_dist_scaled,
               z=~evol_dist_scaled,
               marker = list(color = ~AC.tvar, colorscale = c('#FFE1A1', '#683531'), showscale = TRUE))
fig_AC.tvar <- fig_AC.tvar %>% add_markers()
fig_AC.tvar <- fig_AC.tvar %>% layout(scene = list(xaxis = list(title = 'Taxonomic Rarity'),
                                   yaxis = list(title = 'Functional Rarity'),
                                   zaxis = list(title = 'Evolutionary Rarity')),
                      annotations = list(
                        x = 1.13,
                        y = 1.05,
                        text = 'AC.tvar',
                        xref = 'paper',
                        yref = 'paper',
                        showarrow = FALSE
                      ))
fig_AC.tvar

plotly::export(p = fig, #the graph to export
               file = "graph 1.png")


#_____________



plot_3d_AC.hrange <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                         y=~funct_dist_scaled,
                                                         z=~evol_dist_scaled, alpha = 0.5, color =~AC.hrange  )


plot_3d_AC.hrange


plot_3d_sens <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                            y=~funct_dist_scaled,
                                                            z=~evol_dist_scaled, alpha = 0.5, color =~sens  )


plot_3d_sens


plot_3d_expo <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                       y=~funct_dist_scaled,
                                                       z=~evol_dist_scaled, alpha = 0.5, color =~expo  )


plot_3d_expo

plot_3d_vuln <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_logbio_scaled,
                                                       y=~funct_dist_scaled,
                                                       z=~evol_dist_scaled, alpha = 0.5, color =~vuln , type = "scatter3d",
                                                       mode= 'markers' ) 

|> 
  add_markers() |> 
  layout(scene = list(xaxis = list(title = "Taxon. Scarcity"),
                      yaxis = list(title = "Funct. Dist."),
                      zaxis = list(title = "Evolu. Dist.")))


plot_3d_vuln


# ,
# marker = list(size = ~vuln*3)
# mode   = 'markers',


# select the coordinates from the species that we are working on  
species_coordinates_MAESTRO<- biomass_maestro %>% dplyr::select(genus_sp,  Lon,Lat)
write.csv(species_coordinates_MAESTRO, "species_coordinates_MAESTRO.csv")
write.csv(MAESTRO_species_indices3, "MAESTRO_species_indices3.csv")

MAESTRO_species_indices3$genus_sp<- rownames(MAESTRO_species_indices3)










try<- read.csv("try.csv", h = T )
try$FID<- as.factor(try$FID)


evol<-quantile(try$evol_dist_scaled, c(0.25, 0.5, 0.75), type=1)



try1<- try %>% 
  group_by(FID, species) %>% 
  dplyr::summarise("length.maturity"  = max(length.maturity),
                   "age.maturity"  = max(age.maturity),
                   "growth.coefficient"  = max(growth.coefficient),
                   "tl"  = max(tl),
                   "fecundity"  = max(fecundity),
                   "offspring.size"  = max(offspring.size),
                   "spawning.type"  = unique(spawning.type),
                   "habitat"  = unique(habitat),
                   "feeding.mode"  = unique(feeding.mode))



###### J - use the new AquaMaps list of species 

# Load AquaMaps data for the grid cells from MAESTRO
load("C:/Users/isaac/Desktop/MAESTRO/MAESTRO/Spp_DistsAcrossMAESTRO_FlatFile_BINARY.RData")
# File name - nrb

str(nrb)

#summary(nrb)

sort(unique(nrb$class))

# bony fish
aqua_bony_fish <- subset(nrb, class == "Actinopterygii"  )
aqua_bony_fish1 <- aqua_bony_fish[, c("lon", "lat", "cell", "latname")]

# cartilaginous fish
aqua_elasmo_fish <- subset(nrb, class == "Elasmobranchii"  )
aqua_elasmo_fish1 <- aqua_elasmo_fish[, c("lon", "lat", "cell", "latname")]


# join files
aqua_all<- rbind(aqua_bony_fish1,aqua_elasmo_fish1 )

aqua_all_species_list<-  unique(aqua_all$latname)

library(rfishbase)

# check the environment where the species occur, is there a species that occur only in freshwater systems?

aqua_all_species_fishbase_table<-  species(aqua_all_species_list)

# the list of species from AquaMaps does not present species that occur only in freshwater
# we can keep all the species

colnames(aqua_all)[colnames(aqua_all) == 'latname'] <- 'Species'

write.csv(aqua_all, "aqua_all.csv")

length(unique(aqua_all$Species))
# [1] 1336
length(unique(aqua_all$cell))
# [1] 12649

aqua_all$Species <- gsub(" ", "_", aqua_all$Species)

# aqua_species11$Species <- gsub(" ", "_", aqua_species11$Species)

# considering that all species occur in the sea we do not need to subset
# aqua_all2<- merge(aqua_all,aqua_species11, by = "Species", all.x = TRUE )

aqua_all2<- aqua_all

aqua_all2$cell<- as.factor(aqua_all2$cell)
aqua_all2$ID<-row.names(aqua_all2)
str(aqua_all2)
aqua_all2$ID<- as.factor(aqua_all2$ID)

aqua_all3<- aqua_all2 %>% dplyr::select(Species, ID,cell, lon, lat)
# aqua_all3$ID<-row.names(aqua_all3)

aqua_all3<- aqua_all3[,-1]

####### MAESTRO biomass data #######
# some of the points fell over the continents 
# below I removed a few of those points
# I also noticed the presence of a single species that occur only in frehswater
# Maybe it was a species identification error?
# the name of the species is Alosa agone

# prepare file to bring to ArcGIS, that is, ID lat and long

biomass_maestro6<- biomass_maestro %>% dplyr::select(ID,Lon, Lat )

# bring both files to ArcGIS
# subset the aquamaps file using the file from MAESTRO
# MAESTRO file is the list of coordinates where sampling occurred
# AquaMaps file is the list of coordinates with species occurrences
# Arnaud wants to use only grid cells where sampling occurred

write.csv(aqua_all3, "aqua_all3.csv")
write.csv(biomass_maestro6, "biomass_maestro6.csv")

# load the files saved from ArcGIS sub setting
# load the subset from MAESTRO -
# a few points fell over continents, here are only the points that are on the ocean
MAESTRO_subset <- read.csv("MAESTRO_subset_marineonly1.txt", header=TRUE, sep = ',')
MAESTRO_subset <- MAESTRO_subset[, -1]
MAESTRO_subset <- MAESTRO_subset[, -1]

colnames_maestro<- c("ID", "Lon", "Lat", "lon1", "lat1", "ID_shape")

colnames(MAESTRO_subset)  <- colnames_maestro
MAESTRO_subset$cell<- paste(MAESTRO_subset$lon1, MAESTRO_subset$lat1, sep="_")
# load the subset from AquaMaps - 
# as Arnaud suggested, the aim is to use the same grid cells that we have abundance data
# both files now have the same the grid cell IDs
AquaMaps_subset <- read.csv("AquaMaps_new_subset.txt", header=TRUE, sep = ',')
str(AquaMaps_subset)
AquaMaps_subset$ID<- as.factor(AquaMaps_subset$ID)

AquaMaps_subset1<- merge(AquaMaps_subset,aqua_all2, by = "ID", all.x = TRUE )

AquaMaps_species_list<- unique(AquaMaps_subset1$Species)

# create a file with species and cells
# this file is important
# it will be used to compute taxonomic restrictedness
AquaMaps_subset12<- AquaMaps_subset1 %>% dplyr::select(cell.x,Species )
AquaMaps_subset12<- unique(AquaMaps_subset12)
AquaMaps_subset12_colnames<-c("cell", "genus_sp")
colnames(AquaMaps_subset12)<- AquaMaps_subset12_colnames
write.csv(AquaMaps_subset12, "AquaMaps_pres_absc_long.csv")

######### MAESTRO file sub setting starts here

# merge the subset file with the biomass file again
MAESTRO_subset1<- merge(MAESTRO_subset ,biomass_maestro, by = "ID", all.x = TRUE )
# this step above was done to make sure all MAESTRO species are marine only

# get the new list of species from MAESTRO
MAESTRO_species_list<- unique(MAESTRO_subset1$genus_sp)
MAESTRO_species_list <- gsub("_", " ", MAESTRO_species_list)

library(rfishbase)
fb_tables<- fb_tables(server = c("fishbase"), version = "latest")
sort(fb_tables)
MAESTRO_estimate<- estimate(MAESTRO_species_list)
#MAESTRO_family2<-  rfishbase::fishbase(MAESTRO_species_list)
MAESTRO_ecology2<-  ecology(MAESTRO_species_list)
# check the presence of species that are restricted to freshwater using the information from FishBase
MAESTRO_species2<-  species(MAESTRO_species_list)
# Alosa agone 
# FishBASE distribution INFO about the species
# Europe: Lakes Como, Garda, Orta, Maggiore, Lugano and Iseo (northern Italy, Switzerland). 
# Introduced in Lakes Bolsena, Bracciano and Vico (central Italy).

# Removing species from the MAESTRO file here
MAESTRO_subset2<- subset(MAESTRO_subset1 , !genus_sp == "Alosa_agone" )  
# Acipenser sturio 
MAESTRO_subset2<- subset(MAESTRO_subset2 , !genus_sp == "Acipenser_sturio" )  

species_families_colnames<- c("genus_sp", "Family", "Class", "Order")
colnames(species_families)<-species_families_colnames 

str(MAESTRO_subset2)

MAESTRO_subset2<- merge(MAESTRO_subset2,species_families, by = "genus_sp", all.x = TRUE )

no_family1 <-   MAESTRO_subset2[MAESTRO_subset2$Class == "",]

MAESTRO_subset9<- MAESTRO_subset2[!is.na(MAESTRO_subset2$Family),]

MAESTRO_subset2<- MAESTRO_subset9

MAESTRO_subset2222<- MAESTRO_subset2 %>% 
  group_by(Class, Order, Family) %>% 
  dplyr::summarise("species_count"  = length(unique(genus_sp)),
                   "occurrences"  = length(genus_sp)
                                      )
# removing jawless fish 
# we need to keep only Teleosts and Elasmobranchs
MAESTRO_subset2<- subset(MAESTRO_subset2 , !Class == "Myxini" )  
MAESTRO_subset2<- subset(MAESTRO_subset2 , !Class == "Holocephali" ) 
MAESTRO_subset2<- subset(MAESTRO_subset2 , !Class == "Petromyzonti" ) 

MAESTRO_subset3333<- MAESTRO_subset2 %>% 
  group_by(Class, Order, Family) %>% 
  dplyr::summarise("species_count"  = length(unique(genus_sp)),
                   "occurrences"  = length(genus_sp)
  )

# removing other freshwater based order
# Cypriniformes
# freshwater_species<- subset(MAESTRO_subset2, Order== "Cypriniformes")
MAESTRO_subset2<- subset(MAESTRO_subset2 , !Order== "Cypriniformes") 

MAESTRO_subset3333<- MAESTRO_subset2 %>% 
  group_by(Class, Order, Family) %>% 
  dplyr::summarise("species_count"  = length(unique(genus_sp)),
                   "occurrences"  = length(genus_sp)
  )

MAESTRO_species_list6<- unique(MAESTRO_subset2$genus_sp)
MAESTRO_species_list6 <- gsub("_", " ", MAESTRO_species_list6)
MAESTRO_species2<-  species(MAESTRO_species_list6)
write.csv(MAESTRO_subset2, "MAESTRO_biomass_assemblage.csv")
# MAESTRO_species33<- subset(MAESTRO_species22, Brack == 0 ) 
# aqua_species11 <- aqua_species1[, c("Species", "BodyShapeI", "DemersPelag")]

######### MAESTRO file subseting ends here


## subset our AquaMaps tree here

AquaMaps_species_list <- gsub(" ", "_", AquaMaps_species_list)
subtree_AquaMaps <- castor::get_subtree_with_tips(jawlessfish_sharks_bonyfish_MAESTRO,
                                                  only_tips=AquaMaps_species_list)$subtree

#ggtree(subtree, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Jawless, Bony and Cartilaginous fish")

is.binary(subtree_AquaMaps)
#[1] TRUE
subtree_AquaMaps <- phytools::force.ultrametric(subtree_AquaMaps)
is.ultrametric(subtree_AquaMaps)
#[1] TRUE

saveRDS(subtree_AquaMaps, "subtree_AquaMaps.RDS")
## subset our MAESTRO tree here
MAESTRO_species_list6 <- gsub(" ", "_", MAESTRO_species_list6)
subtree_MAESTRO6 <- castor::get_subtree_with_tips(jawlessfish_sharks_bonyfish_MAESTRO, 
                                                  only_tips=MAESTRO_species_list6)$subtree

#ggtree(subtree, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Jawless, Bony and Cartilaginous fish")

is.binary(subtree_MAESTRO6)
#[1] TRUE
subtree_MAESTRO6 <- phytools::force.ultrametric(subtree_MAESTRO6)
is.ultrametric(subtree_MAESTRO6)
#[1] TRUE

saveRDS(subtree_MAESTRO6, "subtree_MAESTRO6.RDS")
# lets compare the list of species between MAESTRO and AquaMaps

AQUAMAPS_species_in_phylogeny <- subtree_AquaMaps$tip.label  
MAESTRO_species_in_phylogeny <- subtree_MAESTRO6$tip.label  

#display items that are in both vectors
species_in_common<- intersect(AQUAMAPS_species_in_phylogeny, MAESTRO_species_in_phylogeny )

#display items that are only in first vector, but not in second vector
AQUAMAPS_only<- setdiff(AQUAMAPS_species_in_phylogeny, MAESTRO_species_in_phylogeny )
MAESTRO_only<- setdiff( MAESTRO_species_in_phylogeny, AQUAMAPS_species_in_phylogeny)

write.csv(MAESTRO_only, "MAESTRO_only.txt")
write.csv(AQUAMAPS_only, "AQUAMAPS_only.txt")
write.csv(species_in_common, "species_in_common.txt")

library(VennDiagram)       

draw.pairwise.venn(area1=length(AQUAMAPS_species_in_phylogeny), 
                   area2=length(MAESTRO_species_in_phylogeny),
                   cross.area=length(species_in_common), 
                   category=c("AquaMaps","MAESTRO"),fill=c("blue","green"))

draw.pairwise.venn(length(AQUAMAPS_species_in_phylogeny), 
                   length(MAESTRO_species_in_phylogeny), 
                   length(species_in_common),
                   category = c("AquaMaps","MAESTRO"), 
                   lty = rep("blank", 2), ext.text = FALSE,
                   fill = c("blue","green"), 
                   alpha = rep(0.5, 2), 
                   cat.dist = rep(0.025, 2))
###### PLOT TREES

aqua_summary_sp<- aqua_all2 %>% 
  group_by(Species) %>% 
  dplyr::summarise("number_grids"  = length(unique(cell)))


aqua_summary_grids<- aqua_all2 %>% 
  group_by(cell) %>% 
  dplyr::summarise("number_species"  = length(unique(Species)))


aqua_grid_hist<- hist(aqua_summary_sp$number_grids,xlab = "Half Degree Grids", main = "Number of grid cells per species", breaks = seq(min(aqua_summary_sp$number_grids), max(aqua_summary_sp$number_grids), length.out = 30))
aqua_sp_hist<- hist(aqua_summary_grids$number_species,xlab = "Number of Species", main = "Number of species per grid cells", breaks = seq(min(aqua_summary_grids$number_species), max(aqua_summary_grids$number_species), length.out = 30))

aqua_few<- subset(aqua_summary_grids, number_species < 10 ) 

subtitle1<- length(subtree_AquaMaps$tip.label)

aqua_tree1 <-ggtree(subtree_AquaMaps, color="grey", size=0.001, ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "AquaMaps Tree",subtitle = "1245 Species")

aqua_tree1_circular <-ggtree(subtree_AquaMaps, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "AquaMaps Tree",subtitle ="1245 Species" )

aqua_tree1
aqua_tree1_circular

subtitle11<- length(subtree_MAESTRO6$tip.label)

MAESTRO_tree1 <-ggtree(subtree_MAESTRO6, color="grey", size=0.001, ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO Tree",subtitle = "591 Species")

MAESTRO_tree1_circular <-ggtree(subtree_MAESTRO6, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO Tree",subtitle = "591 Species")


MAESTRO_tree1
MAESTRO_tree1_circular







# evolutionary distinctiveness
aqua_tree.cm <- clade.matrix(aqua_tree)
aqua_Evol_Di <- ed.calc(aqua_tree.cm)
head(Evol_Di$spp)
most_distinct_species<-aqua_Evol_Di$spp[order(aqua_Evol_Di$spp$ED,decreasing=T),];head(most_distinct_species)
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
most_distinct_species$dist_scaled <- scale_values(most_distinct_species$ED)

write.csv(most_distinct_species, "aqua_most_distinct_species1.csv")


