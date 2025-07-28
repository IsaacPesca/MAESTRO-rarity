###### J - use the new AquaMaps list of species 

# load MAESTRO file
biomass_maestro <- read.csv("datras_medits_biomass_withtraits.csv", dec = ",", sep = ";")
length(unique(biomass_maestro$genus_sp))
# 620 species
# add a column with a new ID to use at the ArcGIS
biomass_maestro$ID_m<-row.names(biomass_maestro)
species_maestro<- unique(biomass_maestro$genus_sp)
species_maestro<- as.character(species_maestro)

# Load AquaMaps data for the grid cells from MAESTRO
# first file shared by Daniel Boyce
load("C:/Users/isaac/Desktop/MAESTRO/MAESTRO/Spp_DistsAcrossMAESTRO_FlatFile_BINARY.RData")
# second file shared by Daniel Boyce, with probabilities
load("C:/Users/isaac/Desktop/MAESTRO/MAESTRO/AquaMaps_MAESTROSpp_Probs.RData")


# File name - nrb
str(nrb)

nrb_sp_id<- nrb %>% dplyr::select(speciesid, latname, class) 

nrb_sp_id<- unique(nrb_sp_id)

# File name - nr
str(nr)
# now lets subset occurrences with probability higher than 0.9
AquaMaps_high_probability<- subset(nr, prob > 0.7)

nrb2<- merge(AquaMaps_high_probability,nrb_sp_id, by = "speciesid", all.x = TRUE )

# bony fish
aqua_bony_fish <- subset(nrb2, class == "Actinopterygii"  )
# aqua_bony_fish1 <- aqua_bony_fish[, c("lon", "lat", "cell", "latname")]

# cartilaginous fish
aqua_elasmo_fish <- subset(nrb2, class == "Elasmobranchii"  )
# aqua_elasmo_fish1 <- aqua_elasmo_fish[, c("lon", "lat", "cell", "latname")]


# now we can unload the nrb file because it is occupy a lot of space on the global environment
# nrb<- NULL

# join files
aqua_all<- rbind(aqua_bony_fish,aqua_elasmo_fish )

aqua_all_species_list<-  unique(aqua_all$latname)

library(rfishbase)

# check the environment where the species occur, is there a species that occur only in freshwater systems?

aqua_all_species_fishbase_table<-  species(aqua_all_species_list)

# the list of species from AquaMaps does not present species that occur only in freshwater
# we can keep all the species

colnames(aqua_all)[colnames(aqua_all) == 'latname'] <- 'Species'

# write.csv(aqua_all, "aqua_all.csv")

length(unique(aqua_all$Species))
# [1] 513
length(unique(aqua_all$cell))
# [1] 11096

aqua_all$Species <- gsub(" ", "_", aqua_all$Species)

# aqua_species11$Species <- gsub(" ", "_", aqua_species11$Species)

# considering that all species occur in the sea we do not need to subset
# aqua_all2<- merge(aqua_all,aqua_species11, by = "Species", all.x = TRUE )

aqua_all2<- aqua_all

aqua_all2$cell<- as.factor(aqua_all2$cell)
aqua_all2$ID_a<-row.names(aqua_all2)
str(aqua_all2)
# aqua_all2$ID<- as.factor(aqua_all2$ID)

aqua_all3<- aqua_all2 %>% dplyr::select(ID_a, cell, lon, lat) # 
# aqua_all3$ID<-row.names(aqua_all3)


####### MAESTRO biomass data #######
# some of the points fell over the continents 
# below I removed a few of those points
# I also noticed the presence of a single species that occur only in freshwater
# Maybe it was a species identification error?
# the name of the species is Alosa agone

# prepare file to bring to ArcGIS, that is, ID lat and long

biomass_maestro6<- biomass_maestro %>% dplyr::select(ID_m,Lon, Lat )#  

# bring both files to ArcGIS
# subset the aquamaps file using the file from MAESTRO
# MAESTRO file is the list of coordinates where sampling occurred
# AquaMaps file is the list of coordinates with species occurrences
# Arnaud wants to use only grid cells where sampling occurred
# aqua_all3<- aqua_all3[,-2]
write.csv(aqua_all3, "aqua_all5.csv")
write.csv(biomass_maestro6, "biomass_maestro7.csv")

# load the files saved from ArcGIS sub setting
# load the subset from MAESTRO -
# a few points fell over continents, here are only the points that are on the ocean
MAESTRO_subset <- read.csv("MAESTRO_ID_m_ID_f.txt", header=TRUE, sep = ',')
MAESTRO_subset <- MAESTRO_subset[, -1]
MAESTRO_subset <- MAESTRO_subset[, -1]
# this file has the cell identification for MAESTRO biomass and AquaMaps grid cells
# it is useful to combine with the file that has the grid cells identification from the fishing efforts file
# atencao aqui com o ID
# colnames_maestro<- c("ID", "Lon", "Lat", "lon1", "lat1", "ID_shape")

# colnames(MAESTRO_subset)  <- colnames_maestro
# MAESTRO_subset$cell<- paste(MAESTRO_subset$lon1, MAESTRO_subset$lat1, sep="_")
#MAESTRO_subset<- MAESTRO_subset[,-6]
#MAESTRO_subset<- MAESTRO_subset[,-1]
# load the subset from AquaMaps - 
# as Arnaud suggested, the aim is to use the same grid cells that we have abundance data
# both files now have the same the grid cell IDs
#AquaMaps_subset <- read.csv("AquaMaps_new_subset.txt", header=TRUE, sep = ',')

# THIS FILE BELOW HAS THE CELL1 FROM THE FISHING EFFORT FILE
# THIS FILE HAS THE SUBSET FROM THE FISHING EFFORT FILE AND AQUAMAPS SUBSET
AquaMaps_subset <- read.csv("AquaMaps_subset_ID_a_ID_f.txt", header=TRUE, sep = ',')
AquaMaps_subset$cell<- paste(AquaMaps_subset$XCoord, AquaMaps_subset$YCoord, sep="_")
AquaMaps_subset <- AquaMaps_subset[, -1]
AquaMaps_subset <- AquaMaps_subset[, -1]


# cell1 is the identifier from fishing efforts grid cells
# cell1 is the identifier from the smaller AquaMaps grid cells
# each grid cell from fishing efforts file is 4 times larger than a grid cell from AquaMaps
# the subset from fishing effort grids was done using ArcGIS

# fisheffort_colnames<- c("Lat", "Lon", "cell1", "cell")
# colnames(AquaMaps_subset)<- fisheffort_colnames


AquaMaps_subset$ID_a <- as.factor(AquaMaps_subset$ID_a)
AquaMaps_subset$ID_f <- as.factor(AquaMaps_subset$ID_f)
AquaMaps_subset$cell <- as.factor(AquaMaps_subset$cell)
str(AquaMaps_subset)
# AquaMaps_subset$cell1<- as.factor(AquaMaps_subset$cell1)
# aqua_all2$cell<- as.factor(aqua_all2$cell)
# aqua_all2<-aqua_all2[,-2]
# aqua_all2<-aqua_all2[,-1]

# AquaMaps_subset<- unique(AquaMaps_subset)
# str(AquaMaps_subset)
aqua_all2$ID_a<-as.factor(aqua_all2$ID_a)

str(aqua_all2)

AquaMaps_subset6<-AquaMaps_subset[,-1]
AquaMaps_subset6<-unique(AquaMaps_subset6)
AquaMaps_subset1<- merge(AquaMaps_subset6,aqua_all2, by = "cell", all = TRUE )
library(tidyr)
AquaMaps_subset1<- AquaMaps_subset1 %>% drop_na(ID_f)
write.csv(AquaMaps_subset1, "AquaMaps_subset66666.csv")

sum(is.na(AquaMaps_subset1$ID_f))
AquaMaps_species_list<- unique(AquaMaps_subset1$Species)

# create a file with species and cells
# this file is important
# it will be used to compute taxonomic restrictedness
AquaMaps_subset12<- AquaMaps_subset1 %>% dplyr::select(ID_f,Species )
AquaMaps_subset12<- unique(AquaMaps_subset12)
# AquaMaps_subset12_colnames<-c("cell", "genus_sp")
# colnames(AquaMaps_subset12)<- AquaMaps_subset12_colnames
# AquaMaps_subset12<- subset(AquaMaps_subset12 , !Species == "Acipenser_sturio" )
write.csv(AquaMaps_subset12, "AquaMaps_pres_absc_long.csv")

######### MAESTRO file sub setting starts here

# merge the subset file with the biomass file again
# ATENCAO AQUI COM O ID
MAESTRO_subset$ID_m<- as.factor(MAESTRO_subset$ID_m)
MAESTRO_subset$ID_f<- as.factor(MAESTRO_subset$ID_f)
biomass_maestro$ID_m<- as.factor(biomass_maestro$ID_m)
MAESTRO_subset1<- merge(MAESTRO_subset ,biomass_maestro, by = "ID_m", all.x = TRUE )
# this step above was done to make sure all MAESTRO species are marine only

# get the new list of species from MAESTRO
MAESTRO_species_list<- unique(MAESTRO_subset1$genus_sp)
MAESTRO_species_list <- gsub("_", " ", MAESTRO_species_list)

library(rfishbase)
fb_tables<- fb_tables(server = c("fishbase"), version = "latest")
sort(fb_tables)
#MAESTRO_estimate<- estimate(MAESTRO_species_list)
#MAESTRO_family2<-  rfishbase::fishbase(MAESTRO_species_list)
#MAESTRO_ecology2<-  ecology(MAESTRO_species_list)
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

species_families<- read.csv("species_families.csv",h=T)
species_families<-species_families[,-1]
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

# this file below it is useful for finding where the rare species occur
# it is just a presence absence file in a long format
# it is similar to the AquaMaps presence absence file in a long format
MAESTRO_subset223<- MAESTRO_subset2 %>% dplyr::select(ID_f,genus_sp ) 
MAESTRO_subset2233<- unique(MAESTRO_subset223)
write.csv(MAESTRO_subset2233, "MAESTRO_subset_pres_absc_long.csv")

# this file below has the grid cell based on the fishing effort grid cells
MAESTRO_pres_absc_long1<-read.csv("MAESTRO_subset_pres_absc_long.csv", h=T)
MAESTRO_pres_absc_long1<- MAESTRO_pres_absc_long1[,-1]
# MAESTRO_pres_absc_long1<-MAESTRO_pres_absc_long1[!is.na(MAESTRO_pres_absc_long1$cell),]

library(reshape2)

# transform the long format to a wide format
# this is the presence absence matrix from AquaMap
# correct column name
# cname<- c("ID_f", "genus_sp")
# colnames(MAESTRO_pres_absc_long1)<- cname

# MAESTRO_species_list<- MAESTRO_species_indices$genus_sp
# keep only species with phylogeny information
# create object with list of species without phylogeny information
# MAESTRO_only <- setdiff(MAESTRO_pres_absc_long1$genus_sp, MAESTRO_species_indices$genus_sp )
# use object to remove from the list
# MAESTRO_pres_absc_long16 <- MAESTRO_pres_absc_long1[!MAESTRO_pres_absc_long1$genus_sp %in% MAESTRO_only, ]

# importante
# traits_MAESTRO666_subset <- MAESTRO_pres_absc_long1[rownames(MAESTRO_pres_absc_long1) %in% MAESTRO_most_distinct_species$species, ] 

MAESTRO_pres_absc_wide<- dcast(MAESTRO_pres_absc_long1, ID_f~genus_sp, length)

# make sure that it is a matrix
typeof(MAESTRO_pres_absc_wide)

rownames(MAESTRO_pres_absc_wide)<-MAESTRO_pres_absc_wide$ID_f

MAESTRO_pres_absc_wide1<-MAESTRO_pres_absc_wide[,-1]

MAESTRO_pres_absc_wide2<- data.frame(MAESTRO_pres_absc_wide1)
MAESTRO_pres_absc_wide2<- as.matrix(MAESTRO_pres_absc_wide1)

write.csv(MAESTRO_pres_absc_wide2, "MAESTRO_pres_absc_wide2.csv")


# this file below is useful to measure scarcity using the MAESTRO data
MAESTRO_species_list6<- unique(MAESTRO_subset2$genus_sp)
MAESTRO_species_list6 <- gsub("_", " ", MAESTRO_species_list6)
MAESTRO_species2<-  species(MAESTRO_species_list6)
write.csv(MAESTRO_subset2, "MAESTRO_biomass_assemblage.csv")
# MAESTRO_species33<- subset(MAESTRO_species22, Brack == 0 ) 
# aqua_species11 <- aqua_species1[, c("Species", "BodyShapeI", "DemersPelag")]

######### MAESTRO file subseting ends here


## subset our AquaMaps tree here
# first load the mega tree with jawless bony and cartilaginous fish species
jawlessfish_sharks_bonyfish_MAESTRO<- readRDS("jawlessfish_sharks_bonyfish_MAESTRO.RDS")
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

#### what is the overlap between the species with phylogenies, AquaMaps data and data on MAESTRO??

draw.pairwise.venn(length(AQUAMAPS_species_in_phylogeny), 
                   length(MAESTRO_species_in_phylogeny), 
                   length(species_in_common),
                   category = c("AquaMaps","MAESTRO"), 
                   lty = rep("blank", 2), ext.text = FALSE,
                   fill = c("blue","green"), 
                   alpha = rep(0.5, 2), 
                   cat.dist = rep(0.025, 2))
#### 
#### how many species per grid cell??
#### how many grid cell per species??


par(mfrow = c(2, 2))

aqua_summary_sp<- MAESTRO_subset2233 %>% 
  group_by(genus_sp) %>% 
  dplyr::summarise("number_grids"  = length(unique(ID_f)))


aqua_summary_grids<- MAESTRO_subset2233 %>% 
  group_by(ID_f) %>% 
  dplyr::summarise("number_species"  = length(unique(genus_sp)))


aqua_grid_hist<- hist(aqua_summary_sp$number_grids,xlab = "Grids", main = "MAESTRO - Number of grid cells per species", breaks = seq(min(aqua_summary_sp$number_grids), max(aqua_summary_sp$number_grids), length.out = 30))
aqua_sp_hist<- hist(aqua_summary_grids$number_species,xlab = "Number of Species", main = "MAESTRO - Number of species per grid cells", breaks = seq(min(aqua_summary_grids$number_species), max(aqua_summary_grids$number_species), length.out = 30))


aqua_maps_summary_sp<-  AquaMaps_subset12 %>% 
  group_by(Species) %>% 
  dplyr::summarise("number_grids"  = length(unique(ID_f)))


aqua_maps_summary_grids<-  AquaMaps_subset12 %>% 
  group_by(ID_f) %>% 
  dplyr::summarise("number_species"  = length(unique(Species)))


aqua_maps_grid_hist<- hist(aqua_maps_summary_sp$number_grids,xlab = "Grids", main = "AquaMaps - Number of grid cells per species", breaks = seq(min(aqua_maps_summary_sp$number_grids), max(aqua_maps_summary_sp$number_grids), length.out = 30))
aqua_maps_sp_hist<- hist(aqua_maps_summary_grids$number_species,xlab = "Number of Species", main = "AquaMaps - Number of species per grid cells", breaks = seq(min(aqua_maps_summary_grids$number_species), max(aqua_maps_summary_grids$number_species), length.out = 30))

par(mfrow = c(2, 2))


