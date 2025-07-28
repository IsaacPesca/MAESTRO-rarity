# this file below has the grid cell based on the fishing effort grid cells
MAESTRO_pres_absc_long1<-read.csv("MAESTRO_pres_absc_long.csv", h=T)
MAESTRO_pres_absc_long1<- MAESTRO_pres_absc_long1[,-1]
# MAESTRO_pres_absc_long1<-MAESTRO_pres_absc_long1[!is.na(MAESTRO_pres_absc_long1$cell),]

library(reshape2)

# transform the long format to a wide format
# this is the presence absence matrix from AquaMap
# correct column name
cname<- c("ID_f", "genus_sp")
colnames(MAESTRO_pres_absc_long1)<- cname

MAESTRO_species_list<- MAESTRO_species_indices$genus_sp
# keep only species with phylogeny information
# create object with list of species without phylogeny information
MAESTRO_only <- setdiff(MAESTRO_pres_absc_long1$genus_sp, MAESTRO_species_indices$genus_sp )
# use object to remove from the list
MAESTRO_pres_absc_long16 <- MAESTRO_pres_absc_long1[!MAESTRO_pres_absc_long1$genus_sp %in% MAESTRO_only, ]

# importante
# traits_MAESTRO666_subset <- MAESTRO_pres_absc_long1[rownames(MAESTRO_pres_absc_long1) %in% MAESTRO_most_distinct_species$species, ] 

MAESTRO_pres_absc_wide<- dcast(MAESTRO_pres_absc_long16, ID_f~genus_sp, length)

# make sure that it is a matrix
typeof(MAESTRO_pres_absc_wide)

rownames(MAESTRO_pres_absc_wide)<-MAESTRO_pres_absc_wide$ID_f

MAESTRO_pres_absc_wide1<-MAESTRO_pres_absc_wide[,-1]

MAESTRO_pres_absc_wide2<- data.frame(MAESTRO_pres_absc_wide1)
MAESTRO_pres_absc_wide2<- as.matrix(MAESTRO_pres_absc_wide1)