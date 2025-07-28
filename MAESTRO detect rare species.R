
# Code Used to detect what are the rare species and where those species occur around Europe
# What are the grid cells that those species occur? This code provides this answer, 
# then counts the number of rare species that occur inside each grid cell
# This code also binds the fishing effort information with the matrix that has the number of rare species per grid cell
# Isaac Trindade Santos
# MAESTRO workgroup

### top quartile for both restrictedness and distinctiveness (GDR) and 
### restrictedness and uniqueness (UR). This returns a data frame with the 
### total species per grid cell plus number of UR and GDR species

# load species indices from MAESTRO

# aqui eh novo 2025
MAESTRO_species_indices1<- read.csv("MAESTRO_species_indices123.csv", h= T)
MAESTRO_species_indices1<- MAESTRO_species_indices1[,-1]

str(MAESTRO_species_indices1)
length(unique(MAESTRO_species_indices1$genus_sp))
# 555
MAESTRO_species_indices1$genus_sp<- as.factor(MAESTRO_species_indices1$genus_sp)

MAESTRO_subset_pres_absc_long<-read.csv("MAESTRO_subset_pres_absc_long.csv",h= T )
MAESTRO_subset_pres_absc_long<- MAESTRO_subset_pres_absc_long[,-1]
str(MAESTRO_subset_pres_absc_long)

length(unique(MAESTRO_subset_pres_absc_long$genus_sp))
# 555
MAESTRO_subset_pres_absc_long$genus_sp<- as.factor(MAESTRO_subset_pres_absc_long$genus_sp)

# merge both files
MAESTRO_species_indices666 <- merge(MAESTRO_subset_pres_absc_long,
                                   MAESTRO_species_indices1,
                                   by = "genus_sp", all.x = TRUE )

length(unique(MAESTRO_species_indices666$OID_First))
#[1] 2313
length(unique(MAESTRO_species_indices666$ID_f))
#[1] 1177


ID_f_OID_First<- MAESTRO_species_indices666 %>% dplyr::select(OID_First,ID_f )

# calculate the quantiles from each index

qe<-quantile(MAESTRO_species_indices1$evol_dist_scaled, c(0.25, 0.5, 0.75), type=1)
qf<-quantile(MAESTRO_species_indices1$funct_dist_scaled, c(0.25, 0.5, 0.75), type=1)
qt<-quantile(MAESTRO_species_indices1$scarcity_logbio_scaled, c(0.25, 0.5, 0.75), type=1)
qr<-quantile(MAESTRO_species_indices1$restrictedness_scaled, c(0.25, 0.5, 0.75), type=1)


# subset the species that are rare in the three dimensions at the same time

rare3d_scarce<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]]  & 
                 MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] &
                 MAESTRO_species_indices666$scarcity_logbio_scaled>qt[[3]])

rare3d_restrict<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]]  & 
                 MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] &
                 MAESTRO_species_indices666$restrictedness_scaled>qr[[3]])

rare4d_restrict_scarce<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]]  & 
                          MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] &
                          MAESTRO_species_indices666$restrictedness_scaled>qr[[3]]&
                            MAESTRO_species_indices666$scarcity_logbio_scaled>qt[[3]])


rare3d_scarce_1<- rare3d_scarce[,c(-2,-3)]
rare3d_scarce_1<- unique(rare3d_scarce_1)
write.csv(rare3d_scarce_1, "MAESTRO_list_rare_species_three_dimensions_scarce.csv")


# subset the species that are rare functionally and phylogenetically

rare2d_PD_FD<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]]  & 
                       MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] )
rare2d_PD_FD_1<- rare2d_PD_FD[,c(-2,-3)]
rare2d_PD_FD_1<- unique(rare2d_PD_FD_1)
write.csv(rare2d_PD_FD_1, "MAESTRO_list_rare_species_phylogenetic_functional.csv")

# subset the species that are rare functionally and are scarce

rare2d_FD_TD_scarce<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] &
                       MAESTRO_species_indices666$scarcity_logbio_scaled>qt[[3]])

rare2d_FD_TD_scarce_1<- rare2d_FD_TD_scarce[,c(-2,-3)]
rare2d_FD_TD_scarce_1<- unique(rare2d_FD_TD_scarce_1)
write.csv(rare2d_FD_TD_scarce_1, "MAESTRO_list_rare_species_taxonomic_functional_scarce.csv")

# subset the species that are rare functionally and are restrict

rare2d_FD_TD_restrict<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$funct_dist_scaled>qf[[3]] &
                       MAESTRO_species_indices666$restrictedness_scaled>qr[[3]])

rare2d_FD_TD_restrict_1<- rare2d_FD_TD_restrict[,c(-2,-3)]
rare2d_FD_TD_restrict_1<- unique(rare2d_FD_TD_restrict_1)
write.csv(rare2d_FD_TD_restrict_1, "MAESTRO_list_rare_species_taxonomic_functional_restrict.csv")

# subset the species that are rare phylogenetically and are scarce

rare2d_PD_TD_scarce<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]] &
                       MAESTRO_species_indices666$scarcity_logbio_scaled>qt[[3]])
rare2d_PD_TD_scarce_1<- rare2d_PD_TD_scarce[,c(-2,-3)]
rare2d_PD_TD_scarce_1<- unique(rare2d_PD_TD_scarce_1)
write.csv(rare2d_PD_TD_scarce_1, "MAESTRO_list_rare_species_taxonomic_phylogenetic_scarce.csv")


# subset the species that are rare phylogenetically and are restrict

rare2d_PD_TD_restrict<-subset(MAESTRO_species_indices666, MAESTRO_species_indices666$evol_dist_scaled>qe[[3]] &
                       MAESTRO_species_indices666$restrictedness_scaled>qr[[3]])
rare2d_PD_TD_restrict_1<- rare2d_PD_TD_restrict[,c(-2,-3)]
rare2d_PD_TD_restrict_1<- unique(rare2d_PD_TD_restrict_1)
write.csv(rare2d_PD_TD_restrict_1, "MAESTRO_list_rare_species_taxonomic_phylogenetic_restrict.csv")




## get the species totals per grid cell

##mudei aqui
getNumSp<-as.data.frame(MAESTRO_subset_pres_absc_long %>% group_by(OID_First) %>%
                          summarise(spNum=n_distinct(genus_sp)))


## group by grid cell and total the rare (as above) species in each

## count the number of rare species in each grid cell
##  PD and FD
get_rare2d_PD_FD<-as.data.frame(rare2d_PD_FD %>% group_by(OID_First) %>%
                         summarise(rare2d_PD_FD=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD, PD and FD scarce
get_rare3d_scarce<-as.data.frame(rare3d_scarce %>% group_by(OID_First) %>%
                                  summarise(rare3d_scarce=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD, PD and FD restrict
get_rare3d_restrict<-as.data.frame(rare3d_restrict %>% group_by(OID_First) %>%
                                   summarise(rare3d_restrict=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD and PD restrict 
get_rare2d_PD_TD_restrict<-as.data.frame(rare2d_PD_TD_restrict %>% group_by(OID_First) %>%
                            summarise(rare2d_PD_TD_restrict=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD and PD restrict 
get_rare2d_PD_TD_scarce<-as.data.frame(rare2d_PD_TD_scarce %>% group_by(OID_First) %>%
                                           summarise(rare2d_PD_TD_scarce=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD and FD restrict
get_rare2d_FD_TD_restrict<-as.data.frame(rare2d_FD_TD_restrict %>% group_by(OID_First) %>%
                                  summarise(rare2d_FD_TD_restrict=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD and FD scarce
get_rare2d_FD_TD_scarce<-as.data.frame(rare2d_FD_TD_scarce %>% group_by(OID_First) %>%
                                           summarise(rare2d_FD_TD_scarce=n_distinct(genus_sp)))

## count the number of rare species in each grid cell
##  TD and FD scarce
get_rare4d_restrict_scarce<-as.data.frame(rare4d_restrict_scarce %>% group_by(OID_First) %>%
                                         summarise(rare4d_restrict_scarce=n_distinct(genus_sp)))



## create joins between the total species and number of rare species 
## per grid cell, convert NA to zeros and add a field for region name

add_rare  <-left_join(getNumSp, get_rare2d_PD_FD, by="OID_First")
add_rare1 <-left_join(add_rare, get_rare3d_restrict, by="OID_First")
add_rare2 <-left_join(add_rare1, get_rare3d_scarce, by="OID_First")
add_rare3 <-left_join(add_rare2, get_rare2d_FD_TD_restrict, by="OID_First")
add_rare4 <-left_join(add_rare3, get_rare2d_FD_TD_scarce, by="OID_First")
add_rare5 <-left_join(add_rare4, get_rare4d_restrict_scarce, by="OID_First")
add_rare6 <-left_join(add_rare5, get_rare2d_PD_TD_restrict, by="OID_First")
add_rare7 <-left_join(add_rare6, get_rare2d_PD_TD_scarce, by="OID_First")
add_rare7[is.na(add_rare7)]<-0

## now load the file with the fishing effor per grid cell

Fisingeffort_subset_ID_f2 <- read.csv("Fisingeffort_subset_ID_f2.csv", header=TRUE)
Fisingeffort_subset_ID_f2 <- Fisingeffort_subset_ID_f2[,-1]

## now summarise fishing effort using group by the grid cells (ID_f)
Fisingeffort_subset_ID_f3<-as.data.frame(Fisingeffort_subset_ID_f2 %>% 
                                           group_by(ID_f) %>%
                            summarise(
                              NomActive= sum(NomActive),
                              EffActive= sum(EffActive),
                              NV= sum(NV),
                              NV= sum(NV),
                              P= sum(P),
                              GT= sum(GT),
                              NomActiveHours= sum(NomActiveHours),
                              EffActiveHours= sum(EffActiveHours),
                              xcoord= unique(XCord),
                              ycoord= unique(YCord)
                              ))

Fisingeffort_subset_ID_f3$ID_f<- as.factor(Fisingeffort_subset_ID_f3$ID_f)

ID_f_OID_First <- unique(ID_f_OID_First)

MAESTRO_rarity_and_fishing_effort <-left_join(add_rare7, ID_f_OID_First, by="OID_First")

str(MAESTRO_rarity_and_fishing_effort)

MAESTRO_rarity_and_fishing_effort$ID_f<- as.factor(MAESTRO_rarity_and_fishing_effort$ID_f)

MAESTRO_rarity_and_fishing_effort1 <-left_join(MAESTRO_rarity_and_fishing_effort, Fisingeffort_subset_ID_f3, by="ID_f")

str(MAESTRO_rarity_and_fishing_effort)

MAESTRO_rarity_and_fishing_effort1<- data.frame(MAESTRO_rarity_and_fishing_effort1)

write.csv(MAESTRO_rarity_and_fishing_effort1, "MAESTRO_rarity_and_fishing_effort_new666.csv")
