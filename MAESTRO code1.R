
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

library(psych)
library(phytools)
library(tidyverse)
library(gt)
library(glue)
library(paleotree)
library(dispRity)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree) # use the lines above if this package ggtree below is not installed on your library
library(ape)
library(castor)
library(rotl)
library(gridExtra)
library(rfishbase)
library(dplyr)
library(tidyr)
library(missForest)

######## A - build a tree for MAESTRO species  #########


#___________________FIRST PART___________________________________
####### load JAWLESS fish phylogeny

jawless_fish1<-read.tree("jawlessFishTree.tre")
class(jawless_fish1)<- "phylo"
library(RRphylo)
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


sharks_bonyfish_MAESTRO<- multi2di(sharks_bonyfish_MAESTRO)
is.ultrametric(sharks_bonyfish_MAESTRO)

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


#___________________FIFTH PART___________________________________
####
## load species list from MAESTRO 
####

biomass_maestro <- read.csv("datras_medits_biomass_withtraits.csv", dec = ",", sep = ";")
length(unique(biomass_maestro$genus_sp))
# 620 species
species_maestro<- unique(biomass_maestro$genus_sp)
species_maestro<- as.character(species_maestro)

## subset our MAESTRO tree here

subtree <- castor::get_subtree_with_tips(jawlessfish_sharks_bonyfish_MAESTRO, only_tips=species_maestro)$subtree

#ggtree(subtree, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Jawless, Bony and Cartilaginous fish")

is.binary(subtree)
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
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "604 Species")


MAESTRO_circular <-ggtree(subtree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "604 Species")

####
grid.arrange(jaw2,bf2,shr2,MAESTRO1,MAESTRO2, ncol=5)
####
## save final tree for MAESTRO project
subtree$node.label<-c(1:604)
is.binary(subtree)
is.ultrametric(subtree)
subtree$Nnode
saveRDS(subtree,"MAESTRO_subtree_604_species_Jun2024.RDS")

MAESTRO_circular <-ggtree(subtree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "604 Species")

MAESTRO_circular

MAESTRO_tree_species<- subtree$tip.label
ALL_species<- unique(biomass_maestro$genus_sp)

species_not_in_the_tree<- MAESTRO_tree_species[!ALL_species %in% MAESTRO_tree_species ]
species_not_in_the_tree
# list of species from MAESTRO that are missing in the tree

# [1] "Scymnodon_ringens"          "Chauliodus_sloani"          "Bathysaurus_ferox"         
# [4] "Hoplostethus_mediterraneus" "Trachyscorpia_cristulata"   "Chelidonichthys_lastoviza" 
# [7] "Eutrigla_gurnardus"         "Chelidonichthys_lucerna"    "Pachycara_crassiceps"      
# [10] "Melanostigma_atlanticum"    "Myoxocephalus_quadricornis" "Liparis_montagui"          
# [13] "Crystallogobius_linearis"   "Lesueurigobius_suerii"      "Gobius_auratus"            
# [16] "Remora_brachyptera" 

######## B - calculate evolutionary distinctiveness ########
## phylogenetic dis
phy_tree <- subtree
phy_tree$tip.label <- str_replace_all(phy_tree$tip.label, c("_" = " "))

# evolutionary distinctiveness
phy_tree.cm <- clade.matrix(phy_tree)
Evol_Di <- ed.calc(phy_tree.cm)
head(Evol_Di$spp)
most_distinct_species<-Evol_Di$spp[order(Evol_Di$spp$ED,decreasing=T),];head(most_distinct_species)
scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
most_distinct_species$dist_scaled <- scale_values(most_distinct_species$ED)

write.csv(most_distinct_species, "most_distinct_species.csv")

###### C - prepare trait matrix ######


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

traits_MAESTRO1<-traits_MAESTRO1[-576,]
traits_MAESTRO1<-traits_MAESTRO1[-179,]
traits_MAESTRO1<-traits_MAESTRO1[-180,]
traits_MAESTRO1<-traits_MAESTRO1[-251,]

row.names(traits_MAESTRO1)<- traits_MAESTRO1$genus_sp

# load all traits from fishbase with imputation
all_species_traits_fishbase_imputed <- read.csv("all_species_traits_fishbase_imputed.csv", header=TRUE, stringsAsFactors=FALSE)
colnames(all_species_traits_fishbase_imputed)[colnames(all_species_traits_fishbase_imputed) == 'X'] <- 'genus_sp'
all_species_traits_fishbase_imputed1<- all_species_traits_fishbase_imputed %>% dplyr::select(genus_sp, SwimMode_fam)

traits_MAESTRO2 <- merge(traits_MAESTRO1, all_species_traits_fishbase_imputed1, by = "genus_sp", all.x = TRUE)
row_names<- traits_MAESTRO2$genus_sp
row_names<- as.matrix(row_names)
row_names<- data.frame(row_names)
traits_MAESTRO2<- traits_MAESTRO2 %>% dplyr::select(-one_of('genus_sp' )) 
traits_MAESTRO2 <- as.data.frame(unclass(traits_MAESTRO2),stringsAsFactors=TRUE)
varClass(traits_MAESTRO2) 

traits_MAESTRO3<- missForest(traits_MAESTRO2, maxiter = 10, ntree = 100, variablewise = TRUE,
                             decreasing = FALSE, verbose = TRUE,
                             mtry = floor(sqrt(ncol(traits_MAESTRO2))), replace = TRUE,
                             classwt = NULL, cutoff = NULL, strata = NULL,
                             sampsize = NULL, nodesize = NULL, maxnodes = NULL,
                             xtrue = NA, parallelize = c('no', 'variables', 'forests'))

#missForest iteration 1 in progress...done!
#  estimated error(s): 1134.435 0 0 0 9.540333e+15 8438.098 0.08169935 0.4304207 0.3382353 0.4039735 
#difference(s): 0.02540801 0.003629032 
#time: 2.06 seconds

#missForest iteration 2 in progress...done!
#  estimated error(s): 1059.262 0 0 0 1.010004e+16 8912.219 0.07189542 0.4012945 0.3235294 0.3625828 
#difference(s): 0.005920514 0.0004032258 
#time: 2.04 seconds

#missForest iteration 3 in progress...done!
#  estimated error(s): 1102.458 0 0 0 9.601647e+15 8837.728 0.08006536 0.433657 0.3398693 0.3675497 
#difference(s): 0.005058489 0.0008064516 
#time: 1.89 seconds

#missForest iteration 4 in progress...done!
#  estimated error(s): 1058.547 0 0 0 9.92128e+15 8686.656 0.08006536 0.4061489 0.3267974 0.3692053 
#difference(s): 0.006061819 0.0008064516 
#time: 1.97 seconds
traits_MAESTRO666_error_estimates<- traits_MAESTRO3$OOBerror
traits_MAESTRO33<- traits_MAESTRO3$ximp
traits_MAESTRO33$genus_sp<- row_names$row_names

traits_MAESTRO666<- traits_MAESTRO3$ximp
row.names(traits_MAESTRO666)<- row_names$row_names
write.csv(traits_MAESTRO666, "traits_MAESTRO666.csv")

###### D - calculate functional distinctiveness ######
library(funrar)
vspec$latname <- sub(" ", "_", vspec$latname)
vspec_sp<- vspec  %>% dplyr::select(latname )
vspec_sp<- unique(vspec_sp)
colnames(vspec_sp)[colnames(vspec_sp) == 'latname'] <- 'genus_sp'

MAESTRO_species_indices <- merge(vspec_sp, traits_MAESTRO666, by = "genus_sp", all.x = TRUE)




species_vulnerability<- unique(vspec$latname)

traits_MAESTRO_vulnerability <-traits_MAESTRO666[rownames(species_vulnerability),]

traitT<-compute_dist_matrix(traits_MAESTRO_vulnerability, metric="gower", center=FALSE, scale=FALSE)

functional_distinctiveness<-distinctiveness_global(traitT, di_name="global_di")
colnames(functional_distinctiveness)[colnames(functional_distinctiveness) == 'species'] <- 'genus_sp'
MAESTRO_species_indices<- most_distinct_species
colnames(MAESTRO_species_indices)[colnames(MAESTRO_species_indices) == 'dist_scaled'] <- 'evol_dist_scaled'
colnames(MAESTRO_species_indices)[colnames(MAESTRO_species_indices) == 'species'] <- 'genus_sp'

MAESTRO_species_indices$genus_sp <- chartr(" ", "_", MAESTRO_species_indices$genus_sp )
MAESTRO_species_indices <- merge(MAESTRO_species_indices, functional_distinctiveness, by = "genus_sp", all.x = TRUE)

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_species_indices$funct_dist_scaled <- scale_values(MAESTRO_species_indices$global_di)

###### E - prepare biomass matrix #########
biomass_matrix<- biomass_maestro %>% dplyr::select(genus_sp, Year, biomass)

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
most_distinct_species$species <- chartr(" ", "_", most_distinct_species$species )
row_names1<- most_distinct_species$species

traits_MAESTRO6666 <- traits_MAESTRO666[rownames(traits_MAESTRO666) %in% row_names1, ]
traits_distance_matrix2<-compute_dist_matrix(traits_MAESTRO6666, metric="gower", center=FALSE, scale=FALSE)

funrar_MAESTRO<- funrar(biomass_matrix5, traits_distance_matrix2, rel_abund = TRUE)

scarcity_MAESTRO<- funrar_MAESTRO$Si

library(tidyr)

scarcity_MAESTRO1<- colMeans(scarcity_MAESTRO, na.rm=TRUE)
scarcity_MAESTRO1<- as.matrix(scarcity_MAESTRO1)
scarcity_MAESTRO1<- as.data.frame(scarcity_MAESTRO1)
colnames(scarcity_MAESTRO1)[colnames(scarcity_MAESTRO1) == 'V1'] <- 'scarcity_average'

scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_species_indices$scarcity_average <- scarcity_MAESTRO1$scarcity_average
MAESTRO_species_indices$scarcity_scaled <- scale_values(MAESTRO_species_indices$scarcity_average)


MAESTRO_species_indices1<- MAESTRO_species_indices %>% dplyr::select(genus_sp,evol_dist_scaled, funct_dist_scaled,scarcity_scaled )
row.names(MAESTRO_species_indices1)<- MAESTRO_species_indices1$genus_sp
MAESTRO_species_indices1<- MAESTRO_species_indices1[,-1]
#create pairs plot
pairs.panels(MAESTRO_species_indices1, main= "Biodiversity Facet`s Correlations", stars = TRUE)
library(ppcor)
pcor(MAESTRO_species_indices1)

###### G - build a 3D plot with taxonomic, functional and phylogenetic indices ######

phy_tree <- subtree #### mudei aqui
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

p2 <- gheatmap(p, MAESTRO_species_indices3[,c(1:3)], width=0.2, hjust='left', colnames_angle=-10, font.size=3)  +
  # new_scale_fill() +
  #scale_fill_gradient2 # if we have negative values in the variable, use scale_fill_gradient2 
  scale_fill_gradient2(low="blue", high="red", guide="colorbar") 
theme_tree()     

open_tree(p2, 20) %>% rotate_tree(70)                


library(rgl)
plot3d(MAESTRO_species_indices1$scarcity_scaled,MAESTRO_species_indices1$funct_dist_scaled,MAESTRO_species_indices1$evol_dist_scaled, xlab="Scarcity", ylab="Functional Distinctiveness", zlab="Evolutionary Distinctiveness", pch=5 )
#library(ggplot2)
#library(plotly)
#ggplot(MAESTRO_species_indices1, aes(x=scarcity_scaled, y=funct_dist_scaled, z=evol_dist_scaled)) +   #, color=Species
#  theme_void() 
# rgl::axes_3D() +
#  stat_3D()



library(plotly)

plot_3d <- MAESTRO_species_indices1 %>% plot_ly(x=~scarcity_scaled,
                                                y=~funct_dist_scaled,
                                                z=~evol_dist_scaled, alpha = 0.5)


plot_3d




q.evol<- quantile(MAESTRO_species_indices1$evol_dist_scaled, c(0.25, 0.5, 0.75), type=1)
#25%        50%        75% 
#0.08712192 0.12865583 0.18725516 

q.funct<-quantile(MAESTRO_species_indices1$funct_dist_scaled, c(0.25, 0.5, 0.75), type=1)
#25%       50%       75% 
#0.1051994 0.1828711 0.2750250 

q.scarcity<-quantile(MAESTRO_species_indices1$scarcity_scaled, c(0.25, 0.5, 0.75), type=1)
#25%       50%       75% 
#0.9782451 0.9986295 0.9999192
MAESTRO_species_indices2$abundant_scaled<- 1- MAESTRO_species_indices2$scarcity_scaled
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
  mutate(scarcity_scaled_quartile = case_when(scarcity_scaled <= q.scarcity[[1]] ~ "Quartile1",
                                              scarcity_scaled >= q.scarcity[[1]]  & scarcity_scaled <= q.scarcity[[2]] ~ "Quartile2",
                                              scarcity_scaled >= q.scarcity[[2]]  & scarcity_scaled <= q.scarcity[[3]] ~ "Quartile3",
                                              scarcity_scaled >= q.scarcity[[3]] ~ "Quartile4")) 

MAESTRO_species_indices3<- MAESTRO_species_indices2 %>% 
  mutate(RARE_fun_evol = case_when(funct_dist_scaled <= q.funct[[1]] &  evol_dist_scaled <= q.evol[[1]]~ "Common2D",
                                   funct_dist_scaled >= q.funct[[3]] & evol_dist_scaled >= q.evol[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_fun_evol= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_tax_evol = case_when(scarcity_scaled <= q.scarcity[[1]] &  evol_dist_scaled <= q.evol[[1]]~ "Common2D",
                                   scarcity_scaled >= q.scarcity[[3]] & evol_dist_scaled >= q.evol[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_tax_evol= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_tax_funct = case_when(funct_dist_scaled <= q.funct[[1]]  &  scarcity_scaled <= q.scarcity[[1]]~ "Common2D",
                                    funct_dist_scaled >= q.funct[[3]] & scarcity_scaled >= q.scarcity[[3]]~ "Rare2D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_tax_funct= "Mid2D"))

MAESTRO_species_indices3<- MAESTRO_species_indices3 %>% 
  mutate(RARE_3D = case_when(funct_dist_scaled <= q.funct[[1]]  &  scarcity_scaled <= q.scarcity[[1]]&  evol_dist_scaled <= q.evol[[1]]~ "Common3D",
                             funct_dist_scaled >= q.funct[[3]] & scarcity_scaled >= q.scarcity[[3]]&  evol_dist_scaled >= q.evol[[3]]~ "Rare3D"))

MAESTRO_species_indices3 <- MAESTRO_species_indices3 %>% replace_na(list(RARE_3D= "Mid3D"))

ggplot(MAESTRO_species_indices3, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=RARE_fun_evol)) + geom_point(shape=1, alpha= 0.3 )+ theme_bw() +
  geom_point(aes(size=abundant_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]])+
  geom_hline(yintercept=q.evol[[3]])+ 
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")

ggplot(MAESTRO_species_indices3, aes(x=funct_dist_scaled, y=evol_dist_scaled, color=RARE_3D)) + geom_point(shape=1, alpha= 0.3 )+ theme_bw() +
  geom_point(aes(size=abundant_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.funct[[3]])+
  geom_hline(yintercept=q.evol[[3]])+ 
  xlab("Functional Distinctiveness")+
  ylab("Evolutinary Distinctiveness")


ggplot(MAESTRO_species_indices3, aes(x=scarcity_scaled, y=evol_dist_scaled, color=RARE_tax_evol)) + geom_point(shape=1, alpha= 0.3 )+ theme_bw() +
  geom_point(aes(size=funct_dist_scaled))+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.scarcity[[3]])+
  geom_hline(yintercept=q.evol[[3]])+ 
  xlab("Taxonomic Rarity")+
  ylab("Evolutinary Distinctiveness")

ggplot(MAESTRO_species_indices3, aes(x=scarcity_scaled, y=funct_dist_scaled, color=RARE_fun_evol)) + geom_point(shape=1, alpha= 0.03 )+ theme_bw() +
  geom_point(aes(size=evol_dist_scaled), alpha= 0.5)+
  #  geom_smooth(method=lm)+
  geom_vline(xintercept=q.scarcity[[3]])+
  geom_hline(yintercept=q.funct[[3]])+ 
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
                                                                      evol_dist_scaled,scarcity_scaled, funct_dist_scaled
)
MAESTRO_species_indices26 <- merge(MAESTRO_species_indices4,vspec_MAESTRO26, by = "genus_sp" )
MAESTRO_species_indices85 <- merge(MAESTRO_species_indices4,vspec_MAESTRO85, by = "genus_sp" )

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


library(plotly)

plot_3d_S.TSMr <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_scaled,
                                                         y=~funct_dist_scaled,
                                                         z=~evol_dist_scaled, alpha = 0.5, color =~S.TSMr  )


plot_3d_S.TSMr

plot_3d_AC.hrange <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_scaled,
                                                            y=~funct_dist_scaled,
                                                            z=~evol_dist_scaled, alpha = 0.5, color =~AC.hrange  )


plot_3d_AC.hrange


plot_3d_sens <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_scaled,
                                                       y=~funct_dist_scaled,
                                                       z=~evol_dist_scaled, alpha = 0.5, color =~sens  )


plot_3d_sens


plot_3d_expo <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_scaled,
                                                       y=~funct_dist_scaled,
                                                       z=~evol_dist_scaled, alpha = 0.5, color =~expo  )


plot_3d_expo

plot_3d_vuln <- MAESTRO_species_indices226 %>% plot_ly(x=~scarcity_scaled,
                                                       y=~funct_dist_scaled,
                                                       z=~evol_dist_scaled, alpha = 0.5, color =~vuln  )


plot_3d_vuln

# select the coordinates from the species that we are working on  
species_coordinates_MAESTRO<- biomass_maestro %>% dplyr::select(genus_sp,  Lon,Lat)
write.csv(species_coordinates_MAESTRO, "species_coordinates_MAESTRO.csv")
write.csv(MAESTRO_species_indices3, "MAESTRO_species_indices3.csv")

MAESTRO_species_indices3$genus_sp<- rownames(MAESTRO_species_indices3)
