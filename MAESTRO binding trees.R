## MAESTRO
## Phylogenetic tree for MAESTRO research question
## Bind trees - Jawless fish, cartilaginous, bony fish
## Subset MAESTRO species tree

library(phytools)
library(tidyverse)
library(gt)
library(glue)
library(paleotree)
library(dispRity)
#if (!requireNamespace("BiocManager", quietly = TRUE))
#  install.packages("BiocManager")
# BiocManager::install("ggtree")
library(ggtree) # use lines bove if this package here is not installed on your library
library(ape)
library(castor)
library(rotl)
library(gridExtra)

#___________________FIRST PART___________________________________
####### load JAWLESS fish phylogeny

jawless_fish1<-read.tree("jawlessFishTree.tre")
class(jawless_fish1)<- "phylo"
library(RRphylo)
jawless_fish1<-fix.poly(jawless_fish1,type="resolve")
jawless_fish1$Nnode
#[1] 16
jawless_fish1$node.label<- c(1:16)

joutgroup<-c("Myxine_glutinosa", "Petromyzon_marinus","Geotria_australis", "Lampetra_fluviatilis", "Ichthyomyzon_greeleyi", "Ichthyomyzon_unicuspis" )

rooted_tree = root(jawless_fish1, outgroup = joutgroup, resolve.root = TRUE)
is.rooted(jawless_fish1)

class(jawless_fish1)

ggtree(jawless_fish1, layout='circular') + theme_tree2() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

jawless_fish1<- multi2di(jawless_fish1)
max(branching.times(jawless_fish1))
plot(jawless_fish1)
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

ggtree(shark_consensus_MAESTRO, layout='circular') + theme_tree2() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

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
# fish tree of life provides two options
# first option - a tree with 12 thousand species mostly marine
# I tried to use this tree for the MAESTRO project...but there was too many missing species at the tree
# second option - 100 trees with 31 thousand species marine and freshwater
# I decided to use this one because it had a larger match in species from MAESTRO and the tree
# those are the 100 trees with 31 thousand species
# later on we can run some sensitivity analysis with this file
fish_tree_file_MAESTRO32 <- readRDS("treefish32_ALL_updated_corrected.RDS") 
# lets pick one by random
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

#sharks_bonyfish_MAESTRO<-fix.poly(sharks_bonyfish_MAESTRO,type="resolve")

ggtree(sharks_bonyfish_MAESTRO, layout='circular') + theme_tree2()


sharks_bonyfish_MAESTRO<- multi2di(sharks_bonyfish_MAESTRO)
#sharks_bonyfish_MAESTRO <- phytools::force.ultrametric(sharks_bonyfish_MAESTRO)
# Error in dist.nodes(x) : tree too big
is.ultrametric(sharks_bonyfish_MAESTRO)
# [1] FALSE

#MAESTRO1 <-ggtree(sharks_bonyfish_MAESTRO, color="grey", size=0.001,) +
#  theme_tree2() +
#  ggtitle("Bony Fish and Sharks Trees")

#grid.arrange(bf2,shr2,MAESTRO1, ncol=3)

#plot.phylo(sharks_bonyfish_MAESTRO ,cex=0.15,edge.width=0.1, show.tip.label=F)
species<- sharks_bonyfish_MAESTRO$tip.label
length(species)
#[1] 32235
#write.csv(species, "fish_shark_species.csv")
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


#ggtree(sharks_bonyfish_MAESTRO, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Bony and cartilaginous fish")


#ggtree(jawlessfish_sharks_bonyfish_MAESTRO, color="grey", size=0.1,) + 
#  theme_tree2() + 
#  ggtitle("Jawless, Bony and Cartilaginous fish")
# save big tree if needed
#saveRDS(sharks_bonyfish_MAESTRO,"sharks_bonyfish_MAESTRO.RDS")

#___________________FIFTH PART___________________________________
################
## load species list from MAESTRO
# I corrected the species names following FishBase newest version

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
#########################################################
#dev.off()
#########################################################
# shark tree
length(shark_consensus_MAESTRO$tip.label)
#[1] 1154
shr2<-ggtree(shark_consensus_MAESTRO, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label="Sharks - consensus tree", subtitle = "1154 Species")

#########################################################
# jawless fish
length(jawless_fish1$tip.label)
#[1] 1154
jaw2<-ggtree(jawless_fish1, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label="Jawless Fish", subtitle = "17 Species")
#########################################################
# bony fish tree
length(species32$tip.label)
#[1] 31081
bf2<-ggtree(species32, color="grey", size=0.001,) +
  theme_tree2() + 
  ggtitle(label = "Bony fish - fish tree of life", subtitle = "31081 Species")
#########################################################
# binded tree - bony + shark
length(jawlessfish_sharks_bonyfish_MAESTRO$tip.label)
# [1] 32235
MAESTRO1<- ggtree(jawlessfish_sharks_bonyfish_MAESTRO, color="grey", size=0.1,) + 
  theme_tree2() + 
  ggtitle(label = "Jawless, Bony and Shark Fish Tree", subtitle = "32252 Species")
#########################################################
#MAESTRO species tree
length(subtree$tip.label) 
# [1] 609
#MAESTRO all species tree
length(unique(biomass_maestro$genus_sp))
#  [1] 620

MAESTRO_tree_species<- subtree$tip.label
ALL_species<- unique(biomass_maestro$genus_sp)

species_not_in_the_tree<- MAESTRO_tree_species[!ALL_species %in% MAESTRO_tree_species ]

# list of species from MAESTRO that are missing in the tree

# [1] "Scymnodon_ringens"          "Chauliodus_sloani"          "Bathysaurus_ferox"         
# [4] "Hoplostethus_mediterraneus" "Trachyscorpia_cristulata"   "Chelidonichthys_lastoviza" 
# [7] "Eutrigla_gurnardus"         "Chelidonichthys_lucerna"    "Pachycara_crassiceps"      
# [10] "Melanostigma_atlanticum"    "Myoxocephalus_quadricornis" "Liparis_montagui"          
# [13] "Crystallogobius_linearis"   "Lesueurigobius_suerii"      "Gobius_auratus"            
# [16] "Remora_brachyptera" 

MAESTRO2 <-ggtree(subtree, color="grey", size=0.001, ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "604 Species")

ggtree(jawless_fish1, layout='circular') + theme_tree2() + geom_text2(aes(subset=!isTip, label=node), hjust=-.3)

MAESTRO_circular <-ggtree(subtree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "604 Species")



#########################################################
grid.arrange(jaw2,bf2,shr2,MAESTRO1,MAESTRO2, ncol=5)
#########################################################
## save final tree for MAESTRO project
setwd("C:/Users/isaac-trindade/Desktop/01 - OIST - POSTDOC/9 - MAESTRO files")
subtree$node.label<-c(1:5264)
is.binary(subtree)
is.ultrametric(subtree)
subtree$Nnode
# [1] 5264

saveRDS(subtree,"MAESTRO_subtree_5265_species_Jan2024.RDS")
