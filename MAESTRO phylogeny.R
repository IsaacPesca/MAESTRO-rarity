
# Code Used to build the Mega Tree used 
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


#  ape::write.tree(sharks_bonyfish_MAESTRO, file='sharks_bonyfish_MAESTRO_tree.txt')



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
