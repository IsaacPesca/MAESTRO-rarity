# MAESTRO phylogenetic tree - plot tree with 3D rarity indices
# Isaac Trindade Santos 

subtree_MAESTRO<- readRDS("subtree_MAESTRO.rds")
MAESTRO_species_indices3d<- read.csv("MAESTRO_species_indices1.csv", h=T )
rownames(MAESTRO_species_indices3d)<- MAESTRO_species_indices3d$X
MAESTRO_species_indices3d<- MAESTRO_species_indices3d[, -1]

######  build a simple circular tree with MAESTRO species tree

MAESTRO_circular <-ggtree(subtree, color="grey", size=0.001,, layout='circular', ignore.negative.edge=TRUE) +
  theme_tree2() +
  ggtitle(label = "MAESTRO - Subset Tree",subtitle = "602 Species")

MAESTRO_circular

###### build a circular tree with taxonomic, functional and phylogenetic indices ######

phy_tree <- subtree_MAESTRO #### mudei aqui
phy_tree$tip.label <- str_replace_all(phy_tree$tip.label, c(" " = "_"))

# plot phylogenetic tree
# Define colors
cols <- rev(jet(10))

# add several variables                

p <- ggtree(phy_tree, layout='circular', color="grey50", size=0.3) +
  annotate('text', x=0, y=40, label='', family='mono', size=16)# +

p2 <- gheatmap(p, MAESTRO_species_indices3d[,c(1:3)], width=0.2, hjust='left', colnames_angle=-10, font.size=3)  +
  scale_fill_gradient2(low="blue", high="red", guide="colorbar") 
theme_tree()     

open_tree(p2, 20) %>% rotate_tree(70)                


