
###### G - build a 3D plot with taxonomic, functional and phylogenetic indices ######
#saveRDS(subtree, "subtree_MAESTRO.rds")
#subtree_MAESTRO<- readRDS("subtree_MAESTRO.rds")

subtree_MAESTRO666 <- readRDS("subtree_MAESTRO6.RDS")
phy_tree <- subtree_MAESTRO666 #### 
phy_tree$tip.label <- str_replace_all(phy_tree$tip.label, c(" " = "_"))

MAESTRO_species_indices1<- read.csv("MAESTRO_species_indices1.csv", h= T)
str(MAESTRO_species_indices1)
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

p2 <- gheatmap(p, MAESTRO_species_indices1  [,c("funct_dist_scaled","evol_dist_scaled","scarcity_logbio_scaled")], width=0.2, hjust='left', colnames_angle=-10, font.size=3)  +
   new_scale_fill() +
  
  scale_fill_gradient2(low="blue", high="red", guide="colorbar") +
theme_tree()     
#scale_fill_gradient2 + # if we have negative values in the variable, use scale_fill_gradient2 
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
