##LIBRARIES

library(sf)
library(ggplot2)
library(viridis)
library(gridExtra)

##BASEMAP
setwd("C:/Users/isaac/Desktop/MAESTRO/MAESTRO")
wcm<-st_read('ne_50m_land.shp')
GRIDS<-st_read('GRID COPY.shp')

MAESTRO_rarity_and_fishing_effort<- read.csv( "MAESTRO_rarity_and_fishing_effort.csv", h=T)
MAESTRO_rarity_and_fishing_effort<-MAESTRO_rarity_and_fishing_effort[-1347,]
MAESTRO_rarity_and_fishing_effort1 <- st_as_sf(MAESTRO_rarity_and_fishing_effort,
                    coords = c("xcoord", "ycoord"),
                    crs = st_crs(GRIDS))

AquaMaps_rarity_and_fishing_effort<- read.csv( "AquaMaps_rarity_and_fishing_effort.csv", h=T)
#AquaMaps_rarity_and_fishing_effort<-AquaMaps_rarity_and_fishing_effort[-1347,]
AquaMaps_rarity_and_fishing_effort1 <- st_as_sf(AquaMaps_rarity_and_fishing_effort,
                                               coords = c("xcoord", "ycoord"),
                                               crs = st_crs(GRIDS))

MAESTRO_rarity_and_fishing_effort<- MAESTRO_rarity_and_fishing_effort[,-1]
AquaMaps_rarity_and_fishing_effort<- AquaMaps_rarity_and_fishing_effort[,-1]
#AquaMaps_rarity_and_fishing_effort<- AquaMaps_rarity_and_fishing_effort[,-1]
AquaMaps_rarity_and_fishing_effort2<- AquaMaps_rarity_and_fishing_effort[,-c(7:15)]

MAESTRO_AquaMaps<- merge(MAESTRO_rarity_and_fishing_effort[,-c(7:15)], AquaMaps_rarity_and_fishing_effort2)
library(corrplot)
library(psych)
library(ggpmisc)

corrplot(cor(MAESTRO_AquaMaps[,-1]))
pairs.panels(MAESTRO_AquaMaps[,-c(1, 4:6,9:11)])
pairs.panels(AquaMaps_rarity_and_fishing_effort[,c(2:3,12:13)])

ggplot(MAESTRO_AquaMaps,aes(rare2dPDFDAqua, rare2d_PD_FD)) +
   
  geom_smooth(method='lm', formula= y~x) +
  xlab("AquaMaps - Numer of Rare Species") +
  ylab("MAESTRO - Numer of Rare Species") +
  ggtitle("Functional and Phylogenetic Rarity")





#ggplot() +
#  geom_sf(data = GRIDS, fill = "grey10") +
#  geom_sf(data = MAESTRO_rarity_and_fishing_effort1, colour = "darkorange", alpha = 0.5, size = 0.25)

MAESTRO_join <- st_join(GRIDS, MAESTRO_rarity_and_fishing_effort1, left = FALSE)
AquaMaps_join <- st_join(GRIDS, AquaMaps_rarity_and_fishing_effort1, left = FALSE)
write.csv(MAESTRO_join, "MAESTRO_gridcell_rarity_fishing_effort.csv")
write.csv(AquaMaps_join, "AquaMaps_gridcell_rarity_fishing_effort.csv")

# ggplot() +
#   geom_sf(data=MAESTRO_join) +
#   geom_sf(data=MAESTRO_join, aes(fill=rare2d_PD_FD), lwd=0.4, alpha=0.4) +
#   theme_bw()+
#   scale_fill_viridis(option = "C", "Rarity PD and FD")+
#   ggtitle("Hotspots of Species Rare Phylogenetically and Functionally")



p <- ggplot(MAESTRO_join) +
  geom_sf(aes(fill = rare2d_PD_FD), linewidth = 0, alpha = 0.9) +
  theme_void() +
  scale_fill_viridis_c(
    trans = "log", breaks = c(1, 2, 4, 8, 12, 16),
    name = "Number of Rare Species",
    guide = guide_legend(
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      label.position = "bottom",
      title.position = "top",
      nrow = 1
    )
  ) +
  labs(
    title = "Rarity Hotspots (MAESTRO)",
    subtitle = "High Phylogenetic and Functional Diversity",
    caption = "Data: INSEE | Creation: Yan Holtz | r-graph-gallery.com"
  ) +
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(
      size = 20, hjust = 0.01, color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.4, l = 2,
        unit = "cm"
      )
    ),
    plot.subtitle = element_text(
      size = 15, hjust = 0.01,
      color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.43, l = 2,
        unit = "cm"
      )
    ),
    plot.caption = element_text(
      size = 10,
      color = "#4e4d47",
      margin = margin(
        b = 0.3, r = -99, t = 0.3,
        unit = "cm"
      )
    ),
    legend.position = c(0.7, 0.5)
  )

p




a <- ggplot(AquaMaps_join) +
  geom_sf(aes(fill = rare2dPDFDAqua), linewidth = 0, alpha = 0.9) +
  theme_void() +
  scale_fill_viridis_c(
    trans = "log", breaks = c(1, 10, 20, 30, 40, 50),
    name = "Number of Rare Species (AquaMaps)",
    guide = guide_legend(
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      label.position = "bottom",
      title.position = "top",
      nrow = 1
    )
  ) +
  labs(
    title = "Rarity Hotspots (AquaMaps)",
    subtitle = "High Phylogenetic and Functional Diversity",
    caption = "Data: INSEE | Creation: Yan Holtz | r-graph-gallery.com"
  ) +
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(
      size = 20, hjust = 0.01, color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.4, l = 2,
        unit = "cm"
      )
    ),
    plot.subtitle = element_text(
      size = 15, hjust = 0.01,
      color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.43, l = 2,
        unit = "cm"
      )
    ),
    plot.caption = element_text(
      size = 10,
      color = "#4e4d47",
      margin = margin(
        b = 0.3, r = -99, t = 0.3,
        unit = "cm"
      )
    ),
    legend.position = c(0.7, 0.5)
  )

a










scale_values <- function(x){(x-min(x))/(max(x)-min(x))}
MAESTRO_join$EffActiveHours_log <- log10(MAESTRO_join$EffActiveHours) +1.53


pp <- ggplot(MAESTRO_join) +
  geom_sf(aes(fill = EffActiveHours), linewidth = 0, alpha = 0.9) +
  theme_void() +
  scale_fill_viridis_c(
    trans = scales::pseudo_log_trans(sigma = 0.001), breaks = c(0,2,3,4,6,8 ),
    name = "Fishing Effort",
    guide = guide_legend(
      keyheight = unit(3, units = "mm"),
      keywidth = unit(12, units = "mm"),
      label.position = "bottom",
      title.position = "top",
      nrow = 1
    )
  ) +
  labs(
    title = "Fishing Effort Hotspots",
    subtitle = "Effective fishing effort, in kW × days at sea",
    caption = "Data: INSEE | Creation: Yan Holtz | r-graph-gallery.com"
  ) +
  theme(
    text = element_text(color = "#22211d"),
    plot.background = element_rect(fill = "#f5f5f2", color = NA),
    panel.background = element_rect(fill = "#f5f5f2", color = NA),
    legend.background = element_rect(fill = "#f5f5f2", color = NA),
    plot.title = element_text(
      size = 20, hjust = 0.01, color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.4, l = 2,
        unit = "cm"
      )
    ),
    plot.subtitle = element_text(
      size = 15, hjust = 0.01,
      color = "#4e4d47",
      margin = margin(
        b = -0.1, t = 0.43, l = 2,
        unit = "cm"
      )
    ),
    plot.caption = element_text(
      size = 10,
      color = "#4e4d47",
      margin = margin(
        b = 0.3, r = -99, t = 0.3,
        unit = "cm"
      )
    ),
    legend.position = c(0.7, 0.5)
  )

pp




grid.arrange(p, a, pp)
