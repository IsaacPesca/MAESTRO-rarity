
# Code Used to get the subset of Fishing Effort grid cells
# Part of this process was done using ArcGIS, see comments below
#  Arnaud code adapted by Isaac Trindade Santos
# MAESTRO workgroup

library(dplyr)
if(!require(metR)){install.packages("metR"); library(metR)}
if(!require(patchwork)){install.packages("patchwork"); library(patchwork)}
if(!require(tidyverse)){install.packages("tidyverse"); library(tidyverse)}
if(!require(wesanderson)){install.packages("wesanderson"); library(wesanderson)}
if(!require(metR)){install.packages("metR"); library(metR)}
if(!require(gridExtra)){install.packages("gridExtra"); library(gridExtra)}
if(!require(ggplot2)){install.packages("ggplot2"); library(ggplot2)}
if(!require(mapdata)){install.packages("mapdata"); library(mapdata)}
if(!require(sf)){install.packages("sf"); library(sf)}
if(!require(lwgeom)){install.packages("lwgeom"); library(lwgeom)}
if(!require(grDevices)){install.packages("grDevices"); library(grDevices)}
if(!require(scales)){install.packages("scales"); library(scales)}
if(!require(marmap)){install.packages("marmap"); library(marmap)}
if(!require(data.table)){install.packages("data.table"); library(data.table)}
if(!require(readr)){install.packages("readr"); library(readr)}
if(!require(rnaturalearth)){install.packages("rnaturalearth"); library(rnaturalearth)}
if(!require(reshape2)){install.packages("reshape2"); library(reshape2)}
if(!require(mapplots)){install.packages("mapplots"); library(mapplots)}
if(!require(RColorBrewer)){install.packages("RColorBrewer"); library(RColorBrewer)}
if(!require(classInt)){install.packages("classInt"); library(classInt)}
if(!require(maps)){install.packages("maps"); library(maps)}
if(!require(maptools)){install.packages("maptools"); library(maptools)}
if(!require(viridis)){install.packages("viridis"); library(viridis)}
if(!require(stats)){install.packages("stats"); library(stats)}


# upload fishing effort data (Rousseau et al., 2024)
fishing <- system.time(Fishingefforts <- fread("GriddedEffortby_FishingCountry.csv", showProgress = TRUE))
head(Fishingefforts)
dim(Fishingefforts)
# subset only in NeA and Med areas
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lat>30),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lat<65),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lon>-20),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lon<40),]
dim(Fishingefforts)

# subset only the period of MAESTRO analyses : 1983 to 2019 (no fishing data after 2017, so we keep the period 1983 to 2017)
Fishingefforts<-Fishingefforts[which(Fishingefforts$Year>1982),]
# Fishingefforts$cell1<- paste(Fishingefforts$Lon, Fishingefforts$Lat, sep="_")

Fishingefforts$ID_f<- rownames(Fishingefforts)

Fishingefforts_lat_long<- Fishingefforts %>% dplyr::select(ID_f, Lat, Lon)

write.csv(Fishingefforts_lat_long, "Fishingefforts_lat_long12025.csv")
head(Fishingefforts_lat_long)

# load the subset of fishing effort
# the subset was done in ArcGIS
Fisingeffort_subset_ID_f <- read.csv("Fisingeffort_subset_ID_f1.txt", header=TRUE, sep = ',')

fish_col<- c("XCord", "YCord", "ID_f")
colnames(Fisingeffort_subset_ID_f)<-  fish_col
Fishingefforts6<- data.frame(Fishingefforts) 
Fishingefforts6$ID_f<- as.factor(Fishingefforts6$ID_f)

Fisingeffort_subset_ID_f1 <- unique(Fisingeffort_subset_ID_f )

Fisingeffort_subset_ID_f1$ID_f<- as.factor(Fisingeffort_subset_ID_f1$ID_f)

Fisingeffort_subset_ID_f2 <-left_join( Fisingeffort_subset_ID_f1, 
                                       Fishingefforts6, by="ID_f")

write.csv(Fisingeffort_subset_ID_f2, "Fisingeffort_subset_ID_f2.csv")















