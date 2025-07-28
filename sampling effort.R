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
if(!require(stats)){install.packages("stats"); library(fread)}

# keep only ICES rect of the MAESTRO area
MAESTRO_area<-fread("table_large_with_env_interanualok_24 sept 2024_yearly.csv",h=T,showProgress = TRUE)
dim(MAESTRO_area)
MAESTRO_rect<-cbind(MAESTRO_area$Lon,MAESTRO_area$Lat)
colnames(MAESTRO_rect)<-c("lon","lat")
head(MAESTRO_rect)
MAESTRO_rect<-data.frame(MAESTRO_rect)
rect_MAESTRO<-ices.rect2(MAESTRO_rect$lon,MAESTRO_rect$lat)
MAESTRO_area$rect<-rect_MAESTRO

listsurveys<-levels(as.factor(MAESTRO_area$mySurvey))
outplistrect_per_survey<-vector("list", length(listsurveys))
for (j in 1:length(listsurveys)) {
  subsurvey<-MAESTRO_area[which(MAESTRO_area$mySurvey==listsurveys[j]),]
  outplistrect_per_survey[[j]]<-levels(as.factor(subsurvey$rect))
  names(outplistrect_per_survey)<-listsurveys
}
outplistrect_per_survey

rect_MAESTRO<-levels(as.factor(rect_MAESTRO))


# upload fishing effort data (Rousseau et al., 2024)
fishing <- system.time(Fishingefforts <- fread("GriddedEffortby_FishingCountry.csv", showProgress = TRUE))
head(Fishingefforts)

# subset only in NeA and Med areas
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lat>30),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lat<65),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lon>-20),]
Fishingefforts<-Fishingefforts[which(Fishingefforts$Lon<40),]
dim(Fishingefforts)

# subset only the period of MAESTRO analyses : 1983 to 2019 (no fishing data after 2017, so we keep the period 1983 to 2017)
Fishingefforts<-Fishingefforts[which(Fishingefforts$Year>1982),]

rect_fishingdata<-ices.rect2(Fishingefforts$Lon,Fishingefforts$Lat)
Fishingefforts<-data.frame(Fishingefforts)
Fishingefforts<-cbind(Fishingefforts,rect_fishingdata)
head(Fishingefforts)
Fishingefforts<-Fishingefforts[Fishingefforts$rect_fishingdata %in% rect_MAESTRO,] 
dim(Fishingefforts)
head(Fishingefforts)
final<-Fishingefforts


# simple map
world <- ne_countries(scale = "medium", returnclass = "sf")
head(final)

final_moy <- final %>%
  group_by(rect_fishingdata) %>%
  summarise(
    effort = mean(EffActiveHours, na.rm = TRUE),
    Lon    = mean(Lon, na.rm = TRUE),
    Lat    = mean(Lat, na.rm = TRUE),
    .groups = "drop"
  )
final_moy<-data.frame(final_moy)
summary(final_moy)
head(final_moy)

# we can log if needed
final_moy$effort<-log10(final_moy$effort+1)

hist(final_moy$effort)

# we calculate deciles to better visualize spatial patterns of fishing effort 
deciles <- quantile(final_moy$effort, probs = seq(0, 1, 0.1), na.rm = TRUE)
labels_deciles <- paste0(format(round(deciles[-length(deciles)], 2), nsmall = 2),
                         " – ",format(round(deciles[-1], 2), nsmall = 2))

# Allocate each point to the corresponding decile gp
final_moy <- final_moy %>%
  mutate(effort_decile = cut(effort,
                             breaks = deciles,
                             include.lowest = TRUE,
                             labels = labels_deciles))

# final map
ggplot() +
  geom_point(data = final_moy, aes(x = Lon, y = Lat, fill = effort_decile), size = 3, shape = 21) +
  geom_sf(data = world, color = 'black', fill = 'grey40', size = 0.01) +
  scale_fill_viridis_d(
    name = "Effective fihing effort (KW*hours; log10 transformed)",
    option = "D",
    guide = guide_legend(reverse = TRUE)) +
  theme_classic() +
  coord_sf(xlim = c(-18, 36), ylim = c(33.2, 65), expand = TRUE)