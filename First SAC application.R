###### E - prepare biomass matrix #########
###### run the SAC approach ######

library(dplyr)
library(reshape2)
library(vegan)
library(scales)
library(mgcv)
library(ggplot2)
library(purrr)
library(tibble)
library(viridis)
library(tidyr)


# load biomass novo aqui
MAESTRO_biomass_assemblage<- read.csv("MAESTRO_biomass_assemblage.csv", h=T)

MAESTRO_biomass_assemblage$ID_f<- as.factor(MAESTRO_biomass_assemblage$ID_f)
MAESTRO_biomass_assemblage<- as.data.frame(MAESTRO_biomass_assemblage)

MAESTRO_new_gridcells_025<- read.csv("biomass_2025_new.txt", h=T)
col_names4<- c("Long", "Lat",  "ID_m", "ID_f", "OID_First")

colnames(MAESTRO_new_gridcells_025) <- col_names4

MAESTRO_new_gridcells_025$ID_f<- as.factor(MAESTRO_new_gridcells_025$ID_f) 
MAESTRO_new_gridcells_025$ID_m<- as.factor(MAESTRO_new_gridcells_025$ID_m) 
MAESTRO_new_gridcells_025$OID_First<- as.factor(MAESTRO_new_gridcells_025$OID_First) 
MAESTRO_new_gridcells_025<-MAESTRO_new_gridcells_025[!grepl("NULL", MAESTRO_new_gridcells_025$ID_f),]
MAESTRO_new_gridcells_025<- as.data.frame(MAESTRO_new_gridcells_025)

# double check for duplicated values before proceeding with the merge function
n_occur <- data.frame(table(MAESTRO_biomass_assemblage$ID_m))
n_occur[n_occur$Freq > 1,]
# zero duplicate
n_occur <- data.frame(table(MAESTRO_new_gridcells_025$ID_m))
n_occur[n_occur$Freq > 1,]
# zero duplicate - good to go

MAESTRO_biomass_assemblage123 <- merge(MAESTRO_biomass_assemblage, 
                                       MAESTRO_new_gridcells_025, 
                                       by = "ID_m", 
                                       all.x = TRUE)


####################################
####################################
###### SAMPLING EFFORT ISSUE #######
####################################
####################################

# following Maureaud et al., 2019    DOI: 10.1098/rspb.2019.1189
# following Coulon et al., 2022      DOI: 10.1111/geb.13731

# SELECT THE NECESSARY COLUMNS TO PROCEED
MAESTRO_biomass_assemblage4<- MAESTRO_biomass_assemblage123 %>% 
  dplyr::select(OID_First, 
                ID_m,
                ID_f.x,
                Survey,
                Ship,
                haul_number,
                Rect,
                hauling_duration,
                Year,
                month,
                day,
                Quarter,
                biomass,
                genus_sp,
                Lon,
                Lat.x
  )

# what is the mean of the haul duration?
mean(MAESTRO_biomass_assemblage123$hauling_duration, na.rm=TRUE)
#  [1] 40.29864

# How many grid cells??
length(unique(MAESTRO_biomass_assemblage4$OID_First))
# [1] 3843

# SAC should be computed for each grid cell
# what is the average of number of hauls per grid cell?

hauls_per_grid <- MAESTRO_biomass_assemblage4 %>% 
  group_by(OID_First) %>% 
  dplyr::summarise("Hauls_frequency"  = length(unique(haul_number)))

hauls_per_grid1 <- MAESTRO_biomass_assemblage4 %>% 
  group_by(OID_First, haul_number ) %>% 
  dplyr::summarise(
    "Species_richness" = length(unique(genus_sp)                           
                                
    ))
write.csv(hauls_per_grid1, "hauls_per_grid1.csv")


# Histogram showing the hauls frequency per grid cell
hist(hauls_per_grid$Hauls_frequency)
# Mean number of haul frequency per grid cell
mean(hauls_per_grid$Hauls_frequency)
# 15.7338
# maximum number
max(hauls_per_grid$Hauls_frequency)
#  225
# minimum number
#  1

# hauls without the duration
sum(is.na(MAESTRO_biomass_assemblage4$hauling_duration))
#[1] 697420

# lets subset grid cells that has at least more than 

# do some subset here

# select the grid cells with less than 4 hauls
# this needs to be done because I got errors when calculating SAC
# its necessary to keep only grid cells that had at least 4 hauls

hauls_per_grid_low_haul_numbers <- subset(hauls_per_grid, Hauls_frequency < 5)
hauls_per_grid_low_haul_numbers1<- unique(hauls_per_grid_low_haul_numbers$OID_First)

#hauls_per_grid_low_haul_numbers2 <- append(hauls_per_grid_low_haul_numbers1, as.factor("27298"))
#hauls_per_grid_low_haul_numbers3 <- append(hauls_per_grid_low_haul_numbers2, as.factor("27057"))


hauls_per_grid_low_haul_numbers2<- as.factor(hauls_per_grid_low_haul_numbers1)

# problematic grid cell, only one species found in four different hauls

str(hauls_per_grid_low_haul_numbers2)

str(MAESTRO_biomass_assemblage4$OID_First)
#class(hauls_per_grid_low_haul_numbers5)

# exclude from the matrix the grid cells with less than 5 hauls
# this is done because the SAC approach needs a minimum amount of communities and species
# for example, if the grid cell has less than 5 hauls at least, it will give an error
MAESTRO_biomass_assemblage55 <- MAESTRO_biomass_assemblage4 %>% 
  dplyr::filter(!OID_First %in% hauls_per_grid_low_haul_numbers2)

#MAESTRO_biomass_assemblage55 <- MAESTRO_biomass_assemblage55 %>% 
#  dplyr::filter(!OID_First %in% c("27057"))

#MAESTRO_biomass_assemblage55 <- MAESTRO_biomass_assemblage55 %>% 
#  dplyr::filter(!OID_First %in% c("30318"))

# MAESTRO_biomass_assemblage44 <- MAESTRO_biomass_assemblage4 %>% dplyr::filter(OID_First != hauls_per_grid_low_haul_numbers5)

hauls_per_grid3 <- MAESTRO_biomass_assemblage55 %>% 
  group_by(OID_First,haul_number) %>% 
  dplyr::summarise("Species_richness"  = length(unique(genus_sp)))

more_than_one<- subset(hauls_per_grid3, Species_richness < 2 )

MAESTRO_biomass_assemblage55 <- MAESTRO_biomass_assemblage55 %>% 
  dplyr::filter(!OID_First %in% more_than_one$OID_First)


# save the grid cells IDs in a vector to use in the loop
grid_cell_list <- unique(MAESTRO_biomass_assemblage55$OID_First)
SAC_per_grid_cell_list<- list()
SAC_per_grid_cell_vector<- c()

# build a community matrix using the number of hauls in each grid cells to perform the SAC framework

haul_matrix<- MAESTRO_biomass_assemblage55 %>% dplyr::select(OID_First,genus_sp ,haul_number)


sp1_list <- list()
sp2_random_list <- list()
sp3_rarefaction_list <- list()
mod1_list<- list()
mod2_random_list <- list()
mod3_rarefaction_list <- list()
coef_mod1_list <- list()
fitted_mod1_list <- list()
hauls_and_species_per_grid = data.frame(grid_cell_number = NA, species_richness = NA, haul_numbers= NA)
species_per_haul_grid <- list()

# first, test the loop using just a few grid cells
# loop is working well
# grid_cell_list1 <- c(25329,25061, 26755)
# grid_cell_list1<- as.factor(grid_cell_list1)

for ( grid_index in grid_cell_list) {
  
  # subset by grid cell
  print(grid_index)
  haul_matrix_sub<-subset(haul_matrix, OID_First == grid_index ) #  25329  grid_index
  
  grid_cell_number<- grid_index
  species_richness<- length(unique(haul_matrix_sub$genus_sp))
  haul_numbers<-length(unique(haul_matrix_sub$haul_number))
  hauls_and_species_per_grid1 <- data.frame(grid_cell_number = NA, species_richness = NA, haul_numbers= NA)
  
  hauls_and_species_per_grid1$grid_cell_number <-  grid_cell_number
  hauls_and_species_per_grid1$species_richness <-  species_richness
  hauls_and_species_per_grid1$haul_numbers <-  haul_numbers
  
  haul_matrix_sub1_wide<- dcast(haul_matrix_sub, haul_number~genus_sp, length)
  
  # make sure that it is a matrix
  
  rownames(haul_matrix_sub1_wide)<-haul_matrix_sub1_wide$haul_number
  haul_matrix_sub1_wide<-haul_matrix_sub1_wide[,-1]
  haul_matrix_sub1_wide<- data.frame(haul_matrix_sub1_wide)
  haul_matrix_sub1_wide<- as.matrix(haul_matrix_sub1_wide)
  
  # apply sac with using the community matrix (species per hauls) from each grid cell
  sp1 <- specaccum(haul_matrix_sub1_wide)
  ## sp2_random <- specaccum(haul_matrix_sub1_wide, "random")
  sp3_rarefaction <- specaccum(haul_matrix_sub1_wide, "rarefaction")
  
  ## Fit "michaelis-menten" model 
  
  mod1 <- fitspecaccum(sp1, "michaelis-menten")
  mod1_list[[grid_index]] <- mod1 
  
  ##mod2_random <- fitspecaccum(sp2_random, "michaelis-menten")
  ##mod2_random_list[[ grid_index)]] <- mod2_random 
  
  mod3_rarefaction <- fitspecaccum(sp3_rarefaction, "michaelis-menten")
  mod3_rarefaction_list[[grid_index]] <- mod3_rarefaction 
  
  coef <- coef(mod1)
  coef_mod1_list[[grid_index]] <- coef
  
  fitted <- fitted(mod1)
  fitted_mod1_list[[grid_index]] <- fitted
  
  hauls_and_species_per_grid<- rbind(hauls_and_species_per_grid, hauls_and_species_per_grid1)
  
}

hauls_and_species_per_grid<- hauls_and_species_per_grid[-1,]
saveRDS(hauls_and_species_per_grid, "hauls_and_species_per_grid.RDS")
saveRDS(mod1_list, "mod1_list.RDS")

# # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # # # # # # 

sac_objects<- mod1_list
haul_numbers <- hauls_and_species_per_grid

# # # # #  both lines here
# # # # # # # # # # # # # 

# Load data (replace paths)

# Calculate asymptotic richness and 65% threshold per grid
grid_summary <- data.frame(
  GridID = names(sac_objects),
  AsympRichness = sapply(sac_objects, function(x) max(x$fitted)),
  ActualHauls = haul_numbers$haul_numbers[match(names(sac_objects), haul_numbers$grid_cell_number)]
) |>
  mutate(
    Threshold65 = 0.65 * AsympRichness,
    # Find hauls needed to reach 65% threshold
    HaulsToThreshold = sapply(sac_objects, function(sac) {
      sac_df <- data.frame(Sites = sac$sites, Richness = sac$fitted)
      sac_df |> 
        filter(Richness >= Threshold65[1]) |> 
        slice(1) |> 
        pull(Sites)
    }),
    NeedsMoreHauls = HaulsToThreshold > ActualHauls
  )

# Combine SACs into a single data frame
sac_data <- lapply(seq_along(sac_objects), function(i) {
  grid_id <- names(sac_objects)[i]
  hauls <- haul_numbers$haul_numbers[haul_numbers$grid_cell_number == grid_id]
  
  # Repeat haul numbers to match number of sites
  if (length(hauls) < length(sac_objects[[i]]$sites)) {
    hauls <- rep(hauls, length.out = length(sac_objects[[i]]$sites))
  }
  
  data.frame(
    GridID = grid_id,
    Hauls = hauls,
    Sites = sac_objects[[i]]$sites,
    Richness = sac_objects[[i]]$fitted
  )
}) |> bind_rows()

# Convert Sites to numeric (if not already)
sac_data$Sites <- as.numeric(sac_data$Sites)


# Check structure of sac_data
str(sac_data)

# Check structure of grid_summary
str(grid_summary)

# Convert Sites to numeric in sac_data
sac_data$Sites <- as.numeric(as.character(sac_data$Sites))

# Convert HaulsToThreshold to numeric in grid_summary
grid_summary$HaulsToThreshold <- as.numeric(as.character(grid_summary$HaulsToThreshold))

# Remove any resulting NA values
sac_data <- sac_data[!is.na(sac_data$Sites), ]
grid_summary <- grid_summary[!is.na(grid_summary$HaulsToThreshold), ]



##### funciona aqui

# calculate and smooth the mean SAC curve (but don't plot it)
mean_sac <- sac_data %>%
  group_by(Sites) %>%
  summarize(Richness = mean(Richness, na.rm = TRUE)) %>%
  mutate(SmoothRichness = predict(gam(Richness ~ s(Sites), data = .)))

# calculate asymptotic richness from smoothed mean
asymp_richness <- max(mean_sac$SmoothRichness, na.rm = TRUE)
threshold_65 <- 0.65 * asymp_richness

# find where smoothed mean crosses 65% threshold
haul_to_threshold <- mean_sac %>%
  filter(SmoothRichness >= threshold_65) %>%
  slice(1) %>%
  pull(Sites)

# create the plot without mean line
ggplot() +
  # Original individual curves
  geom_line(
    data = sac_data,
    aes(x = Sites, y = Richness, group = GridID, color = Hauls),
    alpha = 0.4,
    linewidth = 0.6
  ) +
  # Vertical line at hauls needed for 65% threshold (red dashed)
  geom_vline(
    xintercept = haul_to_threshold,
    color = "red",
    linetype = "dashed",
    linewidth = 1.2
  ) +
  # Horizontal line at 65% threshold (black dashed)
  geom_hline(
    yintercept = threshold_65,
    color = "black",
    linetype = "dashed",
    linewidth = 1
  ) +
  # Color scale
  scale_color_gradientn(
    colors = c("blue", "green", "orange", "red"),
    values = rescale(c(1, 75, 150, 225)),
    name = "Hauls per Grid"
  ) +
  # Axis settings
  scale_x_continuous(
    limits = c(0, 225),
    breaks = seq(0, 225, 25),
    expand = c(0, 0)
  ) +
  scale_y_continuous(
    limits = c(0, 142),
    breaks = seq(0, 140, 20),
    expand = c(0, 0)
  ) +
  # llabels and theme
  labs(
    x = "Number of Hauls",
    y = "Species Richness",
    title = "Species Accumulation Curves",
    subtitle = paste("65% threshold reached at", haul_to_threshold, "hauls")
  ) +
  theme_minimal() +
  theme(legend.position = "bottom")

# print the threshold information
message("65% threshold value: ", round(threshold_65, 1))
message("Hauls needed to reach threshold: ", haul_to_threshold)
##### funciona aqui



# final subset here
hauls_and_species_per_grid_KEEP<- subset(hauls_and_species_per_grid, haul_numbers > 18)

hauls_and_species_per_grid_EXCLUDE<- subset(hauls_and_species_per_grid, haul_numbers < 19)
hauls_and_species_per_grid_EXCLUDE1 <- unique(hauls_and_species_per_grid_EXCLUDE$grid_cell_number)


MAESTRO_biomass_assemblage666 <- MAESTRO_biomass_assemblage4 %>% 
  dplyr::filter(!OID_First %in% hauls_and_species_per_grid_EXCLUDE1)

saveRDS(MAESTRO_biomass_assemblage666, "MAESTRO_biomass_assemblage666.RDS")
####################################
####################################
###### SAMPLING EFFORT ISSUE #######
####################################
####################################

