#for this analysis we will be using the dataframe 
#created in analyzing_morpho_data_cleaned: fixed_field_data_processed
#so to run this file you should first run "analyzing_morpho_data_cleaned"

#### Loading libraries and relevant data####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the Ripley's K function: Kest


fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv")

#### Creating fixed_field_data_processed dataframes for each population ####

LM_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "SD")

#### Importing Shapefile #### 

#turn the BCS polygon into a shapefile and visualize its outline
BCS_polygon <- st_as_sf(BCS_polygon)
plot(BCS_polygon$geometry)

fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                       coords = c("long", "lat"), crs = 4326)

fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

#### Plot the Baja Polygons and Creating Shapefiles for Each Populations ####

#finding minimum and maximum lat and long values
min_all_locality_long <- min(fixed_field_data_processed$long)*1.0002
max_all_locality_long <- max(fixed_field_data_processed$long) - (max(fixed_field_data_processed$long) *.0002)
min_all_locality_lat <- min(fixed_field_data_processed$lat)*1.02
max_all_locality_lat <- max(fixed_field_data_processed$lat) - (max(fixed_field_data_processed$lat)*.02)

#plotting the BCS polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf, aes(color = Locality)) + 
  coord_sf(xlim = c(min_all_locality_long, max_all_locality_long), 
           ylim = c(min_all_locality_lat, max_all_locality_lat))+
  theme_classic()

#finding minimum and maximum lat and long values for lM
LM_min_all_locality_long <- min(LM_fixed_field_data_processed$long)#*1.0002
LM_max_all_locality_long <- max(LM_fixed_field_data_processed$long)# - (max(LM_fixed_field_data_processed$long) *.0002)
LM_min_all_locality_lat <- min(LM_fixed_field_data_processed$lat)#*1.002
LM_max_all_locality_lat <- max(LM_fixed_field_data_processed$lat) #- (max(LM_fixed_field_data_processed$lat)*.002)

#plotting the BCS LM polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf, aes(color = Locality)) + 
  coord_sf(xlim = c(LM_min_all_locality_long, LM_max_all_locality_long), 
           ylim = c(LM_min_all_locality_lat, LM_max_all_locality_lat))+
  theme_classic()

#creating LM boundary shapefile
LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sfc()

LM_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_bbox %>%
  st_as_sfc()

#finding minimum and maximum lat and long values for LC
LC_min_all_locality_long <- min(LC_fixed_field_data_processed$long)#*1.0002
LC_max_all_locality_long <- max(LC_fixed_field_data_processed$long)# - (max(LM_fixed_field_data_processed$long) *.0002)
LC_min_all_locality_lat <- min(LC_fixed_field_data_processed$lat)#*1.002
LC_max_all_locality_lat <- max(LC_fixed_field_data_processed$lat) #- (max(LM_fixed_field_data_processed$lat)*.002)

#plotting the BCS LC polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf, aes(color = Locality)) + 
  coord_sf(xlim = c(LC_min_all_locality_long, LC_max_all_locality_long), 
           ylim = c(LC_min_all_locality_lat, LC_max_all_locality_lat))+
  theme_classic()

#creating LC boundary shapefile
LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sfc()

LC_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_bbox %>%
  st_as_sfc()


#finding minimum and maximum lat and long values for SD
SD_min_all_locality_long <- min(SD_fixed_field_data_processed$long)#*1.0002
SD_max_all_locality_long <- max(SD_fixed_field_data_processed$long)# - (max(LM_fixed_field_data_processed$long) *.0002)
SD_min_all_locality_lat <- min(SD_fixed_field_data_processed$lat)#*1.002
SD_max_all_locality_lat <- max(SD_fixed_field_data_processed$lat) #- (max(LM_fixed_field_data_processed$lat)*.002)

#plotting the BCS SD polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf, aes(color = Locality)) + 
  coord_sf(xlim = c(SD_min_all_locality_long, SD_max_all_locality_long), 
           ylim = c(SD_min_all_locality_lat, SD_max_all_locality_lat))+
  theme_classic()

#creating SD boundary shapefile
SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sfc()

SD_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_bbox %>%
  st_as_sfc()



#### Ripley's K ####

#running the Ripley's K analysis for all points
#creating the ppp for the entire extent of points
BCS_ppp <- ppp(x = fixed_field_data_processed$long, y = fixed_field_data_processed$lat, window = W) #creating the poisson point pattern
plot(BCS_ppp)
K_BCS <- Kest(BCS_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_BCS <- Kest(BCS_ppp) #includes the edge corrections
plot(K_BCS, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_BCS, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot

ggplot(data = LM_fixed_field_data_processed_sf)+
  geom_sf()
View(LM_fixed_field_data_processed_sf)


#### Ripley's K for LM ####
LM_win <- as.owin(LM_fixed_field_data_processed_box)
LM_ppp <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win) #creating the poisson point pattern for lm
plot(LM_ppp, pch = 16, cex = 0.5)
LM_k <- Kest(LM_ppp, correction = "Ripley")
plot(LM_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot


#creating the ppp for LC
LC_ppp <- ppp(x = LC_fixed_field_data_processed$long, y = LC_fixed_field_data_processed$lat, window = W) #creating the poisson point pattern for lm
plot(LC_ppp)

#creating the ppp for SD
SD_ppp <- ppp(x = SD_fixed_field_data_processed$long, y = SD_fixed_field_data_processed$lat, window = W) #creating the poisson point pattern for lm
plot(SD_ppp)


#running the Ripley's K analysis for LC
K_LC <- Kest(LC_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_LC <- Kest(LC_ppp) #includes the edge corrections
plot(K_LC, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_LC, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot

#running the Ripley's K analysis for SD
K_SD <- Kest(SD_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_SD <- Kest(SD_ppp) #includes the edge corrections
plot(K_SD, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_SD, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot


#### Ripley's L ####

#Ripley's L
L <- Lest(BCS_ppp, main=NULL)
L <- Lest(BCS_ppp, main=NULL, correction = "Ripley")
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))


#### ANN Analysis (test for clustering/dispersion) ####

ann.p

#### Moran's I ####


