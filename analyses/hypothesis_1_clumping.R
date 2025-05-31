# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei is more clumped or dispersed than by random%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#This script is for investigating whether Quercus brandegeei seeds are predominantly dispersed either by heavy rainfall or gravity, 
#which may impact both the spatial distribution and the genetic structure of a given population. 
  #If they were dispersed by heavy rainfall, we would expect more dispersal than at random. 
  #If they were dispersed by gravity, we would expect more clumping of trees than at random. 
#To test this, we used a Ripley's K and compared the known tree locations to randomly generated locations
#produced in convex hulls, the river shapefiles, and buffers around the river 
#(all three were attempted seperately and compared) to determine whether the 
#tree points were more clustered or dispersed than we would expect at random. 
#We then used an Average Nearest Neighbor Analysis (ANN) to support with p-values whether the points 
#seem more clustered or dispersed than at random. 
#Finally, we used Poisson Point Models, to see if models that take into account the effect that the inverse 
#distance to the river of the points has on the placement of the trees better explains the distribution of 
#the points than if they were distributed at random.

# The Ripley's K and ANN both test whether the trees are more clumped or dispersed than at random, whereby 
#the ANN is the only test of the two that provides a p-value of the two. 
#The PPM test is the only test of the three that sees whether the river itself has an influence on the 
#points distribution.

# The script is broken into sections of 
        #1) loading and processing the packages and spatial data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations, 
        #2) running the Ripley's K analysis, 
        #3) running the Average Nearest Neighber (ANN) analysis,
        #4) running the Poisson Point Model Analysis. 


#### Loading libraries and relevant data ####

library(tidyverse) 
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the Ripley's K function: Kest
library(stars) # for sf_rasterize function
library(raster) #to use crop
library(starsExtra) #to use dist_to_nearest
library(geostatsp) # To successfully use as.im
library(tmaptools)

# loading in the tree data (size, elevation, lat/lon, ID, size/shape)

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

# creating the point shapefiles of the tree locations for each population in UTM 12 N

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#creating a transformed point shapefile with UTM 12 N an equal area projection
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912)

#storing point shapefiles for the trees by population

LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sf()

LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sf()

SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sf()

# Importing BCS and River Shapefiles

#turning the Baja California Sur polygon into a shapefile, to be able to visualize the point locations
BCS_polygon <- read_sf("./data/Shapefiles/BCS_Shapefile/bcs_entidad.shp")
BCS_polygon <- st_as_sf(BCS_polygon) #ensures foreign will be an sf object (in this case a multipolygon)
plot(BCS_polygon$geometry) #plotting the polygon

#Loading in ArcGIS river shapefile and storing out polygons for each population

#Las Matancitas (LM)
river_LM <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LM River/LM_Rivers_Final.shp") 
river_LM  <- river_LM$geometry[1]

#La Cobriza (LC)
river_LC  <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LC River/LC_Rivers_Final.shp")
river_LC  <- river_LC$geometry[1]

#San Dionisio (SD)
river_SD <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/SD River/SD_Rivers_Final.shp")
river_SD <- river_SD$geometry[1]


#transforming the coordinate reference system of the river polygons to be equal area projection (UTM 12N), 
#uses meters as distance measurement
river_LM_trans <- st_as_sf(st_transform(river_LM, crs = 26912))
river_LC_trans <- st_as_sf(st_transform(river_LC, crs = 26912))
river_SD_trans <- st_as_sf(st_transform(river_SD, crs = 26912))


#plotting points with river shapefiles

#LM
ggplot()+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed_sf) +
  theme_light()

#LC
ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_sf) +
  theme_light()

#SD
ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_sf) +
  theme_light()


#creating buffers around the rivers and visualizing the buffers with the river/trees
#buffers 

#LM
river_buffer_LM <- st_buffer(river_LM_trans, 100) #100 m buffer, sf object
ggplot()+
  geom_sf(data = river_buffer_LM)+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed_sf)

#LC
river_buffer_LC<- st_buffer(river_LC_trans, 100) #100 m buffer, sf object
ggplot()+
  geom_sf(data = river_buffer_LC)+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_sf)

#SD
river_buffer_SD <- st_buffer(river_SD_trans, 70) #70 m buffer, sf object
ggplot()+
  geom_sf(data = river_buffer_SD)+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_sf)

## Creating fixed_field_data_processed dataframes for each population ##

LM_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "SD")

#### Ripley's K Analysis (version with box, convex hull, and 20 m buffer around river) ####

## Plotting the Baja/tree Polygons and Creating Bounding Boxes using tree points of each population ##

# all trees

#finding minimum and maximum lat and long values for visualizing the tree locations with *1.002 wiggle room
min_all_locality_long <- min(fixed_field_data_processed$long)*1.0002
max_all_locality_long <- max(fixed_field_data_processed$long) - (max(fixed_field_data_processed$long) *.0002)
min_all_locality_lat <- min(fixed_field_data_processed$lat)*1.02
max_all_locality_lat <- max(fixed_field_data_processed$lat) - (max(fixed_field_data_processed$lat)*.02)

#plotting the BCS polygon with all of the tree points, colored by locality
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf_transformed, aes(color = Locality)) + 
  coord_sf(xlim = c(min_all_locality_long, max_all_locality_long), 
           ylim = c(min_all_locality_lat, max_all_locality_lat))+
  theme_classic()

#creating BCS boundary shapefile, turning sf of all points into sfc
fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  st_as_sfc()

#creating a boundary box with the UTM 12 N min and max lat lon values and then 
#turning it into a simple feature geometry
#will be useful for the Ripley's K
fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  st_bbox %>%
  st_as_sfc()

# LM

#finding minimum and maximum lat and long values for LM
LM_min_all_locality_long <- min(LM_fixed_field_data_processed$long)
LM_max_all_locality_long <- max(LM_fixed_field_data_processed$long)
LM_min_all_locality_lat <- min(LM_fixed_field_data_processed$lat)
LM_max_all_locality_lat <- max(LM_fixed_field_data_processed$lat) 

#plotting the BCS LM polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf_transformed, aes(color = Locality)) + 
  coord_sf(xlim = c(LM_min_all_locality_long, LM_max_all_locality_long), 
           ylim = c(LM_min_all_locality_lat, LM_max_all_locality_lat))+
  theme_classic()

#creating LM boundary shapefile, turning sf of all points into sfc
LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sfc()

#creating a boundary box of LM with the UTM 12 N min and max lat lon values and then 
#turning it into a simple feature geometry
LM_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_bbox %>%
  st_as_sfc()

# plotting the points, river, and bounding box
ggplot(LM_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = river_LM_trans) +
  geom_sf(data = LM_fixed_field_data_processed_sf)

# LC

#finding minimum and maximum lat and long values for LC
LC_min_all_locality_long <- min(LC_fixed_field_data_processed$long)
LC_max_all_locality_long <- max(LC_fixed_field_data_processed$long)
LC_min_all_locality_lat <- min(LC_fixed_field_data_processed$lat)
LC_max_all_locality_lat <- max(LC_fixed_field_data_processed$lat) 

#plotting the BCS LC polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf_transformed, aes(color = Locality)) + 
  coord_sf(xlim = c(LC_min_all_locality_long, LC_max_all_locality_long), 
           ylim = c(LC_min_all_locality_lat, LC_max_all_locality_lat))+
  theme_classic()

#creating LC boundary shapefile, turning sf of all points into sfc
LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sfc()

#creating a boundry box of LC with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
LC_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_bbox %>%
  st_as_sfc()

# plotting the points, river, and bounding box
ggplot(LC_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = river_LC_trans) +
  geom_sf(data = LC_fixed_field_data_processed_sf)

# SD

#finding minimum and maximum lat and long values for SD
SD_min_all_locality_long <- min(SD_fixed_field_data_processed$long)
SD_max_all_locality_long <- max(SD_fixed_field_data_processed$long)
SD_min_all_locality_lat <- min(SD_fixed_field_data_processed$lat)
SD_max_all_locality_lat <- max(SD_fixed_field_data_processed$lat) 

#plotting the BCS SD polygon with the tree points
ggplot(data = BCS_polygon) +
  geom_sf() +
  geom_sf(data = fixed_field_data_processed_sf_transformed, aes(color = Locality)) + 
  coord_sf(xlim = c(SD_min_all_locality_long, SD_max_all_locality_long), 
           ylim = c(SD_min_all_locality_lat, SD_max_all_locality_lat))+
  theme_classic()

#creating SD boundary shapefile, turning sf of all points into sfc
SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sfc()

#creating a boundry box of SD with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
SD_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_bbox %>%
  st_as_sfc()

# plotting the poins, river, and bounding box
ggplot(SD_fixed_field_data_processed_box)+
  geom_sf() +
  geom_sf(data = river_SD_trans) +
  geom_sf(data = SD_fixed_field_data_processed_sf)

## Creating Convex Hulls using tree points of each population ##

# creating one geometry with the tree points using st_union and then creating a convex hull with st_convex_hull

# LM
river_LM_convex_hull <- st_convex_hull(st_union(LM_fixed_field_data_processed_sf))  
ggplot(river_LM_convex_hull)+
  geom_sf() +
  geom_sf(data = river_LM_trans) +
  geom_sf(data = LM_fixed_field_data_processed_sf)

# LC
river_LC_convex_hull <- st_convex_hull(st_union(LC_fixed_field_data_processed_sf)) 
ggplot(river_LC_convex_hull)+
  geom_sf() +
  geom_sf(data = river_LC_trans) +
  geom_sf(data = LC_fixed_field_data_processed_sf)

# SD
river_SD_convex_hull <- st_convex_hull(st_union(SD_fixed_field_data_processed_sf)) 
ggplot(river_SD_convex_hull)+
  geom_sf() +
  geom_sf(data = river_SD_trans) +
  geom_sf(data = SD_fixed_field_data_processed_sf)

## Calculating the Ripley's K ##

# For all Ripley's K:
  # 1) Turn the bounding box into a window object
  # 2) Create the Poisson Point Pattern within the window (a randomized version of the point locations)
  # 3) Calculate the Ripley's K function (Kest), using the Ripley's correction (edge correction)
  # 4) Plot the Ripley's K function output. Dotted red line is the randomized point locations. 
        # r is the distance from a center point (radius)
        # K(r) is the expected number of random points within the radius of a point
        # If the K(r) of our known points (black line) is lower than the randomized points (red line) than the points are more 
              #dispersed than expected at random
        # If the K(r) of our known points (black line) is higher than the randomized points (red line) than the points are more 
              #clustered than expected at random

# All points

win <- as.owin(fixed_field_data_processed_box) #turning the box into a window
ppp <- as.ppp(st_coordinates(fixed_field_data_processed_sf), W = win) #creating the poisson point pattern for LM
plot(ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
K <- Kest(ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) # plotting Ripley's output

# Ripley's K for LM 

# Bounding Box
LM_win <- as.owin(LM_fixed_field_data_processed_box) #turning the box into a window
LM_ppp <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win) #creating the poisson point pattern for LM
plot(LM_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
LM_k <- Kest(LM_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LM_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Convex Hull
LM_win_convex <- as.owin(river_LM_convex_hull)  #turning the convex hull into a window
LM_ppp_convex <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_convex) #creating the poisson point pattern for lm
plot(LM_ppp_convex, pch = 16, cex = 0.5) # plotting the randomized point pattern
LM_k_convex <- Kest(LM_ppp_convex, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LM_k_convex, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Buffer River (100 m)
LM_win_buffer <- as.owin(river_buffer_LM) #turning the buffer into a window
LM_ppp_buffer <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_buffer) #creating the poisson point pattern for lm
plot(LM_ppp_buffer, pch = 16, cex = 0.5) # plotting the randomized point pattern
LM_k_buffer <- Kest(LM_ppp_buffer, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LM_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Making a nicer version of the plot for papers/presentations
plot(LM_k_buffer, main=NULL, las=1, ylab = "", xlab = "", legendargs=list(cex=0.8, xpd=TRUE), 
     yaxp = c(0, 220000, 10))#legend inside of the plot
par(mar = c(6,6,6,6))
title(ylab = bquote(italic("K(r), Las Matancitas trees")), cex = 1.2, line = 4)
title(xlab = bquote(italic("r (m)")),)


#Ripley's K for LC 

# Bounding Box
LC_win <- as.owin(LC_fixed_field_data_processed_box) #turning the box into a window
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
LC_k <- Kest(LC_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Convex Hull
LC_win_convex <- as.owin(river_LC_convex_hull) #turning the convex hull into a window
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_convex) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
LC_k <- Kest(LC_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Buffer River
LC_win_buffer <- as.owin(river_buffer_LC) #turning the buffer into a window
LC_ppp_buffer <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_buffer) #creating the poisson point pattern for lm
plot(LC_ppp_buffer, pch = 16, cex = 0.5) # plotting the randomized point pattern
LC_k_buffer <- Kest(LC_ppp_buffer, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(LC_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

# Making a nicer version of the plot for papers/presentations
plot(LC_k_buffer, main=NULL, las=1, ylab = "", xlab = "", legendargs=list(cex=0.8, xpd=TRUE), 
     yaxp = c(0, 120000, 10))  # plotting Ripley's output
par(mar = c(6,6,6,6))
title(ylab = bquote(italic("K(r), La Cobriza trees")), cex = 1.2, line = 4)
title(xlab = bquote(italic("r (m)")),)

#Ripley's K for SD
SD_win <- as.owin(SD_fixed_field_data_processed_box) #turning the box into a window
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
SD_k <- Kest(SD_ppp, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output

#Ripley's K for SD with Convex Hull
SD_win_convex <- as.owin(river_SD_convex_hull) #turning the convex hull into a window
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_convex) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5) # plotting the randomized point pattern
SD_k <- Kest(SD_ppp, correction = "Ripley") #Ripley's K function
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE))  # plotting Ripley's output
# Making a nicer version of the plot for papers/presentations
plot(SD_k, main=NULL, las=1, ylab = "", xlab = "", legendargs=list(cex=0.8, xpd=TRUE), 
     yaxp = c(0, 300000, 10), cex.axis = 1) # plotting Ripley's output
par(mar = c(8,8,8,8))
title(ylab = bquote(italic("K(r), San Dionisio trees")), cex.lab = 1.1, line = 4)
title(xlab = bquote(italic("r (m)")), cex.lab = 1.1)


# Buffer River
SD_win_buffer <- as.owin(river_buffer_SD) #turning the buffer into a window
SD_ppp_buffer <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_buffer) #creating the poisson point pattern for lm
plot(SD_ppp_buffer, pch = 16, cex = 0.5) # plotting the randomized point pattern
SD_k_buffer <- Kest(SD_ppp_buffer, correction = "Ripley") # Ripley's K function, using the isotropic/Ripley correction
# Making a nicer version of the plot for papers/presentations
plot(SD_k_buffer, main=NULL, las=1, ylab = "", xlab = "", legendargs=list(cex=0.8, xpd=TRUE), 
     yaxp = c(0, 300000, 10), cex.axis = 1) # plotting Ripley's output
par(mar = c(8,8,8,8))
title(ylab = bquote(italic("K(r), San Dionisio trees")), cex.lab = 1.1, line = 4)
title(xlab = bquote(italic("r (m)")), cex.lab = 1.1)

#### ANN Analysis (test for clustering/dispersion) ####

# For all ANN Analyses (for each population and for convex hulls/buffers/and river shapefiles):
  # 1) Find Average Nearest Neighbor Value for the trees of a population
  # 2) Simulate a randomized distribution of points, calculate and store the average nearest neighbor, 599 times
  # 3) Make sure the randomized points have the correct coordinate reference system for mapping
  # 4) Plot the Randomized Points, convex hukk, and river shapefile
  # 5) Plot the histogram of simulated randomized average nearest neighbor values with a line for the actuall ANN value
  # 6) Calculate and store the p-value 

# The convex hull windows were used because they produced the most close to the original point distances as possible.

# In this analysis, we do not control for the presence of the river as something that would influence the distribution of points

# When the average nearest neighbor is to the left of the histogram the points are more clustered than expected at random
# When the average nearest neighbor is to the right of the histogram the points are more dispersed than expected at random

## Convex Hull

# LM

#calculating the average nearest neighbor value for the entire population of trees
ann.p_LM <- mean(nndist(LM_ppp, k=1))
ann.p_LM

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), win = river_LM_convex_hull) # generating the random points within the convex hull window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the number of points at LM, it assigns a random point within the convex hull window
plot(rand.p) #plotting the random points

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and the river
ggplot()+ 
  geom_sf(data=river_LM_trans)+ #plotting the river
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points


#creating a histogram of the ANN Simulation Results
as_tibble(ann.r) %>%  #turns the list of ann values from the simulations of random points and turns it into a tibble/dataframe
  ggplot()+
    geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) +
    xlim(range(ann.p_LM, ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
    geom_vline(xintercept=ann.p_LM, col = "red") + #adds a verticle line of our tree'\s' ANN
    xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p_LM){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN (p-value)

# LC

#assigning average nearest neighbor values for the entire population of trees
ann.p_LC <- mean(nndist(LC_ppp, k=1))
ann.p_LC

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf), win = river_LC_convex_hull) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the convex hull window
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and the river
ggplot()+ 
  geom_sf(data=river_LC_trans)+ #plotting the river 
  geom_sf(data=LC_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points


#creating a histogram of the ANN Simulation Results
as_tibble(ann.r) %>%  #turns the list of ann values from the simulations of random points and turns it into a tibble/dataframe
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) +
  xlim(range(ann.p_LC, ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=ann.p_LC, col = "red") + #adds a verticle line of our tree's ANN
  xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p_LC){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

# SD

#assigning average nearest neighbor values for the entire population of trees
ann.p_SD <- mean(nndist(SD_ppp, k=1))
ann.p_SD

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(SD_fixed_field_data_processed_sf), win = river_SD_convex_hull) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the convex hull window
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and the river
ggplot()+ 
  geom_sf(data=river_SD_trans)+ #plotting the river edge raster
  geom_sf(data=SD_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#creating a histogram of the ANN Simulation Results
as_tibble(ann.r) %>%  #turns the list of ann values from the simulations of random points and turns it into a tibble/dataframe
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) +
  xlim(range(ann.p_SD, ann.r)) + #sets the limit of the x-axis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=ann.p_SD, col = "red") + #adds a verticle line of our tree's ANN
  xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p_SD){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#### ANN Analysis (test for clustering/dispersion) while controlling for the river ####

# The steps for this ANN are the same as previously, except we use three different versions of the windows in which we generate random points
# with varying levels of control for the river to see if the points still seem significantly clustered despite the 
# presence of the rivers (similar to the PPM analysis later).

# The three ways of controlling for the river include 
    # a) controlling for the river border (using a river multipoint raster window), 
    # b) controlling for on, inside, and around the river (using an inverse distance raster window), and 
    # c) controlling for on and inside the river (using a raster of the river window)


# To do this, we add new steps in the beginning 
  # 1) Create rasters of the river shapefile, river buffer, and create a raster with the inverse distances of 
        #each cell to the river shapefile (closer cells are weighted higher)
        #cells within a certain distance of the river equal 1 and the other points equals 1/distance
  # 2) Run the simulations whereby the windows either use the river border raster, the inverse distance raster where 
        #the randomized points are placed more likely based on the raster or the higher cell weights, and the river polygon raster

###  LM

# creating the rasters that will be used for point generation later

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LM_trans_outline <- st_cast(river_LM_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_LM_trans_point_raster <- st_rasterize(river_LM_trans_outline) #create raster of lake edge points
plot(river_LM_trans_point_raster)

#creating a raster of the river buffer polygon within distances can be calculated
river_LM_buffer_trans_outline <- st_cast(river_buffer_LM, "LINESTRING") #turns the polyline of the river into a multipoint object
river_buffer_LM_point_raster <- st_rasterize(river_LM_buffer_trans_outline) #create raster of lake edge points, actually a stars object
plot(river_buffer_LM_point_raster)


#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, LM_fixed_field_data_processed_sf, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run, but it depends on the computer
dist_near_river_buffer_LM_inverse <- 1/dist_near_river_buffer_LM #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LM_inverse)

#creating a raster with assigned values of 1 to cells within 70 m of the river edge and 1/distance to the cells outside to turn the distances into values 0-1
dist_near_river_buffer_LM_inverse <- dist_near_river_buffer_LM %>% #creating a new stars object with new defined values for distance
  st_as_sf() %>% #converting the stars to a shapefile
  mutate(d = case_when(d <= 70 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 70 m away from rivers edge with value of 1 and taking 1/distance for all other cells
  st_rasterize() #convert the shapefile into a raster
plot(dist_near_river_buffer_LM_inverse)


## Version of ANN analysis controlling for the river with just the river multipoint 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), f = as.im(river_LM_trans_point_raster)) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point on top of the river's edge while controlling for the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=river_LM_trans_point_raster)+ #plotting the river edge raster
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs$geom, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_LM, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LM, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_LM){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with inside, on, and outside the river

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf),
                   f = as.im(dist_near_river_buffer_LM_inverse)) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the distance raster while controlling for distance to the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=dist_near_river_buffer_LM_inverse)+ #plotting the distance inverse raster 
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_LM, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LM, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_LM){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with on and inside the river 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), f = as.im(st_rasterize(river_LM_trans))) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the raster while controlling for the river
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and river raster
ggplot()+ 
  geom_stars(data=st_rasterize(river_LM_trans))+ #plotting the river raster 
  geom_sf(data=LM_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_LM, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LM, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_LM){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

### LC

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LC_trans_points <- st_cast(river_LC_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_LC_trans_point_raster <- st_rasterize(river_LC_trans_points) #create raster of lake edge points
plot(river_LC_trans_point_raster)

river_buffer_LC_points <- st_cast(river_buffer_LC, "LINESTRING") #turns the polyline of the river buffer into a multipoint object in stars
river_buffer_LC_point_raster <- st_rasterize(river_buffer_LC_points) #create raster of lake edge points
plot(river_buffer_LC_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LC_point_raster[is.na(river_buffer_LC_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LC <- dist_to_nearest(river_buffer_LC_point_raster, river_LC_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
dist_near_river_buffer_LC_inverse <- 1/dist_near_river_buffer_LC #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LC_inverse)

#creating a raster with assigned values of 1 to cells within 30 m of the river edge and 1/distance to the cells outside to turn the distances into values 0-1
dist_near_river_buffer_LC_inverse <- dist_near_river_buffer_LC %>% #creating a new stars object with new defined values for distance
  st_as_sf() %>% #converting the stars to a shapefile
  mutate(d = case_when(d <= 20 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 20 m away from rivers edge with value of 1 and taking 1/distance for all other cells
  st_rasterize() #convert the shapefile into a raster
plot(dist_near_river_buffer_LC_inverse)

## Version of ANN analysis controlling for the river with just the river multipoint 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf), f = as.im(river_LC_trans_point_raster)) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point on top of the river's edge while controlling for the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=river_LC_trans_point_raster)+ #plotting the river edge raster
  geom_sf(data=LC_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_LC, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LC, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with inside, on, and outside the river

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf),
                   f = as.im(dist_near_river_buffer_LC_inverse)) # generating the random points within the window
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the distance raster while controlling for distance to the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=dist_near_river_buffer_LC_inverse)+ #plotting the distance inverse raster 
  geom_sf(data=LC_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_LC, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LC, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with on and inside the river 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf), f = as.im(st_rasterize(river_LC_trans))) #dist_near_river_buffer_LM_inverse 
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the raster while controlling for the river
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and river raster
ggplot()+ 
  geom_stars(data=st_rasterize(river_LC_trans))+ #plotting the river raster 
  geom_sf(data=LC_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(ann.p_LC, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_LC, col = "red", size = 1.2) + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor (ANN)") +
  theme_classic()+
  theme(axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15))

#calculating pseudo p-value for the disperse mean ANN value
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_LC){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

###test for SD

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_SD_trans_points <- st_cast(river_SD_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_SD_trans_point_raster <- st_rasterize(river_SD_trans_points) #create raster of lake edge points
plot(river_SD_trans_point_raster)

river_buffer_SD_points <- st_cast(river_buffer_SD, "LINESTRING") #turns the polyline of the river buffer into a multipoint object
river_buffer_SD_point_raster <- st_rasterize(river_buffer_SD_points) #create raster of lake edge points
plot(river_buffer_SD_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_SD_point_raster[is.na(river_buffer_SD_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_SD <- dist_to_nearest(river_buffer_SD_point_raster, river_SD_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
dist_near_river_buffer_SD_inverse <- 1/dist_near_river_buffer_SD #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_SD_inverse)

#creating a raster with assigned values of 1 to cells within 50 m of the river edge and 1/distance to the cells outside to turn the distances into values 0-1
dist_near_river_buffer_SD_inverse <- dist_near_river_buffer_SD %>% #creating a new stars object with new defined values for distance
  st_as_sf() %>% #converting the stars to a shapefile
  mutate(d = case_when(d <= 50 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 50 m away from rivers edge with value of 1 and taking 1/distance for all other cells
  st_rasterize() #convert the shapefile into a raster
plot(dist_near_river_buffer_SD_inverse)

## Version of ANN analysis controlling for the river with just the river multipoint 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(SD_fixed_field_data_processed_sf), f = as.im(river_SD_trans_point_raster)) 
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point on top of the river's edge while controlling for the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#assigning average nearest neighbor values for the entire population of trees
ann.p_SD <- mean(nndist(SD_ppp, k=1))
ann.p_SD

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=river_SD_trans_point_raster)+ #plotting the river edge raster
  geom_sf(data=SD_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_SD, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_SD, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_SD){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with inside, on, and outside the river

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(SD_fixed_field_data_processed_sf),
                   f = as.im(dist_near_river_buffer_SD_inverse)) #dist_near_river_buffer_LM_inverse #forcewin = T, win=as.owin(river_LM_convex_hull)
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the distance raster while controlling for distance to the river's edge
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=dist_near_river_buffer_SD_inverse)+ #plotting the distance inverse raster 
  geom_sf(data=SD_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA)+ #plotting the random points
  labs(color = "Trees", fill = "Inverse Distance (m)")
  
  
#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p_SD, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_SD, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_SD){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN


## Version of ANN analysis controlling for the river with on and inside the river 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(SD_fixed_field_data_processed_sf), f = as.im(st_rasterize(river_SD_trans))) #assigns a random point for the number of trees in SD favoring placements in the river raster
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the raster while controlling for the river
plot(rand.p)

#adding the UTM 12 crs to rand.p
rand.p.crs <- rand.p %>% 
  st_as_sf()%>%
  st_set_crs(26912)

#plotting the randomly generated points, tree points, and river raster
ggplot()+ 
  geom_stars(data=st_rasterize(river_SD_trans))+ #plotting the river raster 
  geom_sf(data=SD_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points


as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(ann.p_SD, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p_SD, col = "red", size = 1.2) + #plotting our tree's mean ANN
  xlab("Average Nearest Neighbor (ANN)") +
  theme_classic()+
  theme(axis.text=element_text(size=15),  axis.title.x =element_text(size= 15),
        axis.title.y =element_text(size= 15))

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p_SD){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#### PPM analysis ####

# for every PPM analysis 

  # 1) generate an image object of inverse distance raster 
  # 2) create the Poisson Point Model using ppm() function with the river influencing the location of the points (Alternative Hypothesis Model)
  # 3) create the Poisson Point Model using ppm() function with the river not influencing the location of the points (Null Hypothesis Model)
  # 4) Use an ANOVA likelihood Ratio Test to compare the Alternate and Null hypotheses
  # 5) Plot the influence of the river as the distance to the river decreases (inverse distance)

#Test for LM

#creating the image of the distance to river stars
dist_near_river_buffer_LM_inverse_im <- as.im(dist_near_river_buffer_LM_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(LM_fixed_field_data_processed_sf) ~ dist_near_river_buffer_LM_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(LM_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LM_inverse_im", se.fit = TRUE), main = "Distance to River of Las Matancitas",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)

#Test for LC

#creating the image of the distance to river stars
dist_near_river_buffer_LC_inverse_im <- as.im(dist_near_river_buffer_LC_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(LC_fixed_field_data_processed_sf) ~ dist_near_river_buffer_LC_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(LC_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LC_inverse_im", se.fit = TRUE), main = "Distance to River of La Cobriza",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)




#Test for SD

#creating the image of the distance to river stars
dist_near_river_buffer_SD_inverse_im <- as.im(dist_near_river_buffer_SD_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(SD_fixed_field_data_processed_sf) ~ dist_near_river_buffer_SD_inverse_im) 
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(SD_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_SD_inverse_im", se.fit = TRUE), main = "Distance to River of San Dionisio",
     ylab = "Quercus brandegeei Trees", xlab = "Inverse Distance to River", legend = FALSE)


# making examples of random point distributions vs. points only along the river's edge for presentation 

points_box = sf::st_sample(SD_box_sf, size=50) #randomizing points onlu in population bbox
points_river = sf::st_sample(river_SD_trans_points, size=50) #randomizing points along river's edge

points <- st_as_sf(points, crs = 26912) #projecting the points

#plotting the randomized box points
ggplot()+
  geom_sf(data = river_SD_trans)+
  #geom_sf(data = SD_box_sf)+
  geom_sf(data = points_box, size = 2)+
  theme_classic()

#plotting the randomized river's edge points
ggplot()+
  geom_sf(data = river_SD_trans)+
  #geom_sf(data = SD_box_sf)+
  geom_sf(data = points_river, size = 2)+
  theme_classic()

