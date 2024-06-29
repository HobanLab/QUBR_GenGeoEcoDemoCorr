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
library(geostatsp)
library(tmaptools)

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#upload river shapefile and filter out polygons for each population
rivers <- st_read("./data/QUBR Rivers and Trees.kml", "Rivers", crs = 4326)
rivers_2d <- st_zm(rivers, drop = T) #we had a z dimension with max and min, so we got rid of it because it was giving us weird errors and disrupting later statistics
river_LC <- filter(rivers_2d, Name == "River LC")
river_SD <- filter(rivers_2d, Name == "River SD")
river_LM <- filter(rivers_2d, Name == "LM River")

#changing the coordinate reference system of the river polygons to be equal area projection (UTM 12N), uses meters as distance measurement 
river_LM_trans <- st_transform(river_LM, crs = 26912) 
river_LC_trans <- st_transform(river_LC, crs = 26912)
river_SD_trans <- st_transform(river_SD, crs = 26912)

#### Creating buffer around river polygons ####

river_buffer_LM<- st_buffer(river_LM_trans, 200) #200 m buffer
ggplot(river_buffer_LM)+
  geom_sf()

river_buffer_LC<- st_buffer(river_LC_trans, 230) #230 m buffer
ggplot(river_buffer_LC)+
  geom_sf()

river_buffer_SD<- st_buffer(river_SD_trans, 120) #120 m buffer
ggplot(river_buffer_SD)+
  geom_sf()

#### Creating fixed_field_data_processed dataframes for each population ####

LM_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "SD")

#### Importing Shapefile #### 

#turn the BCS polygon into a shapefile and visualize its outline
BCS_polygon <- read_sf("./data/Shapefiles/BCS_Polygon/bcs_entidad.shp")
BCS_polygon <- st_as_sf(BCS_polygon)
plot(BCS_polygon$geometry)

#creating a shapefile of all points with lat lon coordinates in WGS 1984
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                       coords = c("long", "lat"), crs = 4326)

#creating a transformed shapefile with UTM 12 N an equal area projection
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) 

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


#creating BCS boundary shapefile, turning sf of all points into sfc
fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  st_as_sfc()

#creating a boundry box with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  st_bbox %>%
  st_as_sfc()

#finding minimum and maximum lat and long values for LM
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

#creating LM boundary shapefile, turning sf of all points into sfc
LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sfc()

#creating a boundry box of LM with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
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

#creating LC boundary shapefile, turning sf of all points into sfc
LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sfc()

#creating a boundry box of LC with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
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

#creating SD boundary shapefile, turning sf of all points into sfc
SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sfc()

#creating a boundry box of SD with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
SD_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_bbox %>%
  st_as_sfc()


#### Creating Convex Hulls using tree points of each population ####
river_LM_convex_hull <- st_convex_hull(st_union(LM_fixed_field_data_processed_sf)) #LM_fixed_field_data_processed_sf
ggplot(river_LM_convex_hull)+
  geom_sf()

river_LC_convex_hull <- st_convex_hull(st_union(LC_fixed_field_data_processed_sf)) #LM_fixed_field_data_processed_sf
ggplot(river_LC_convex_hull)+
  geom_sf()

river_SD_convex_hull <- st_convex_hull(st_union(SD_fixed_field_data_processed_sf)) #LM_fixed_field_data_processed_sf
ggplot(river_SD_convex_hull)+
  geom_sf()

#### Ripley's K Analysis (version with box, convex hull, and 20 m buffer around river) ####

#Ripley's K for all points 
win <- as.owin(fixed_field_data_processed_box) #turning the box into a window
ppp <- as.ppp(st_coordinates(fixed_field_data_processed_sf), W = win) #creating the poisson point pattern for lm
plot(ppp, pch = 16, cex = 0.5)
K <- Kest(ppp, correction = "Ripley") #Ripley's K function
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM 
LM_win <- as.owin(LM_fixed_field_data_processed_box) #turning the box into a window
LM_ppp <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win) #creating the poisson point pattern for lm
plot(LM_ppp, pch = 16, cex = 0.5)
LM_k <- Kest(LM_ppp, correction = "Ripley") #Ripley's K function
plot(LM_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM with Convex Hull
LM_win_convex <- as.owin(river_LM_convex_hull)  #turning the convex hull into a window
LM_ppp_convex <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_convex) #creating the poisson point pattern for lm
plot(LM_ppp_convex, pch = 16, cex = 0.5)
LM_k_convex <- Kest(LM_ppp_convex, correction = "Ripley") #Ripley's K function
plot(LM_k_convex, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM with Buffer River 
LM_win_buffer <- as.owin(river_buffer_LM) #turning the buffer into a window
LM_ppp_buffer <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_buffer) #creating the poisson point pattern for lm
plot(LM_ppp_buffer, pch = 16, cex = 0.5)
LM_k_buffer <- Kest(LM_ppp_buffer, correction = "Ripley") #Ripley's K function
plot(LM_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC 
LC_win <- as.owin(LC_fixed_field_data_processed_box) #turning the box into a window
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5)
LC_k <- Kest(LC_ppp, correction = "Ripley") #Ripley's K function
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC with Convex Hull
LC_win_convex <- as.owin(river_LC_convex_hull) #turning the convex hull into a window
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_convex) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5)
LC_k <- Kest(LC_ppp, correction = "Ripley") #Ripley's K function
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC with Buffer River
LC_win_buffer <- as.owin(river_buffer_LC) #turning the buffer into a window
LC_ppp_buffer <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_buffer) #creating the poisson point pattern for lm
plot(LC_ppp_buffer, pch = 16, cex = 0.5)
LC_k_buffer <- Kest(LC_ppp_buffer, correction = "Ripley") #Ripley's K function
plot(LC_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD
SD_win <- as.owin(SD_fixed_field_data_processed_box) #turning the box into a window
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5)
SD_k <- Kest(LC_ppp, correction = "Ripley") #Ripley's K function
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD with Convex Hull
SD_win_convex <- as.owin(river_SD_convex_hull) #turning the convex hull into a window
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_convex) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5)
SD_k <- Kest(SD_ppp, correction = "Ripley") #Ripley's K function
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD with Buffer River
SD_win_buffer <- as.owin(river_buffer_SD) #turning the buffer into a window
SD_ppp_buffer <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_buffer) #creating the poisson point pattern for lm
plot(SD_ppp_buffer, pch = 16, cex = 0.5)
SD_k_buffer <- Kest(SD_ppp_buffer, correction = "Ripley") #Ripley's K function
plot(SD_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#### ANN Analysis (test for clustering/dispersion) ####

#test for LM

#assigning average nearest neighbor values for the entire population of trees
ann.p_LM <- mean(nndist(LM_ppp, k=1))
ann.p_LM

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), win = river_LM_convex_hull) #river_buffer_LM  #river_LM_trans. #river_LM_convex_hull
  ann.r[i] <- mean(nndist(rand.p, k=1)) #for each simulated random distribution of points it calculates the mean ANN across all of the trees
} #for the length of the number of points at LM, it assigns a random point within the convex hull window
plot(rand.p)

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
    xlim(range(ann.p, ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
    geom_vline(xintercept=ann.p, col = "red") + #adds a verticle line of our tree'\s' ANN
    xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#Test for LC

#assigning average nearest neighbor values for the entire population of trees
ann.p_LC <- mean(nndist(LC_ppp, k=1))
ann.p_LC

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf), win = river_LC_convex_hull) #river_buffer_LM  #river_LM_trans. #river_LM_convex_hull
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
  xlim(range(ann.p, ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=ann.p, col = "red") + #adds a verticle line of our tree'\s' ANN
  xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#Test for SD

#assigning average nearest neighbor values for the entire population of trees
ann.p_SD <- mean(nndist(SD_ppp, k=1))
ann.p_SD

#simulation to create a list of ANN from randomly placed points
n <- 566L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(SD_fixed_field_data_processed_sf), win = river_SD_convex_hull) #river_buffer_LM  #river_LM_trans. #river_LM_convex_hull
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
  xlim(range(ann.p, ann.r)) + #sets the limit of the xaxis to encompass the ANN for our trees and histogram of ANNs from the simulation
  geom_vline(xintercept=ann.p, col = "red") + #adds a verticle line of our tree'\s' ANN
  xlab("ANN")+
  theme_classic()

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#### ANN Analysis (test for clustering/dispersion) while controlling for the river ####

###Test for LM

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LM_trans_points <- st_cast(river_LM_trans, "MULTIPOINT") #turns the polyline of the river into a multipoint object
river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #create raster of lake edge points
plot(river_LM_trans_point_raster)

river_buffer_LM_points <- st_cast(river_buffer_LM, "MULTIPOINT") #turns the polyline of the river buffer into a multipoint object
river_buffer_LM_point_raster <- st_rasterize(river_buffer_LM_points) #create raster of lake edge points
plot(river_buffer_LM_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, river_LM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
dist_near_river_buffer_LM_inverse <- 1/dist_near_river_buffer_LM #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LM_inverse)

#creating a raster with assigned values of 1 to cells within 30 m of the river edge and 1/distance to the cells outside to turn the distances into values 0-1
dist_near_river_buffer_LM_inverse <- dist_near_river_buffer_LM %>% #creating a new stars object with new defined values for distance
  st_as_sf() %>% #converting the stars to a shapefile
  mutate(d = case_when(d <= 30 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 30 m away from rivers edge with value of 1 and taking 1/distance for all other cells
  st_rasterize() #convert the shapefile into a raster
plot(dist_near_river_buffer_LM_inverse)

## Version of ANN analysis controlling for the river with just the river multipoint 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), f = as.im(river_LM_trans_point_raster)) 
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
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf),
                   f = as.im(dist_near_river_buffer_LM_inverse)) #dist_near_river_buffer_LM_inverse #forcewin = T, win=as.owin(river_LM_convex_hull)
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
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), f = as.im(st_rasterize(river_LM_trans))) #dist_near_river_buffer_LM_inverse 
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
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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

###test for LC

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LC_trans_points <- st_cast(river_LC_trans, "MULTIPOINT") #turns the polyline of the river into a multipoint object
river_LC_trans_point_raster <- st_rasterize(river_LC_trans_points) #create raster of lake edge points
plot(river_LC_trans_point_raster)

river_buffer_LC_points <- st_cast(river_buffer_LC, "MULTIPOINT") #turns the polyline of the river buffer into a multipoint object in stars
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
  mutate(d = case_when(d <= 30 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 30 m away from rivers edge with value of 1 and taking 1/distance for all other cells
  st_rasterize() #convert the shapefile into a raster
plot(dist_near_river_buffer_LC_inverse)

## Version of ANN analysis controlling for the river with just the river multipoint 

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){ 
  rand.p <- rpoint(n=length(LC_fixed_field_data_processed_sf), f = as.im(river_LC_trans_point_raster)) 
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
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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
                   f = as.im(dist_near_river_buffer_LC_inverse)) #dist_near_river_buffer_LM_inverse #forcewin = T, win=as.owin(river_LM_convex_hull)
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
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
  xlab("ANN") +
  theme_classic()

#calculating pseudo p-value for the disperse mean ANN value
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] > ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

###test for SD

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_SD_trans_points <- st_cast(river_SD_trans, "MULTIPOINT") #turns the polyline of the river into a multipoint object
river_SD_trans_point_raster <- st_rasterize(river_SD_trans_points) #create raster of lake edge points
plot(river_SD_trans_point_raster)

river_buffer_SD_points <- st_cast(river_buffer_SD, "MULTIPOINT") #turns the polyline of the river buffer into a multipoint object
river_buffer_SD_point_raster <- st_rasterize(river_buffer_SD_points) #create raster of lake edge points
plot(river_buffer_SD_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_SD_point_raster[is.na(river_buffer_SD_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_SD <- dist_to_nearest(river_buffer_SD_point_raster, river_SD_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
dist_near_river_buffer_SD_inverse <- 1/dist_near_river_buffer_SD #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_SD_inverse)

#creating a raster with assigned values of 1 to cells within 30 m of the river edge and 1/distance to the cells outside to turn the distances into values 0-1
dist_near_river_buffer_SD_inverse <- dist_near_river_buffer_SD %>% #creating a new stars object with new defined values for distance
  st_as_sf() %>% #converting the stars to a shapefile
  mutate(d = case_when(d <= 60 ~ 1, 
                       d > 1 ~ 1/d)) %>% #assigning cells less than 30 m away from rivers edge with value of 1 and taking 1/distance for all other cells
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

#plotting the randomly generated points, tree points, and probability/distance raster
ggplot()+ 
  geom_stars(data=river_SD_trans_point_raster)+ #plotting the river edge raster
  geom_sf(data=SD_fixed_field_data_processed_sf, aes(col = "red"))+ #plotting the tree points
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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
  geom_sf(data=rand.p.crs, fill = NA) #plotting the random points

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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

#graphing the histogram of simulated ANN values and the mean ANN from our trees
as_tibble(ann.r) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "dodgerblue1", color = "black", bins = 50) + 
  xlim(range(ann.p, ann.r)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=ann.p, col = "red") + #plotting our tree's mean ANN
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

#### PPM analysis ####

#Test for LM

#creating the image of the distance to river stars
dist_near_river_buffer_LM_inverse_im <- as.im(dist_near_river_buffer_LM_inverse)

#Alternative hypothesis, seeing if the distance to the river's edge influences the tree point placement
PPM1 <- ppm(Q = as.ppp(LM_fixed_field_data_processed_sf) ~ dist_near_river_buffer_LM_inverse_im) #as.im(dist_near_river_buffer_LM_inverse))
PPM1

#null hypothesis, no change in the trend of the points
PPM0 <- ppm(as.ppp(LM_fixed_field_data_processed_sf) ~ 1)
PPM0

#using a likelihood ratio test to compare the alternative and null models
anova(PPM0, PPM1, test="LRT")

#plotting the alternative model
plot(effectfun(PPM1, "dist_near_river_buffer_LM_inverse_im", se.fit = TRUE), main = "Distance to River of Las Matancitas",
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)

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
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)

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
     ylab = "Quercus brandegeei Trees", xlab = "Distance to River", legend = FALSE)
