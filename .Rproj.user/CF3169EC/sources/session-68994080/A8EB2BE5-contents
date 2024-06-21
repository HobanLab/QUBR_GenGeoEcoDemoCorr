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


fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#upload river shapefile and filter out polygons for each population
rivers <- st_read("./data/QUBR Rivers and Trees.kml", "Rivers", crs = 4326)
rivers_2d <- st_zm(rivers, drop = T) #we had a z dimension with max and min, so we got rid of it because it was giving us weird errors and disrupting later statistics
river_LC <- filter(rivers_2d, Name == "River LC")
river_SD <- filter(rivers_2d, Name == "River SD")
river_LM <- filter(rivers_2d, Name == "LM River")

river_LM_trans <- st_transform(river_LM, crs = 26912) #equal area projection, uses meters as distance measurement
river_LC_trans <- st_transform(river_LC, crs = 26912)
river_SD_trans <- st_transform(river_SD, crs = 26912)

#### Creating 20 m buffer around river polygons ####

river_buffer_LM<- st_buffer(river_LM_trans, 200)
ggplot(river_buffer_LM)+
  geom_sf()

river_buffer_LC<- st_buffer(river_LC_trans, 230)
ggplot(river_buffer_LC)+
  geom_sf()

river_buffer_SD<- st_buffer(river_SD_trans, 120)
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

#creating entire boundary shapefile
fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  st_as_sfc()

fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  st_bbox %>%
  st_as_sfc()

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


#### Creating Convex Hulls ####
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
win <- as.owin(fixed_field_data_processed_box)
ppp <- as.ppp(st_coordinates(fixed_field_data_processed_sf), W = win) #creating the poisson point pattern for lm
plot(ppp, pch = 16, cex = 0.5)
K <- Kest(ppp, correction = "Ripley")
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM 
LM_win <- as.owin(LM_fixed_field_data_processed_box)
LM_ppp <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win) #creating the poisson point pattern for lm
plot(LM_ppp, pch = 16, cex = 0.5)
LM_k <- Kest(LM_ppp, correction = "Ripley")
plot(LM_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM with Convex Hull
LM_win_convex <- as.owin(river_LM_convex_hull)
LM_ppp_convex <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_convex) #creating the poisson point pattern for lm
plot(LM_ppp_convex, pch = 16, cex = 0.5)
LM_k_convex <- Kest(LM_ppp_convex, correction = "Ripley")
plot(LM_k_convex, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LM with Buffer River 20 m
LM_win_buffer <- as.owin(river_buffer_LM)
LM_ppp_buffer <- as.ppp(st_coordinates(LM_fixed_field_data_processed_sf), W = LM_win_buffer) #creating the poisson point pattern for lm
plot(LM_ppp_buffer, pch = 16, cex = 0.5)
LM_k_buffer <- Kest(LM_ppp_buffer, correction = "Ripley")
plot(LM_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC 
LC_win <- as.owin(LC_fixed_field_data_processed_box)
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5)
LC_k <- Kest(LC_ppp, correction = "Ripley")
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC with Convex Hull
LC_win_convex <- as.owin(river_LC_convex_hull)
LC_ppp <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_convex) #creating the poisson point pattern for lm
plot(LC_ppp, pch = 16, cex = 0.5)
LC_k <- Kest(LC_ppp, correction = "Ripley")
plot(LC_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for LC with Buffer River 20 m
LC_win_buffer <- as.owin(river_buffer_LC)
LC_ppp_buffer <- as.ppp(st_coordinates(LC_fixed_field_data_processed_sf), W = LC_win_buffer) #creating the poisson point pattern for lm
plot(LC_ppp_buffer, pch = 16, cex = 0.5)
LC_k_buffer <- Kest(LC_ppp_buffer, correction = "Ripley")
plot(LC_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD
SD_win <- as.owin(SD_fixed_field_data_processed_box)
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5)
SD_k <- Kest(LC_ppp, correction = "Ripley")
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD with Convex Hull
SD_win_convex <- as.owin(river_SD_convex_hull)
SD_ppp <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_convex) #creating the poisson point pattern for lm
plot(SD_ppp, pch = 16, cex = 0.5)
SD_k <- Kest(SD_ppp, correction = "Ripley")
plot(SD_k, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#Ripley's K for SD with Buffer River 20 m
SD_win_buffer <- as.owin(river_buffer_SD)
SD_ppp_buffer <- as.ppp(st_coordinates(SD_fixed_field_data_processed_sf), W = SD_win_buffer) #creating the poisson point pattern for lm
plot(SD_ppp_buffer, pch = 16, cex = 0.5)
SD_k_buffer <- Kest(SD_ppp_buffer, correction = "Ripley")
plot(SD_k_buffer, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot

#### ANN Analysis (test for clustering/dispersion) ####

#test for LM
ann.p <- mean(nndist(LM_ppp, k=1))
ann.p

n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), win = river_LM_convex_hull) #river_LM_trans
  ann.r[i] <- mean(nndist(rand.p, k=1))
} #for the length of the number of points at LM, it assigns a random point within the convex hull window
plot(rand.p)

hist(ann.r, main = NULL, las=1, breaks = 40, col = "light blue", xlim = range(ann.p, ann.r))
abline(v=ann.p, col = "red")

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#to calculate 1 sided pseudo p-value the way done by Manny, I am unsure why he added 1s
N.greater <- sum(ann.r > ann.p)
p <- min(N.greater + 1, n + 1 - N.greater) / (n +1)
p

#turning river polygon into multipoints and then into a raster
river_LM_trans_points <- st_cast(river_LM_trans, "MULTIPOINT")

river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #create raster of lake edge points
plot(river_LM_trans_point_raster)
plot(as.im(river_LM_trans_point_raster))

#ANN analysis controlling for river
n <- 599L #defines the number of simulations
ann.r <- vector(length = n) #creates the empty object that we can store ANN values in
for (i in 1:n){
  rand.p <- rpoint(n=length(LM_fixed_field_data_processed_sf), win = river_LM_convex_hull,
                   f = as.im(river_LM_trans_point_raster))
  ann.r[i] <- mean(nndist(rand.p, k=1))
} #for the length of the number of points at LM, it assigns a random point within the convex hull window while controlling for the river's edge
plot(rand.p)

Window(rand.p) <- as.owin(river_LM_trans)
plot(rand.p)
ggplot()+ # ask about the ggplot
  geom_sf(LM_fixed_field_data_processed_sf)+
  geom_sf(river_LM_trans)
plot(LM_fixed_field_data_processed_sf)
hist(ann.r, main = NULL, las=1, breaks = 40, col = "light blue", xlim = range(ann.p, ann.r))
abline(v=ann.p, col = "red")

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(ann.r)){
  if (ann.r[i] < ann.p){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(ann.r)) #the proportion of random ANNs that are less than our ANN

#test for LC



#Test for SD

#### PPM analysis ####

PPM1 <- ppm(as.ppp(LM_fixed_field_data_processed_sf) ~ log(river_LM_trans_points))


