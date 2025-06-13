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
library(mgcv) #to use gam function 
library(plotly) #to 3d plot variables
library(MuMIn) #to be able to use dredge
library(visreg) # to be able to plot Aspect/categorical variables with GAM

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#transforming the shapefile of trees from WGS84 into equal area projection UTM 12N
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection


#create dataframe with X and Y UTM coordinates
fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with seperate x and y columns from the UTM 12N transformation
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
  cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

#export the csv of the UTM 12N points for using the file in ArcGIS to make new shapefiles
fixed_field_data_processed_sf_trans_coordinates_download <- write.csv(fixed_field_data_processed_sf_trans_coordinates, "./data/fixed_field_data_processed_sf_trans_coordinates.csv", row.names = F)
View(fixed_field_data_processed_sf_transformed)

#creating shapefiles for each population, turning sf of all points into sfc

LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sfc()

LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sfc()

SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sfc()

#transformations of LM variables (log, square root, ) for linear models

#creating columns with transformations: logged all of the variables
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")

#set elevation as a numeric value
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

#creating a new column in the whole dataset to get rid of  360 m outlier and turn the values in feet into meter
fixed_field_data_processed_sf_trans_coordinates <-  fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. < 700 & Elevation..m. != 360) ~ Elevation..m.,
                                        (Elevation..m. == 360) ~ NA, 
                                        (Elevation..m. > 700) ~ Elevation..m.*0.3048))  #because LM and LC do not have a 360 elevation and SD and LC do have values above 700, this should not effect them


#creating a new elevation column so the values that were mistakenly put in feet are in meters
LM_fixed_field_data_processed <-  LM_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))

#creating a new elevation column so LC, LM, and SD all have this same column, makes it easier for combining the populaiton data frames
LC_fixed_field_data_processed <-  LC_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))



#upload ArcGIS river shapefile and filter out polygons for each population
river_LM <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LM River/LM_Rivers_Final.shp")
river_LM  <- river_LM$geometry[1]
plot(river_LM)

river_LC  <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LC River/LC_Rivers_Final.shp")
river_LC  <- river_LC$geometry[1]
plot(river_LC)

river_SD <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/SD River/SD_Rivers_Final.shp")
river_SD <- river_SD$geometry[1]
plot(river_SD)

#changing the coordinate reference system of the river polygons to be equal area projection (UTM 12N), uses meters as distance measurement 
river_LM_trans <- st_as_sf(st_transform(river_LM, crs = 26912))
river_LC_trans <- st_as_sf(st_transform(river_LC, crs = 26912))
river_SD_trans <- st_as_sf(st_transform(river_SD, crs = 26912))

#creating buffers around the rivers
river_buffer_LM <- st_buffer(river_LM_trans, 100) #100 m buffer
ggplot()+
  geom_sf(data = river_buffer_LM)+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed_sf)

river_buffer_LC<- st_buffer(river_LC_trans, 100) #230 m buffer
ggplot()+
  geom_sf(data = river_buffer_LC)+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_sf)

river_buffer_SD<- st_buffer(river_SD_trans, 70) #70 m buffer
ggplot()+
  geom_sf(data = river_buffer_SD)+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_sf)


#creating bounding boxes for each population

#creating a boundry box of LM with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
LM_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_bbox %>%
  st_as_sfc()

#creating a boundry box of LC with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
LC_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_bbox %>%
  st_as_sfc()

#creating a boundry box of SD with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
SD_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_bbox %>%
  st_as_sfc()

#creating bboxs for all of the river shapefiles for each population
LM_box <- st_bbox(river_LM_trans)
LC_box <- st_bbox(river_LC_trans)
SD_box <- st_bbox(river_SD_trans)


###creating distance

#LM

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LM_trans_points <- st_cast(river_LM_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #create raster of lake edge points
plot(river_LM_trans_point_raster)

river_LM_buffer_trans_outline <- st_cast(river_buffer_LM, "LINESTRING") #turns the polyline of the river into a multipoint object
river_buffer_LM_point_raster <- st_rasterize(river_LM_buffer_trans_outline) #create raster of lake edge points
plot(river_buffer_LM_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, river_LM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
#dist_near_river_buffer_LM_inverse <- 1/dist_near_river_buffer_LM #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LM) #not using inverse distance


#LC

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
#dist_near_river_buffer_LC_inverse <- 1/dist_near_river_buffer_LC #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LC) #not using inverse distance



#SD

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
#dist_near_river_buffer_SD_inverse <- 1/dist_near_river_buffer_SD #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_SD) #not using inverse distance


#making it so points within or overlapping with the river are assigned 1 

#Assigning points within and overlapping with the river to have true

#LM

LM_points_intersects_river <- st_intersects(LM_fixed_field_data_processed, river_LM_trans, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
LM_fixed_field_data_processed_intersects_river <- cbind(LM_fixed_field_data_processed, LM_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_LM_trans)+
  geom_sf(data=LM_fixed_field_data_processed)+
  geom_sf(data=LM_fixed_field_data_processed_intersects_river, aes(color = LM_points_intersects_river))

#LC 

LC_points_intersects_river <- st_intersects(LC_fixed_field_data_processed, river_LC_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
LC_fixed_field_data_processed_intersects_river <- cbind(LC_fixed_field_data_processed, LC_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_LC_trans)+
  geom_sf(data=LC_fixed_field_data_processed)+
  geom_sf(data=LC_fixed_field_data_processed_intersects_river, aes(color = LC_points_intersects_river))

#SD

SD_points_intersects_river <- st_intersects(SD_fixed_field_data_processed, river_SD_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
SD_fixed_field_data_processed_intersects_river <- cbind(SD_fixed_field_data_processed, SD_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_SD_trans)+
  geom_sf(data=SD_fixed_field_data_processed)+
  geom_sf(data=SD_fixed_field_data_processed_intersects_river, aes(color = SD_points_intersects_river))



###extract distance to river for each point 

#all points


#LM
LM_distance_data_pts <- st_extract(dist_near_river_buffer_LM, LM_fixed_field_data_processed) #extracting aspect for each point value
LM_fixed_field_data_processed_distance  <- cbind(LM_fixed_field_data_processed, LM_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


#LC
LC_distance_data_pts <- st_extract(dist_near_river_buffer_LC, LC_fixed_field_data_processed) #extracting aspect for each point value
LC_fixed_field_data_processed_distance  <- cbind(LC_fixed_field_data_processed, LC_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


#SD
SD_distance_data_pts <- st_extract(dist_near_river_buffer_SD, SD_fixed_field_data_processed) #extracting aspect for each point value
SD_fixed_field_data_processed_distance  <- cbind(SD_fixed_field_data_processed, SD_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


### Assigning all points within/overlapping river to distances of 0

#LM
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>% 
  mutate(d = case_when((LM_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LM_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 


#LC
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(d = case_when((LC_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LC_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 


#SD

SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(d = case_when((SD_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (SD_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value




#transformations of variables (log, square root, ) for linear models


#LM

#creating columns with transformations: logged all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

#LC

#creating columns with transformations: logged all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

#SD

#creating columns with transformations: logged all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

plot(CEM_15_utm_SD)

#Importing the cropped rasters for LM, LC, and SD
CEM_15_utm_LM <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LM.tif"))
CEM_15_utm_LC <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LC.tif"))
CEM_15_utm_SD <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_SD.tif"))

#creating the all points raster by merging the LM, LC, and SD rasters
CEM_15_utm_all_points <- raster::merge(CEM_15_utm_LM, CEM_15_utm_LC, CEM_15_utm_SD)

ggplot()+
  geom_raster(data= as.data.frame(CEM_15_utm_all_points, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = fixed_field_data_processed_sf_transformed)

## Extracting the slope 

#all points 

#extracting the slope in degrees, using the queens method (neighbor = 8)
all_points_slope_raster_15 <- terra::terrain(CEM_15_utm_all_points, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(all_points_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = fixed_field_data_processed_sf_trans_coordinates)+
  scale_fill_viridis_c()


#LM

#extracting the slope in degrees, using the queens method (neighbor = 8)
LM_slope_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LM_fixed_field_data_processed)+
  scale_fill_viridis_c()

#LC

#extracting the slope in degrees, using the queens method (neighbor = 8)
LC_slope_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LC_fixed_field_data_processed)+
  scale_fill_viridis_c()

#SD

#extracting the slope in degrees, using the queens method (neighbor = 8)
SD_slope_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = SD_fixed_field_data_processed)+
  scale_fill_viridis_c()

## Extracting the aspect 

#all points 

#extracting the slope in degrees, using the queens method (neighbor = 8)
all_points_aspect_raster_15 <- terra::terrain(CEM_15_utm_all_points, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(all_points_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = fixed_field_data_processed_sf_trans_coordinates)+
  scale_fill_viridis_c()


#LM

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LM_aspect_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LM_fixed_field_data_processed)+
  scale_fill_viridis_c()

#LC

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LC_aspect_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LC_fixed_field_data_processed)+
  scale_fill_viridis_c()

#SD

#extracting the aspect in degrees, using the queens method (neighbor = 8)
SD_aspect_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'aspect')


#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = SD_fixed_field_data_processed)+
  scale_fill_viridis_c()


#creating dataframes for each population and the slope and aspect data by extracting the slope and aspect data fromk each cell for each point and combining the data into a dataframe


#all points
all_points_aspect_raster_15_data_pts <- extract(all_points_aspect_raster_15, fixed_field_data_processed_sf_trans_coordinates) #extracting aspect for each point value
all_points_slope_raster_15_data_pts <- extract(all_points_slope_raster_15, fixed_field_data_processed_sf_trans_coordinates) #extracting slope for each point value
all_points_fixed_field_data_processed_terrain <- cbind(fixed_field_data_processed_sf_trans_coordinates, all_points_aspect_raster_15_data_pts) #bind the aspect data for each point to the LM point dataframe
all_points_fixed_field_data_processed_terrain <- cbind(all_points_fixed_field_data_processed_terrain, all_points_slope_raster_15_data_pts) #bind the slope data for each point to the LM point dataframe

View(all_points_fixed_field_data_processed_terrain)


#LM

LM_aspect_raster_15_data_pts <- extract(LM_aspect_raster_15, LM_fixed_field_data_processed) #extracting aspect for each point value
LM_slope_raster_15_data_pts <- extract(LM_slope_raster_15, LM_fixed_field_data_processed) #extracting slope for each point value
LM_elevation_raster_15_data_pts <- extract(CEM_15_utm_LM, LM_fixed_field_data_processed) #extracting the elevation for each point value
LM_fixed_field_data_processed_terrain <- cbind(LM_fixed_field_data_processed, LM_aspect_raster_15_data_pts) #bind the aspect data for each point to the LM point dataframe
LM_fixed_field_data_processed_terrain <- cbind(LM_fixed_field_data_processed_terrain, LM_slope_raster_15_data_pts) #bind the slope data for each point to the LM point dataframe
LM_fixed_field_data_processed_terrain <- cbind(LM_fixed_field_data_processed_terrain, LM_elevation_raster_15_data_pts) #bind the elevation data for each point to the LM point dataframe


#LC
LC_aspect_raster_15_data_pts <- extract(LC_aspect_raster_15, LC_fixed_field_data_processed) #extracting aspect for each point value
LC_slope_raster_15_data_pts <- extract(LC_slope_raster_15, LC_fixed_field_data_processed) #extracting slope for each point value
LC_elevation_raster_15_data_pts <- extract(CEM_15_utm_LC, LC_fixed_field_data_processed) #extracting the elevation for each point value
LC_fixed_field_data_processed_terrain <- cbind(LC_fixed_field_data_processed, LC_aspect_raster_15_data_pts) #bind the aspect data for each point to the SD point dataframe
LC_fixed_field_data_processed_terrain <- cbind(LC_fixed_field_data_processed_terrain, LC_slope_raster_15_data_pts) #bind the slope data for each point to the SD point dataframe
LC_fixed_field_data_processed_terrain <- cbind(LC_fixed_field_data_processed_terrain, LC_elevation_raster_15_data_pts) #bind the elevation data for each point to the LM point dataframe


View(LC_fixed_field_data_processed_terrain)


#SD
SD_aspect_raster_15_data_pts <- extract(SD_aspect_raster_15, SD_fixed_field_data_processed) #extracting aspect for each point value
SD_slope_raster_15_data_pts <- extract(SD_slope_raster_15, SD_fixed_field_data_processed) #extracting slope for each point value
SD_elevation_raster_15_data_pts <- extract(CEM_15_utm_SD, SD_fixed_field_data_processed) #extracting the elevation for each point value
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed, SD_aspect_raster_15_data_pts) #bind the aspect data for each point to the SD point dataframe
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed_terrain, SD_slope_raster_15_data_pts) #bind the slope data for each point to the SD point dataframe
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed_terrain, SD_elevation_raster_15_data_pts) #bind the elevation data for each point to the LM point dataframe

View(SD_fixed_field_data_processed_terrain)

#recategorizing the aspect data

#setting values of 360 to 0 

#all points
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts = case_when((all_points_aspect_raster_15_data_pts == "360") ~  0,
                                                          (all_points_aspect_raster_15_data_pts != "360")~ all_points_aspect_raster_15_data_pts))


#LM
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts = case_when((LM_aspect_raster_15_data_pts == "360") ~  0,
                                                  (LM_aspect_raster_15_data_pts != "360")~ LM_aspect_raster_15_data_pts))
#LC

LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts = case_when((LC_aspect_raster_15_data_pts == "360") ~  0,
                                                  (LC_aspect_raster_15_data_pts != "360")~ LC_aspect_raster_15_data_pts))


#SD

SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts = case_when((SD_aspect_raster_15_data_pts == "360") ~  0,
                                                  (SD_aspect_raster_15_data_pts != "360")~ SD_aspect_raster_15_data_pts))


# all points

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_8_categorical = case_when((all_points_aspect_raster_15_data_pts > 0 & all_points_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 337.5 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N", #359.99999
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 67.5 & all_points_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 112.5 & all_points_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                        (all_points_aspect_raster_15_data_pts >= 157.5 & all_points_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                        (all_points_aspect_raster_15_data_pts >= 202.5 & all_points_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                        (all_points_aspect_raster_15_data_pts >= 247.5 & all_points_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 292.5 & all_points_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees
all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts[584]
all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical[584]

# North, East, South, West

# the directions are a range of 90 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_4_categorical = case_when((all_points_aspect_raster_15_data_pts >= 0 & all_points_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

# all_points_fixed_field_data_processed_terrain_download <- write.csv(all_points_fixed_field_data_processed_terrain, "/Users/chewbecca/Morton Arboretum REU 2024/Untitled/QUBR_GenGeoEcoDemoCorr/data/all_points_fixed_field_data_processed_terrain.csv", row.names = F)
# View(all_points_fixed_field_data_processed_terrain)

# LM

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_8_categorical = case_when((LM_aspect_raster_15_data_pts > 0 & LM_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LM_aspect_raster_15_data_pts >= 337.5 & LM_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 67.5 & LM_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 112.5 & LM_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (LM_aspect_raster_15_data_pts >= 157.5 & LM_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (LM_aspect_raster_15_data_pts >= 202.5 & LM_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (LM_aspect_raster_15_data_pts >= 247.5 & LM_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 292.5 & LM_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# North, East, South, West

# the directions are a range of 90 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_4_categorical = case_when((LM_aspect_raster_15_data_pts >= 0 & LM_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LM_aspect_raster_15_data_pts >= 315 & LM_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LM_aspect_raster_15_data_pts >= 135 & LM_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LM_aspect_raster_15_data_pts >= 225 & LM_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


# LC

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_8_categorical = case_when((LC_aspect_raster_15_data_pts > 0 & LC_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 337.5 & LC_aspect_raster_15_data_pts < 359.999999999999999) ~ "N", 
                                                                (LC_aspect_raster_15_data_pts >= 22.5 & LC_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 67.5 & LC_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 112.5 & LC_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (LC_aspect_raster_15_data_pts >= 157.5 & LC_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (LC_aspect_raster_15_data_pts >= 202.5 & LC_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (LC_aspect_raster_15_data_pts >= 247.5 & LC_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 292.5 & LC_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees


# North, East, South, West

# the directions are a range of 90 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_4_categorical = case_when((LC_aspect_raster_15_data_pts >= 0 & LC_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 315 & LC_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LC_aspect_raster_15_data_pts >= 45 & LC_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LC_aspect_raster_15_data_pts >= 135 & LC_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LC_aspect_raster_15_data_pts >= 225 & LC_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

#SD
# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_8_categorical = case_when((SD_aspect_raster_15_data_pts > 0 & SD_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 337.5 & SD_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (SD_aspect_raster_15_data_pts >= 22.5 & SD_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 67.5 & SD_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 112.5 & SD_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (SD_aspect_raster_15_data_pts >= 157.5 & SD_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (SD_aspect_raster_15_data_pts >= 202.5 & SD_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (SD_aspect_raster_15_data_pts >= 247.5 & SD_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 292.5 & SD_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# North, East, South, West

# the directions are a range of 90 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_4_categorical = case_when((SD_aspect_raster_15_data_pts >= 0 & SD_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 315 & SD_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (SD_aspect_raster_15_data_pts >= 45 & SD_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (SD_aspect_raster_15_data_pts >= 135 & SD_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (SD_aspect_raster_15_data_pts >= 225 & SD_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


#creating X sequential columns in LC and SD point data which will make it easier to select random points from each grid later

#creating an x_sequential column that is 1 through the number of LM points, which will make it easier to randomly choose one point
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(LM_fixed_field_data_processed))

#creating an x_sequential column that is 1 through the number of LC points, which will make it easier to randomly choose one point
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(LC_fixed_field_data_processed))

#creating an x_sequential column that is 1 through the number of SD points, which will make it easier to randomly choose one point
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(SD_fixed_field_data_processed))


## combining elevation, slope, aspect data with distance data 

#LM
LM_fixed_field_data_processed_terrain_dist <- cbind(LM_fixed_field_data_processed_terrain, LM_fixed_field_data_processed_distance)

#LC
LC_fixed_field_data_processed_terrain_dist <- cbind(LC_fixed_field_data_processed_terrain, LC_fixed_field_data_processed_distance)

#SD
SD_fixed_field_data_processed_terrain_dist <- cbind(SD_fixed_field_data_processed_terrain, SD_fixed_field_data_processed_distance)


### LM ###

# removing the NAs from the explanatory variables, stored in a dataframe
LM_fixed_field_data_processed_terrain_dist_no_NA <- LM_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.FIXED)) %>% #Elevation NAs removed
  filter(!is.na(LM_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(LM_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed

## SCA ##

#removing the spatial geometry to be able to use the gam function
LM_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LM_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LM_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation has significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM with just elevation smoothing splines s()
LM_add.gam_SCA.terrain_dist.just.elevation.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants only s(elevation)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.just.elevation) # AICs
anova(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.just.elevation)  #ANOVA F-Test
#the more complex model with smoothing splines on all three quantitative variables is preferable

# setting up the dredge GAM of the larger equation, just quant variables
LM_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs

#comparing model with smoothing with and without aspect
AIC(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist)
summary(LM_add.gam_SCA.terrain_dist.dredge)
#aspect does not appear to be necessary

#making a gam with just elevation smoothed, no aspect
LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth)
dredge <- dredge(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
#just s(elevation is preferred)

#making a dredge model with just smooth elevation
LM_add.gam_SCA.terrain_dist.dredge.just.elev <- gam(Canopy_short ~ s(Elevation..m.FIXED), 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth)

#comparing the models with smoothed distance, elevation, and slope to one with just smoothed elevation and the other variables, and one with only s(elevation)
AIC(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth, LM_add.gam_SCA.terrain_dist.dredge.just.elev)
#LM_add.gam_SCA.terrain_dist.dredge has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_SCA.terrain_dist.dredge 
summary(LM_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LM_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LM_add.gam_SCA.terrain_dist.dredge

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_SCA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LM_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LM_add.gam_SCA.inter, LM_add.gam_SCA.terrain_dist.dredge)

#overall best model:LM_add.gam_SCA.terrain_dist.dredge



## LCA ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LM_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, elevation, and slope have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_LCA.terrain_dist.no.aspect <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_LCA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation
LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_long ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                    data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth)

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect, LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_LCA.terrain_dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_LCA.terrain_dist.dredge 
summary(LM_add.gam_LCA.terrain_dist)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.terrain_dist) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LM_add.gam_LCA.terrain_dist)

#Chosen model: LM_add.gam_LCA.terrain_dist

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_LCA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LM_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LM_add.gam_LCA.inter, LM_add.gam_LCA.terrain_dist)

#overall best model: LM_add.gam_LCA.terrain_dist


## CA ##

# I am using the logged transformation of canop area to get more normal residuals and less heteroscedasticity

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LM_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_CA.terrain_dist.no.aspect <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                             data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_CA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_area_lg ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect, LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_CA.terrain_dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CA.terrain_dist.no.aspect 
summary(LM_add.gam_CA.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.terrain_dist.no.aspect) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LM_add.gam_CA.terrain_dist.no.aspect)

#Chosen model: LM_add.gam_CA.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_CA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LM_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_CA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LM_add.gam_CA.inter, LM_add.gam_CA.terrain_dist.no.aspect)

#overall best model: LM_add.gam_CA.terrain_dist.no.aspect


## CS ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LM_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, slope, elevation have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_CS.terrain_dist.no.aspect <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_CS.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth <- gam(Crown_spread ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect, LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_CS.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CS.terrain_dist.no.aspect 
summary(LM_add.gam_CS.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.terrain_dist.no.aspect) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LM_add.gam_CS.terrain_dist.no.aspect)

#Chosen model: LM_add.gam_CS.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_CS.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LM_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LM_add.gam_CS.inter, LM_add.gam_CS.terrain_dist.no.aspect)

#overall best model: LM_add.gam_CS.terrain_dist.no.aspect


## DBH ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LM_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only elevation have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_DBH.terrain_dist.no.aspect <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LM_add.gam_DBH.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth <- gam(DBH_ag ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit

#making a model with no smoothing splines
LM_add.gam_DBH.terrain_dist.dredge.no.smooth <- gam(DBH_ag ~ d + Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit


#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect, LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth, LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#LM_add.gam_DBH.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.terrain_dist.dredge.no.smooth) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)

#Chosen model: LM_add.gam_DBH.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_DBH.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_DBH.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LM_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LM_add.gam_DBH.inter, LM_add.gam_DBH.terrain_dist)
#while there is no significant interaction, the model with the interaction performs significantly better

#overall best model: LM_add.gam_DBH.terrain_dist.dredge.no.smooth (a multiple linear regression)


### LC ###

# removing the NAs from the explanatory variables, stored in a dataframe
LC_fixed_field_data_processed_terrain_dist_no_NA <- LC_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.FIXED)) %>% #Elevation NAs removed
  filter(!is.na(LC_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(LC_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed

## SCA ##

#removing the spatial geometry to be able to use the gam function
LC_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LC_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LC_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#none of the variables has significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation) and aspect

# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.just.elevation.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.just.elevation.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_SCA.terrain_dist.just.elevation.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants only s(elevation) and aspect

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth) # AICs
anova(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth)  #ANOVA F-Test
#the more complex model with smoothing splines on all three quantitative variables is preferable

# setting up the dredge model aspect and s(elevation)
LC_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs

#comparing model with smoothing with and without aspect
AIC(LC_add.gam_SCA.terrain_dist.dredge, LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth)
#LC_add.gam_SCA.terrain_dist.dredge has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_SCA.terrain_dist.dredge 
summary(LC_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LC_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LC_add.gam_SCA.terrain_dist.dredge


#plotting the model LC_add.gam_SCA.terrain_dist.dredge
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
visreg(LC_add.gam_SCA.terrain_dist.dredge, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA")

#plotting the model LC_add.gam_SCA.terrain_dist
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_SCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_SCA.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LC_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_SCA.inter)
#there was a significant interaction term 

#interaction plots
plot.gam(LC_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LC_add.gam_SCA.inter, LC_add.gam_SCA.terrain_dist.dredge)

#overall best model: LC_add.gam_SCA.inter


## LCA ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LC_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, elevation, and slope do not have significant non-linear function, aspect is significant linearly for W

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation)

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.no.aspect.slope <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED), 
                                             data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.no.aspect.slope) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_LCA.terrain_dist.no.aspect.slope) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope) # AICs
anova(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope)  #ANOVA F-Test
#the model with aspect and slope is significantly better

#making a dredge model with just smooth elevation
LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth)

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope, LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth)
#none of the models have much better AIC values 

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_LCA.terrain_dist.dredge 
summary(LC_add.gam_LCA.terrain_dist)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.terrain_dist) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LC_add.gam_LCA.terrain_dist)

#Chosen model: LC_add.gam_LCA.terrain_dist

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_LCA.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LC_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LC_add.gam_LCA.inter, LC_add.gam_LCA.terrain_dist)

#overall best model: LC_add.gam_LCA.terrain_dist


## CA ##

# I am using the logged transformation of canop area to get more normal residuals and less heteroscedasticity

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LC_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

 #elevation have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), aspect

# Checking a GAM without s(slope) but slope instead
LC_add.gam_CA.terrain_dist.less.s <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.less.s) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_CA.terrain_dist.less.s) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), aspect

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s) # AICs
anova(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s)  #ANOVA F-Test
#the model with s(slope) is significantly better

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.no.slope <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s, LC_add.gam_CA.terrain_dist.dredge.no.slope)
#LC_add.gam_CA.terrain_dist.dredge.no.slope has lowest AIC, but really none are much better

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CA.terrain_dist.dredge.no.slope 
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)
#but all models seem to do similarly well

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.terrain_dist.dredge.no.slope) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)

#Chosen model: LC_add.gam_CA.terrain_dist.dredge.no.slope

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.terrain_dist.dredge.no.slope, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CA.terrain_dist.dredge.no.slope, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_CA.terrain_dist.dredge.no.slope, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LC_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CA.inter)
#there was is a significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LC_add.gam_CA.inter, LC_add.gam_CA.terrain_dist.no.aspect)

#overall best model: LC_add.gam_CA.terrain_dist.no.aspect


## CS ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LC_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, slope, elevation have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LC_add.gam_CS.terrain_dist.no.aspect <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_CS.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect) # AICs
anova(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth <- gam(Crown_spread ~ d + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts, 
                                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect, LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#LC_add.gam_CS.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CS.terrain_dist.no.aspect 
summary(LC_add.gam_CS.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.terrain_dist.no.aspect) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LC_add.gam_CS.terrain_dist.no.aspect)

#Chosen model: LC_add.gam_CS.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_CS.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LC_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LC_add.gam_CS.inter, LC_add.gam_CS.terrain_dist.no.aspect)

#overall best model: LC_add.gam_CS.terrain_dist.no.aspect


## DBH ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
LC_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#none have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d)

# Checking a GAM without aspect
LC_add.gam_DBH.terrain_dist.only.distance <- gam(DBH_ag ~ s(d), 
                                             data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.only.distance) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(LC_add.gam_DBH.terrain_dist.only.distance) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance) # AICs
anova(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance)  #ANOVA F-Test
#the model with only distance is sufficient

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth <- gam(DBH_ag ~ s(d) + Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth)
#s(distance) is marginally not significantly  non-linear, slope and elevation are not significantly useful with the linear fit

#making a model with no smoothing splines
LC_add.gam_DBH.terrain_dist.dredge.no.smooth <- gam(DBH_ag ~ d + Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.no.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit


#comparing the models 
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance, LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth, LC_add.gam_DBH.terrain_dist.dredge.no.smooth)
#LC_add.gam_DBH.terrain_dist.dredge.no.smooth has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LC_add.gam_DBH.terrain_dist)
#the model with just distance smooth is marginally better, so I will use the one with smoothing splines on all of the quantitative variables

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.terrain_dist) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(LC_add.gam_DBH.terrain_dist)

#Chosen model: LC_add.gam_DBH.terrain_dist

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_DBH.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_DBH.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
LC_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(LC_add.gam_DBH.inter, LC_add.gam_DBH.terrain_dist)
#while there is significant interaction, the model with the interaction performs significantly better

#overall best model: LC_add.gam_DBH.inter



### SD ###

# removing the NAs from the explanatory variables, stored in a dataframe
SD_fixed_field_data_processed_terrain_dist_no_NA <- SD_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.)) %>% #Elevation NAs removed
  filter(!is.na(SD_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(SD_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) %>% #DBH NAs removed
  mutate(Elevation..m. = as.numeric(Elevation..m.))

## SCA ##

#removing the spatial geometry to be able to use the gam function
SD_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(SD_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
SD_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values


#slope has significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM with no aspect
SD_add.gam_SCA.terrain_dist.no.aspect <- gam(Canopy_short ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                         data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#slope of the linear fits are significant

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_SCA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect) # AICs
anova(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model without aspect is better

# setting up the dredge GAM with only slope being smoothed
SD_add.gam_SCA.terrain_dist.slope.smooth <- gam(Canopy_short ~ d + Elevation..m. + s(SD_slope_raster_15_data_pts), 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_SCA.terrain_dist.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants elevation, s(slope)

# setting up the dredge GAM with only slope being smoothed, no distance
SD_add.gam_SCA.terrain_dist.slope.smooth.dredge <- gam(Canopy_short ~ Elevation..m. + s(SD_slope_raster_15_data_pts), 
                                                data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.slope.smooth.dredge)
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_SCA.terrain_dist.slope.smooth.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)


#comparing models 
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect,  SD_add.gam_SCA.terrain_dist.slope.smooth, SD_add.gam_SCA.terrain_dist.slope.smooth.dredge)
summary(SD_add.gam_SCA.terrain_dist.no.aspect)
#SD_add.gam_SCA.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_SCA.terrain_dist.dredge 
summary(SD_add.gam_SCA.terrain_dist.no.aspect)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.terrain_dist.no.aspect) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(SD_add.gam_SCA.terrain_dist.no.aspect)

#Chosen model: SD_add.gam_SCA.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_SCA.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
SD_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(SD_add.gam_SCA.inter, SD_add.gam_SCA.terrain_dist.no.aspect)

#overall best model:SD_add.gam_SCA.terrain_dist.no.aspect


## LCA ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
SD_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function and sw was significant aspect

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(slope) and aspect

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.slope.smooth <- gam(Canopy_long ~ d + Elevation..m. + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.slope.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only west is significant
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_LCA.terrain_dist.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#making a dredge model with just smooth elevation
SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist <- gam(Canopy_long ~ s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)
#s(slope) and NE and SW are significant

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(SD_add.gam_LCA.terrain_dist, SD_add.gam_LCA.terrain_dist.slope.smooth, SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)
#SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist 
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#Chosen model: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist, select=1, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
SD_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(SD_add.gam_LCA.inter, SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#overall best model: SD_add.gam_LCA.inter


## CA ##

# I am using the logged transformation of canop area to get more normal residuals and less heteroscedasticity

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
SD_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
SD_add.gam_CA.terrain_dist.no.aspect <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_CA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect) # AICs
anova(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and slope and no aspect
SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth <- gam(Canopy_area_lg ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#making a dredge model with just smooth elevation and slope
SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist <- gam(Canopy_area_lg ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                                data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit


#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect, SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth, SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist)
#SD_add.gam_CA.terrain_dist.no.aspect  has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CA.terrain_dist.no.aspect 
summary(SD_add.gam_CA.terrain_dist.no.aspect)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.terrain_dist.no.aspect) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(SD_add.gam_CA.terrain_dist.no.aspect)

#Chosen model: SD_add.gam_CA.terrain_dist.no.aspect

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_CA.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
SD_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(SD_add.gam_CA.inter, SD_add.gam_CA.terrain_dist.no.aspect)

#overall best model: SD_add.gam_CA.terrain_dist.no.aspect


## CS ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
SD_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function and Southwest are significant

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope), and aspect

# Checking a GAM without d smoothed
SD_add.gam_CS.terrain_dist.no.dist.smooth <- gam(Crown_spread ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.no.dist.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation and slope) are significant
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_CS.terrain_dist.no.dist.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope), aspect

#making a dredge model with just smooth elevation and no aspect
SD_add.gam_CS.terrain_dist.no.dist <- gam(Crown_spread ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.no.dist)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models
AIC(SD_add.gam_CS.terrain_dist, SD_add.gam_CS.terrain_dist.no.dist.smooth, SD_add.gam_CS.terrain_dist.no.dist)
#SD_add.gam_CS.terrain_dist.no.dist has lowest AIC, but none are significantlt better than the others

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CS.terrain_dist.no.aspect 
summary(SD_add.gam_CS.terrain_dist.no.dist)
#but all models seem to do similarly well

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.terrain_dist.no.dist) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(SD_add.gam_CS.terrain_dist.no.dist)

#Chosen model: SD_add.gam_CS.terrain_dist.no.dist

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CS.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_CS.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CS.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_CS.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
SD_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(SD_add.gam_CS.inter, SD_add.gam_CS.terrain_dist)

#overall best model: SD_add.gam_CS.terrain_dist


## DBH ##

# Checking a GAM with smoothing splines s(), can't put splines on a categorical variable
SD_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only slope have significant non-linear function 

#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

# Checking a GAM without smoothing distance 
SD_add.gam_DBH.terrain_dist.no.dist.smooth <- gam(DBH_ag ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.no.dist.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_DBH.terrain_dist.no.dist.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#making a dredge model with just smooth elevation and no aspect
SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth <- gam(DBH_ag ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#s(slope) is significantly non-linear but not s(elevation), slope and distance are not significantly useful with the linear fit
#dredging the gam to see which variables to keep
dredge <- dredge(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the models 
AIC(SD_add.gam_DBH.terrain_dist, SD_add.gam_DBH.terrain_dist.no.dist.smooth, SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(SD_add.gam_DBH.terrain_dist.dredge.no.smooth)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth) #pretty normal residuals and no heterodescadisticty 

#looking at significance
summary(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)

#Chosen model: SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth

#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth, select=1, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_DBH.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#looking for interaction using tensor interaction to get interaction smooths
SD_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

AIC(SD_add.gam_DBH.inter, SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#the interaction is not significant

#overall best model: SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth 




