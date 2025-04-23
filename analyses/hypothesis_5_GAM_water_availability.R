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






## combining elevation, slope, aspect data with distance data 

#LM
LM_fixed_field_data_processed_terrain_dist <- cbind(LM_fixed_field_data_processed_terrain, LM_fixed_field_data_processed_distance)

#LC
LC_fixed_field_data_processed_terrain_dist <- cbind(LC_fixed_field_data_processed_terrain, LC_fixed_field_data_processed_distance)

#SD
SD_fixed_field_data_processed_terrain_dist <- cbind(SD_fixed_field_data_processed_terrain, SD_fixed_field_data_processed_distance)




