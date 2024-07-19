#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(perm.t.test) #permutation t test 


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

#plotting the tree points by elevation (m)
ggplot()+
  geom_sf(data = LM_fixed_field_data_processed, aes(color = Elevation..m.FIXED))


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

#creating a new elevation column so the values that were mistakenly put in feet are in meters
LM_fixed_field_data_processed <-  LM_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))

##creating a new elevation column so that the 360 m outlier is 460
SD_fixed_field_data_processed <-  SD_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. == 360) ~ 460, 
                                        (Elevation..m. != 360) ~ Elevation..m.))
View(SD_fixed_field_data_processed)
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

#creating the aspect and slope rasters

#loading in xyz ASCII from INEGI on elevation, 5 m resolution rasters of elevation around BCS
continental_relief_elevation_xyz <- read.table("./data/ASCII Elevation Inegi/f12b43b4_ms.xyz")
continental_relief_elevation_raster <- rasterFromXYZ(continental_relief_elevation_xyz) #26912
plot(continental_relief_elevation_raster)

elevation_xyz_f1 <- read.table("./data/ASCII Elevation Inegi/f12b14f1_ms.xyz")
elevation_xyz_f1 <- rasterFromXYZ(elevation_xyz_f1) #26912
plot(elevation_xyz_f1)

elevation_xyz_f2 <- read.table("./data/ASCII Elevation Inegi/f12b14f2_ms.xyz")
elevation_xyz_f2 <- rasterFromXYZ(elevation_xyz_f2) #26912
plot(elevation_xyz_f2)

elevation_xyz_f3 <- read.table("./data/ASCII Elevation Inegi/f12b14f3_ms.xyz")
elevation_xyz_f3 <- rasterFromXYZ(elevation_xyz_f3) #26912
plot(elevation_xyz_f3)

elevation_xyz_f4 <- read.table("./data/ASCII Elevation Inegi/f12b14f4_ms.xyz")
elevation_xyz_f4 <- rasterFromXYZ(elevation_xyz_f4) #26912
plot(elevation_xyz_f4)

elevation_xyz_c4 <- read.table("./data/ASCII Elevation Inegi/f12b24c4_ms.xyz")
elevation_xyz_c4 <- rasterFromXYZ(elevation_xyz_c4) #26912
plot(elevation_xyz_c4)

elevation_xyz_b25a3 <- read.table("./data/ASCII Elevation Inegi/f12b25a3_ms.xyz")
elevation_xyz_b25a3 <- rasterFromXYZ(elevation_xyz_b25a3) #26912
plot(elevation_xyz_b25a3)

elevation_xyz_b34a1 <- read.table("./data/ASCII Elevation Inegi/f12b34a1_ms.xyz")
elevation_xyz_b34a1 <- rasterFromXYZ(elevation_xyz_b34a1) #26912
plot(elevation_xyz_b34a1)

elevation_xyz_b34a2 <- read.table("./data/ASCII Elevation Inegi/f12b34a2_ms.xyz")
elevation_xyz_b34a2 <- rasterFromXYZ(elevation_xyz_b34a2) #26912
plot(elevation_xyz_b34a2)

elevation_xyz_b34a3 <- read.table("./data/ASCII Elevation Inegi/f12b34a3_ms.xyz")
elevation_xyz_b34a3 <- rasterFromXYZ(elevation_xyz_b34a3) #26912
plot(elevation_xyz_b34a3)

elevation_xyz_b34a4 <- read.table("./data/ASCII Elevation Inegi/f12b34a4_ms.xyz")
elevation_xyz_b34a4 <- rasterFromXYZ(elevation_xyz_b34a4) #26912
plot(elevation_xyz_b34a4)

elevation_xyz_b34b1 <- read.table("./data/ASCII Elevation Inegi/f12b34b1_ms.xyz")
elevation_xyz_b34b1 <- rasterFromXYZ(elevation_xyz_b34b1) #26912
plot(elevation_xyz_b34b1)

elevation_xyz_b34b2 <- read.table("./data/ASCII Elevation Inegi/f12b34b2_ms.xyz")
elevation_xyz_b34b2 <- rasterFromXYZ(elevation_xyz_b34b2) #26912
plot(elevation_xyz_b34b2)

elevation_xyz_b34b3 <- read.table("./data/ASCII Elevation Inegi/f12b34b3_ms.xyz")
elevation_xyz_b34b3 <- rasterFromXYZ(elevation_xyz_b34b3) #26912
plot(elevation_xyz_b34b3)

elevation_xyz_b34b4 <- read.table("./data/ASCII Elevation Inegi/f12b34b4_ms.xyz")
elevation_xyz_b34b4 <- rasterFromXYZ(elevation_xyz_b34b4) #26912
plot(elevation_xyz_b34b4)

elevation_xyz_b34c1 <- read.table("./data/ASCII Elevation Inegi/f12b34c1_ms.xyz")
elevation_xyz_b34c1 <- rasterFromXYZ(elevation_xyz_b34c1) #26912
plot(elevation_xyz_b34c1)

elevation_xyz_b34c3 <- read.table("./data/ASCII Elevation Inegi/f12b34c3_ms.xyz")
elevation_xyz_b34c3 <- rasterFromXYZ(elevation_xyz_b34c3) #26912
plot(elevation_xyz_b34c3)

#merging all of the rasters
f1_f4_merged_rasters <- raster::merge(elevation_xyz_f1, elevation_xyz_f2, 
                                      elevation_xyz_f3,  elevation_xyz_f4, 
                                      elevation_xyz_c4,elevation_xyz_b25a3,
                                      elevation_xyz_b34a1, elevation_xyz_b34a2,
                                      elevation_xyz_b34a3, elevation_xyz_b34a4, 
                                      elevation_xyz_b34b1, elevation_xyz_b34b2,
                                      elevation_xyz_b34b3, elevation_xyz_b34b4,
                                      elevation_xyz_b34c1, elevation_xyz_b34c3)
plot(f1_f4_merged_rasters)

#plotting the merged rasters
ggplot()+
  geom_raster(data= as.data.frame(f1_f4_merged_rasters, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = fixed_field_data_processed_NN_UTM)

#mapping cropped 
f1_f4_merged_rasters_cropped <- crop(f1_f4_merged_rasters, extent(c(SD_box[1], SD_box[3], SD_box[2], SD_box[4])))

#plotting the cropped rasters, cropped to cover SD
ggplot()+
  geom_raster(data= as.data.frame(f1_f4_merged_rasters, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = SD_fixed_field_data_processed)

#isolating just the SD 5 m portions
SD_rasters <- raster::merge(elevation_xyz_b34a1, elevation_xyz_b34a2,
                             elevation_xyz_b34a3, elevation_xyz_b34a4, 
                             elevation_xyz_b34b1, elevation_xyz_b34b2,
                             elevation_xyz_b34b3, elevation_xyz_b34b4,
                             elevation_xyz_b34c1, elevation_xyz_b34c3)

plot(SD_rasters)

#very close with the elevation rasters
ggplot()+
  geom_raster(data= as.data.frame(SD_rasters, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = SD_fixed_field_data_processed)

#loading in 50 M resolution elevation rasters for LM and LC
f12b12 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b12me.bil")
plot(f12b12)

f12b13 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b13me.bil")
plot(f12b13)

f12b14 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b14me.bil")
plot(f12b14)

f12b22 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b22me.bil")
plot(f12b22)

f12b23 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b23me.bil")
plot(f12b23)

f12b33 <- raster("./data/ASCII Elevation Inegi/smaller scale elevation/f12b33me.bil")
plot(f12b33)

f12b_merged_rasters <- raster::merge(f12b12, f12b13, 
                                     f12b14,  f12b22, 
                                     f12b23, f12b33)
plot(f12b_merged_rasters)

#projecting the 50 m resolution raster to be in UTM12N 
crs_utm <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
proj4string(f12b_merged_rasters) <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#mapping with cropped for LM 
f12b_merged_rasters_cropped_LM <- crop(f12b_merged_rasters, extent(c(LM_box[1]-200, LM_box[3]+200, LM_box[2]-200, LM_box[4]+200)))

ggplot()+
  geom_raster(data= as.data.frame(f12b_merged_rasters_cropped_LM, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = LM_fixed_field_data_processed)

#mapping with cropped for LC
f12b_merged_rasters_cropped_LC <- crop(f12b_merged_rasters, extent(c(LC_box[1]-200, LC_box[3]+200, LC_box[2]-200, LC_box[4]+200)))

ggplot()+
  geom_raster(data= as.data.frame(f12b_merged_rasters_cropped_LC, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = LC_fixed_field_data_processed)

#mapping with cropped for SD
f12b_merged_rasters_cropped_SD <- crop(f12b_merged_rasters, extent(c(SD_box[1], SD_box[3], SD_box[2], SD_box[4])))
ggplot()+
  geom_raster(data= as.data.frame(f12b_merged_rasters_cropped_SD, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = SD_fixed_field_data_processed)

f1202 <- raster("./data/ASCII Elevation Inegi/larger scale elevation/f1202mde.bil")
plot(f1202)

f1206 <- raster("./data/ASCII Elevation Inegi/larger scale elevation/f1206mde.bil")
plot(f1206)


#50 m resolution elevation for SD
f12B24 <- raster("./data/ASCII Elevation Inegi/larger scale elevation/f12b24me.bil")
plot(f12B24)

crs_utm <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"
proj4string(f12B24) <- "+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs"

#cropping the F12B24 raster to be the extent of SD
f12B24_cropped_SD <- crop(f12B24, extent(c(SD_box[1]-50, SD_box[3]+50, SD_box[2]-50, SD_box[4]+50)))

#Plotting the cropped raster
ggplot()+
  geom_raster(data= as.data.frame(f12B24_cropped_SD, xy = T), aes(x=x, y=y, fill = f12b24me))+
  geom_sf(data = SD_fixed_field_data_processed)

#contour lines map
F12B24_tif <- raster(paste0("./data/ASCII Elevation Inegi/Larger Areas/F12B24/702825702502_t/f12b24.tif"))
plot(F12B24_tif)

F12B24_tif_cropped_SD <- crop(F12B24_tif, extent(c(SD_box[1], SD_box[3], SD_box[2], SD_box[4])))

ggplot()+
  geom_raster(data= as.data.frame(F12B24_tif_cropped_SD, xy = T), aes(x=x, y=y, fill = f12b24))+
  geom_sf(data = SD_fixed_field_data_processed)


## Extracting the slope 

#LM

#extracting the slope in degrees, using the queens method (neighbor = 8)
LM_slope_raster_50 <- terra::terrain(f12b_merged_rasters_cropped_LM, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_slope_raster_50, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LM_fixed_field_data_processed)

#LC

#extracting the slope in degrees, using the queens method (neighbor = 8)
LC_slope_raster_50 <- terra::terrain(f12b_merged_rasters_cropped_LC, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_slope_raster_50, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LC_fixed_field_data_processed)

#SD

#extracting the slope in degrees, using the queens method (neighbor = 8)
SD_slope_raster_50 <- terra::terrain(f12B24_cropped_SD, unit = 'degrees', neighbors = 8, 'slope')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_slope_raster_50, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = SD_fixed_field_data_processed)

## Extracting the aspect 

#LM

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LM_aspect_raster_50 <- terra::terrain(f12b_merged_rasters_cropped_LM, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_aspect_raster_50, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LM_fixed_field_data_processed)

#LC

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LC_aspect_raster_50 <- terra::terrain(f12b_merged_rasters_cropped_LC, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_aspect_raster_50, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LC_fixed_field_data_processed)

#SD

#extracting the aspect in degrees, using the queens method (neighbor = 8)
SD_aspect_raster_50 <- terra::terrain(f12B24_cropped_SD, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_aspect_raster_50, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = SD_fixed_field_data_processed)

#creating dataframes for each population and the slope and aspect data

#LM

LM_aspect_raster_50_data_pts <- extract(LM_aspect_raster_50, LM_fixed_field_data_processed) #extracting aspect for each point value
LM_slope_raster_50_data_pts <- extract(LM_slope_raster_50, LM_fixed_field_data_processed) #extracting slope for each point value
LM_fixed_field_data_processed_terrain <- cbind(LM_fixed_field_data_processed, LM_aspect_raster_50_data_pts) #bind the aspect data for each point to the LM point dataframe
LM_fixed_field_data_processed_terrain <- cbind(LM_fixed_field_data_processed_terrain, LM_slope_raster_50_data_pts) #bind the slope data for each point to the LM point dataframe

View(LM_fixed_field_data_processed_terrain)


#LC
LC_aspect_raster_50_data_pts <- extract(LC_aspect_raster_50, LC_fixed_field_data_processed) #extracting aspect for each point value
LC_slope_raster_50_data_pts <- extract(LC_slope_raster_50, LC_fixed_field_data_processed) #extracting slope for each point value
LC_fixed_field_data_processed_terrain <- cbind(LC_fixed_field_data_processed, LC_aspect_raster_50_data_pts) #bind the aspect data for each point to the SD point dataframe

LC_fixed_field_data_processed_terrain <- cbind(LC_fixed_field_data_processed_terrain, LC_slope_raster_50_data_pts) #bind the slope data for each point to the SD point dataframe

View(LC_fixed_field_data_processed_terrain)


#SD
SD_aspect_raster_50_data_pts <- extract(SD_aspect_raster_50, SD_fixed_field_data_processed) #extracting aspect for each point value
SD_slope_raster_50_data_pts <- extract(SD_slope_raster_50, SD_fixed_field_data_processed) #extracting slope for each point value
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed, SD_aspect_raster_50_data_pts) #bind the aspect data for each point to the SD point dataframe
SD_fixed_field_data_processed_terrain_otherversion <- SD_fixed_field_data_processed %>%  
  add_column(SD_aspect_raster_50_data_pts)
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed_terrain, SD_slope_raster_50_data_pts) #bind the slope data for each point to the SD point dataframe

View(SD_fixed_field_data_processed_terrain)
View(SD_fixed_field_data_processed_terrain_otherversion)

#recategorizing the aspect data

#setting values of 360 to 0 

#LM
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_50_data_pts = case_when((LM_aspect_raster_50_data_pts == "360") ~  0,
                                                  (LM_aspect_raster_50_data_pts != "360")~ LM_aspect_raster_50_data_pts))
View(LM_fixed_field_data_processed_terrain)

#LC

LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_50_data_pts = case_when((LC_aspect_raster_50_data_pts == "360") ~  0,
                                                  (LC_aspect_raster_50_data_pts != "360")~ LC_aspect_raster_50_data_pts))
View(LC_fixed_field_data_processed_terrain)

#SD

SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_50_data_pts = case_when((SD_aspect_raster_50_data_pts == "360") ~  0,
                                                  (SD_aspect_raster_50_data_pts != "360")~ SD_aspect_raster_50_data_pts))
View(SD_fixed_field_data_processed_terrain)

# LM

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_50_data_pts_8_categorical = case_when((LM_aspect_raster_50_data_pts > 0 & LM_aspect_raster_50_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                              (LM_aspect_raster_50_data_pts >= 337.5 & LM_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                              (LM_aspect_raster_50_data_pts >= 22.5 & LM_aspect_raster_50_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                              (LM_aspect_raster_50_data_pts >= 67.5 & LM_aspect_raster_50_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                              (LM_aspect_raster_50_data_pts >= 112.5 & LM_aspect_raster_50_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                              (LM_aspect_raster_50_data_pts >= 157.5 & LM_aspect_raster_50_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                              (LM_aspect_raster_50_data_pts >= 202.5 & LM_aspect_raster_50_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                              (LM_aspect_raster_50_data_pts >= 247.5 & LM_aspect_raster_50_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                              (LM_aspect_raster_50_data_pts >= 292.5 & LM_aspect_raster_50_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees
  
View(LM_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_50_data_pts_4_categorical = case_when((LM_aspect_raster_50_data_pts >= 0 & LM_aspect_raster_50_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                              (LM_aspect_raster_50_data_pts >= 315 & LM_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                              (LM_aspect_raster_50_data_pts >= 22.5 & LM_aspect_raster_50_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                              (LM_aspect_raster_50_data_pts >= 135 & LM_aspect_raster_50_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                              (LM_aspect_raster_50_data_pts >= 225 & LM_aspect_raster_50_data_pts < 315) ~ "W")) #west is between 225 and 315
                                                     
View(LM_fixed_field_data_processed_terrain)

# LC

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_50_data_pts_8_categorical = case_when((LC_aspect_raster_50_data_pts > 0 & LC_aspect_raster_50_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LC_aspect_raster_50_data_pts >= 337.5 & LC_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                                (LC_aspect_raster_50_data_pts >= 22.5 & LC_aspect_raster_50_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (LC_aspect_raster_50_data_pts >= 67.5 & LC_aspect_raster_50_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (LC_aspect_raster_50_data_pts >= 112.5 & LC_aspect_raster_50_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (LC_aspect_raster_50_data_pts >= 157.5 & LC_aspect_raster_50_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (LC_aspect_raster_50_data_pts >= 202.5 & LC_aspect_raster_50_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (LC_aspect_raster_50_data_pts >= 247.5 & LC_aspect_raster_50_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (LC_aspect_raster_50_data_pts >= 292.5 & LC_aspect_raster_50_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

View(LC_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_50_data_pts_4_categorical = case_when((LC_aspect_raster_50_data_pts >= 0 & LC_aspect_raster_50_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LC_aspect_raster_50_data_pts >= 315 & LC_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                                (LC_aspect_raster_50_data_pts >= 45 & LC_aspect_raster_50_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LC_aspect_raster_50_data_pts >= 135 & LC_aspect_raster_50_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LC_aspect_raster_50_data_pts >= 225 & LC_aspect_raster_50_data_pts < 315) ~ "W")) #west is between 225 and 315

View(LC_fixed_field_data_processed_terrain)


#SD
# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_50_data_pts_8_categorical = case_when((SD_aspect_raster_50_data_pts > 0 & SD_aspect_raster_50_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (SD_aspect_raster_50_data_pts >= 337.5 & SD_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                                (SD_aspect_raster_50_data_pts >= 22.5 & SD_aspect_raster_50_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (SD_aspect_raster_50_data_pts >= 67.5 & SD_aspect_raster_50_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (SD_aspect_raster_50_data_pts >= 112.5 & SD_aspect_raster_50_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (SD_aspect_raster_50_data_pts >= 157.5 & SD_aspect_raster_50_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (SD_aspect_raster_50_data_pts >= 202.5 & SD_aspect_raster_50_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (SD_aspect_raster_50_data_pts >= 247.5 & SD_aspect_raster_50_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (SD_aspect_raster_50_data_pts >= 292.5 & SD_aspect_raster_50_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

View(SD_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_50_data_pts_4_categorical = case_when((SD_aspect_raster_50_data_pts >= 0 & SD_aspect_raster_50_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (SD_aspect_raster_50_data_pts >= 315 & SD_aspect_raster_50_data_pts < 359.99999) ~ "N",
                                                                (SD_aspect_raster_50_data_pts >= 45 & SD_aspect_raster_50_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (SD_aspect_raster_50_data_pts >= 135 & SD_aspect_raster_50_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (SD_aspect_raster_50_data_pts >= 225 & SD_aspect_raster_50_data_pts < 315) ~ "W")) #west is between 225 and 315

View(SD_fixed_field_data_processed_terrain)



#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")


ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Elevation..m.FIXED))+
  xlab("Elevation (m)")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation (m)")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation (m)")+
  ylab("Frequency")

#histograms for slope

ggplot(LM_fixed_field_data_processed_terrain) + # Generate the base plot
  geom_histogram(aes(x = LM_slope_raster_50_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + # Generate the base plot
  geom_histogram(aes(x = LC_slope_raster_50_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + # Generate the base plot
  geom_histogram(aes(x = SD_slope_raster_50_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

#barcharts for aspect

# 8 categories of direction

ggplot(LM_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LM_aspect_raster_50_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LC_aspect_raster_50_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = SD_aspect_raster_50_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

# 4 categories of direction

ggplot(LM_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LM_aspect_raster_50_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LC_aspect_raster_50_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = SD_aspect_raster_50_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")



#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)


### Sizes vs. Elevation ###

# LM 

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_elev, aes(x= LM_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_elev, aes(sample = LM_lm_sca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_sca_elev, aes(x = LM_lm_sca_elev$fitted.values, y = LM_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_sca_elev)

#correlation test
cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression

LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_lca_elev, aes(x= LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_lca_elev, aes(sample = LM_lm_lca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_lca_elev, aes(x = LM_lm_lca_elev$fitted.values, y = LM_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CA_elev, aes(x= LM_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CA_elev, aes(sample = LM_lm_CA_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CA_elev)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CS_elev, aes(x= LM_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CS_elev, aes(sample = LM_lm_CS_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CS_elev, aes(x = LM_lm_CS_elev$fitted.values, y = LM_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CS_elev)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_DBH_elev, aes(x= LM_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_DBH_elev, aes(sample = LM_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_DBH_elev, aes(x = LM_lm_DBH_elev$fitted.values, y = LM_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_DBH_elev)



#LC linear models

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with logged transformation of short canopy axis
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short_lg ~ LC_fixed_field_data_processed$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_sca_elev, aes(x= LC_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_elev, aes(sample = LC_lm_sca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_sca_elev, aes(x = LC_lm_sca_elev$fitted.values, y = LC_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(LC_lm_sca_elev)

#correlation test
cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_long_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long ~ LC_fixed_field_data_processed$Elevation..m.)

#linear transformation with logged long canopy axis
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long_lg ~ LC_fixed_field_data_processed$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_lca_elev, aes(x= LC_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_lca_elev, aes(sample = LC_lm_lca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_lca_elev, aes(x = LC_lm_lca_elev$fitted.values, y = LC_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(LC_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of canopy area
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of canopy area
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CA_elev, aes(x= LC_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CA_elev, aes(sample = LC_lm_CA_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(LC_lm_CA_elev)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread ~ LC_fixed_field_data_processed$Elevation..m.)

#linear transformation with logged crown spread
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear transformation with square rooted crown spread
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CS_elev, aes(x= LC_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CS_elev, aes(sample = LC_lm_CS_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_CS_elev, aes(x = LC_lm_CS_elev$fitted.values, y = LC_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(LC_lm_CS_elev)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with logged transformation of aggregated DBH
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of aggregated DBH
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_DBH_elev, aes(x= LC_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_DBH_elev, aes(sample = LC_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_DBH_elev, aes(x = LC_lm_DBH_elev$fitted.values, y = LC_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(LC_lm_DBH_elev)



#SD linear models

#using the fixed elevation 

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#creating the linear regression
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_sca_elev, aes(x= SD_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_elev, aes(sample = SD_lm_sca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_sca_elev, aes(x = SD_lm_sca_elev$fitted.values, y = SD_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_sca_elev)

#correlation test
cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_lca_elev, aes(x= SD_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_lca_elev, aes(sample = SD_lm_lca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_lca_elev, aes(x = SD_lm_lca_elev$fitted.values, y = SD_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CA_elev, aes(x= SD_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CA_elev, aes(sample = SD_lm_CA_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_CA_elev, aes(x = SD_lm_CA_elev$fitted.values, y = SD_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_CA_elev)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CS_elev, aes(x= SD_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CS_elev, aes(sample = SD_lm_CS_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_CS_elev, aes(x = SD_lm_CS_elev$fitted.values, y = SD_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_CS_elev)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_DBH_elev, aes(x= SD_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_DBH_elev, aes(sample = SD_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_DBH_elev, aes(x = SD_lm_DBH_elev$fitted.values, y = SD_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_DBH_elev)


### Sizes vs. Slope ###

# LM 

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x= LM_slope_raster_50_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")


#creating the linear regression
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_slope, aes(x= LM_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_slope, aes(sample = LM_lm_sca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_sca_slope, aes(x = LM_lm_sca_slope$fitted.values, y = LM_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_sca_slope)

#correlation test
cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts, LM_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_50_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression

LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_lca_slope, aes(x= LM_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_lca_slope, aes(sample = LM_lm_lca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_lca_slope, aes(x = LM_lm_lca_slope$fitted.values, y = LM_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_lca_slope)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_50_data_pts, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#linear regression with log transformation of canopy area
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#linear regression with square root transformation of canopy area
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CA_slope, aes(x= LM_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CA_slope, aes(sample = LM_lm_CA_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CA_slope, aes(x = LM_lm_CA_slope$fitted.values, y = LM_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_CA_slope)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_50_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression

LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CS_slope, aes(x= LM_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CS_slope, aes(sample = LM_lm_CS_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CS_slope, aes(x = LM_lm_CS_slope$fitted.values, y = LM_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_CS_slope)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_50_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#creating the linear regression
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#linear regression with logged transformation of aggregated DBH
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)

#linear regression with square root transformation of aggregated DBH
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_DBH_slope, aes(x= LM_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_DBH_slope, aes(sample = LM_lm_DBH_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_DBH_slope, aes(x = LM_lm_DBH_slope$fitted.values, y = LM_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_DBH_slope)



#LC linear models

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_50_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")


#creating the linear regression
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear regression with logged transformation of short canopy axis
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_sca_slope, aes(x= LC_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_slope, aes(sample = LC_lm_sca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_sca_slope, aes(x = LC_lm_sca_slope$fitted.values, y = LC_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LC_lm_sca_slope)

#correlation test
cor.test(LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts, LC_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_50_data_pts, y=Canopy_long_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear transformation with logged long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_lca_slope, aes(x= LC_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_lca_slope, aes(sample = LC_lm_lca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_lca_slope, aes(x = LC_lm_lca_slope$fitted.values, y = LC_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LC_lm_lca_slope)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_50_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear regression with log transformation of canopy area
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear regression with square root transformation of canopy area
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CA_slope, aes(x= LC_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CA_slope, aes(sample = LC_lm_CA_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_CA_slope, aes(x = LC_lm_CA_slope$fitted.values, y = LC_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LC_lm_CA_slope)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_50_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear transformation with logged crown spread
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear transformation with square rooted crown spread
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CS_slope, aes(x= LC_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CS_slope, aes(sample = LC_lm_CS_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_CS_slope, aes(x = LC_lm_CS_slope$fitted.values, y = LC_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LC_lm_CS_slope)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_50_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope")+
  ylab("DBH")

#creating the linear regression
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear regression with logged transformation of aggregated DBH
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)

#linear regression with square root transformation of aggregated DBH
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_DBH_slope, aes(x= LC_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_DBH_slope, aes(sample = LC_lm_DBH_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_lm_DBH_slope, aes(x = LC_lm_DBH_slope$fitted.values, y = LC_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LC_lm_DBH_slope)


#SD linear models

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_50_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#creating the linear regression
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with log transformation of canopy area
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_sca_slope, aes(x= SD_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_slope, aes(sample = SD_lm_sca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_sca_slope, aes(x = SD_lm_sca_slope$fitted.values, y = SD_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_sca_slope)

#correlation test
cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts, SD_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_50_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with log transformation of canopy area
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_lca_slope, aes(x= SD_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_lca_slope, aes(sample = SD_lm_lca_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_lca_slope, aes(x = SD_lm_lca_slope$fitted.values, y = SD_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_lca_slope)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_50_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with log transformation of canopy area
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CA_slope, aes(x= SD_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CA_slope, aes(sample = SD_lm_CA_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_CA_slope, aes(x = SD_lm_CA_slope$fitted.values, y = SD_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_CA_slope)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_50_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with logged transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with square root transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CS_slope, aes(x= SD_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CS_slope, aes(sample = SD_lm_CS_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_CS_slope, aes(x = SD_lm_CS_slope$fitted.values, y = SD_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_CS_slope)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_50_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#creating the linear regression
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with logged transformation of aggregated DBH
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)

#linear regression with square root transformation of aggregated DBH
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_50_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_DBH_slope, aes(x= SD_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_DBH_slope, aes(sample = SD_lm_DBH_slope$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_DBH_slope, aes(x = SD_lm_DBH_slope$fitted.values, y = SD_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_DBH_slope)

## Size vs. Aspect ##

#8 categories for direction

# LM

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_50_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LM_aov_SCA_aspect <- aov(Canopy_short ~ LM_aspect_raster_50_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_SCA_aspect <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_short, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_50_data_pts_8_categorical, p.adj = "bonf")

#long canopy axis
# canopy area
#crown spread
#DBH ag

# LC

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_50_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LC_aov_SCA_aspect <- aov(Canopy_short ~ LC_aspect_raster_50_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_50_data_pts_8_categorical, p.adj = "bonf")

#long canopy axis
# canopy area
#crown spread
#DBH ag

# SD

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_50_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
SD_aov_SCA_aspect <- aov(Canopy_short ~ SD_aspect_raster_50_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_50_data_pts_8_categorical, p.adj = "bonf")


#long canopy axis
# canopy area
#crown spread
#DBH ag

#4 categories for direction

# LM

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_50_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LM_aov_SCA_aspect <- aov(Canopy_short ~ LM_aspect_raster_50_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_SCA_aspect <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_short, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_50_data_pts_4_categorical, p.adj = "bonf")

#long canopy axis
# canopy area
#crown spread
#DBH ag

# LC

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_50_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LC_aov_SCA_aspect <- aov(Canopy_short ~ LC_aspect_raster_50_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_50_data_pts_4_categorical, p.adj = "bonf")

#long canopy axis
# canopy area
#crown spread
#DBH ag

# SD

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_50_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
SD_aov_SCA_aspect <- aov(Canopy_short ~ SD_aspect_raster_50_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_50_data_pts_8_categorical, p.adj = "bonf")


#long canopy axis
# canopy area
#crown spread
#DBH ag
