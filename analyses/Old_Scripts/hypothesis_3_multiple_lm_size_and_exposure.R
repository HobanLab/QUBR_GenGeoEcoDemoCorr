#because hypothesis 3 was becoming a long script of code, 
#I decided to separate out the multiple linear regression code into this script

#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(car) #to create added variable plots and to run levene's test for checking ANOVA conditions
library(stars) # to convert raster into stars
library(gdalUtilities) #to be able to use gdalwarp
library(npreg) #to use the gsm function for the generalized additive models
library(starsExtra) #to use dist_to_nearest
library(MuMIn) #to use the dredge function
library(rpart) #to use the function rpart to check recurissive binary


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
#fixed_field_data_processed_sf_trans_coordinates_download <- write.csv(fixed_field_data_processed_sf_trans_coordinates, "/Users/chewbecca/Morton Arboretum REU 2024/Untitled/QUBR_GenGeoEcoDemoCorr/data/fixed_field_data_processed_sf_trans_coordinates.csv", row.names = F)
#View(fixed_field_data_processed_sf_transformed)

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

#set elevation as a numeric value
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))


#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")

LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
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

#plotting the tree points by elevation (m)
ggplot()+
  geom_sf(data = LM_fixed_field_data_processed, aes(color = Elevation..m.FIXED))

ggplot()+
  geom_sf(data = fixed_field_data_processed_sf_trans_coordinates, aes(color = Elevation..m.FIXED))


##creating a new elevation column so that the 360 m outlier is 460
SD_fixed_field_data_processed <-  SD_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. == 360) ~ NA, 
                                        (Elevation..m. != 360) ~ Elevation..m.))

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

#creating the aspect and slope rasters

#elevation data from INEGI 15 m, so we can calculate slope and aspect


#BECAUSE THE ELEVATION RASTERS WERE TOO BIG TO DOWNLOAD DIRECTLY, FROM GOOGLE DRIVE, OR OPEN FROM A ZIP, 
#After LOADING IN THE ORIGINAL DATA, WE CROPPED IT TO FIT OUR POPULATIONS, EXPORTED THOSE FILES, AND THEN DOWNLOADED THOSE
#HERE IS A LINK TO A GOOGLE DRIVE WITH THE INGEI 15 m continuous elevation model DATA: https://drive.google.com/drive/folders/17RxjebifsRFFS4iEucDQtMFaqjzRI-Ss?usp=sharing 

#SO WE COMMENTED OUT CODE THAT IS HOW WE LOADED IN THE ORIGINAL RASTER AND CREATED AND EXPORTED THE CROPPED RASTERS FOR EACH POPULATION

# #projecting the INGEI 15 m continuous elevation model into UTM 12N 
# gdalwarp(srcfile = './data/15 m Elevation Raster/CEM bcs 15 m INEGI/CEM_V3_20170619_R15_E03_TIF/BajaCaliforniaS_15m.tif',  
#          dstfile = './data/15 m Elevation Raster/CEM_15_utm.tif', 
#          s_srs = '+proj=longlat +ellps=GRS80 +no_defs', 
#          t_srs= '+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 
#          tr = c(15, 15), overwrite=T)
# 
# #loading in the projected file
# CEM_15_utm <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm.tif"))
# 
# plot(CEM_15_utm) # just to visualize the raster
# 
# #cropping the rasters for each population
# 
# #all points
# 
# #mapping cropped 
# CEM_15_utm_all_points <- crop(CEM_15_utm, extent((c(LM_box[1]-100, SD_box[3]+100, SD_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the all points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_all_points, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = fixed_field_data_processed_sf_trans_coordinates)
# 
# #LM
# 
# #mapping cropped 
# CEM_15_utm_LM <- crop(CEM_15_utm, extent((c(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LM, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LM_fixed_field_data_processed)
# 
# #LC
# 
# #mapping cropped 
# CEM_15_utm_LC <- crop(CEM_15_utm, extent((c(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LC, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LC_fixed_field_data_processed)
# 
# #SD
# 
# #mapping cropped 
# CEM_15_utm_SD <- crop(CEM_15_utm, extent((c(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_SD, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = SD_fixed_field_data_processed)

# #exporting this cropped rasters as a tif
# writeRaster(CEM_15_utm_LM$CEM_15_utm,'./data/15 m Elevation Raster/CEM_15_utm_LM.tif') # sending the raster to the data folder and then to the 15 m elevation raster folder
# writeRaster(CEM_15_utm_LC$CEM_15_utm,'./data/15 m Elevation Raster/CEM_15_utm_LC.tif')
# writeRaster(CEM_15_utm_SD$CEM_15_utm,'./data/15 m Elevation Raster/CEM_15_utm_SD.tif')

# #projecting the INGEI 14 m continuous elevtion model into UTM 12N 
# gdalwarp(srcfile = "./data/CEM bcs 15 m INEGI/CEM_V3_20170619_R15_E03_TIF/BajaCaliforniaS_15m.tif", 
#          dstfile = "./data/CEM bcs 15 m INEGI/CEM_V3_20170619_R15_E03_TIF/CEM_15_utm.tif", 
#          s_srs = '+proj=longlat +ellps=GRS80 +no_defs', 
#          t_srs= '+proj=utm +zone=12 +ellps=GRS80 +towgs84=0,0,0,0,0,0,0 +units=m +no_defs', 
#          tr = c(15, 15), overwrite=TRUE)
# 
# #loading in the projected file
# CEM_15_utm <- raster(paste0("./data/CEM bcs 15 m INEGI/CEM_V3_20170619_R15_E03_TIF/CEM_15_utm.tif"))
# 
# plot(CEM_15_utm)

#HERE IS THE IMPORTATION OF THE CROPPED RASTERS

#cropping the rasters for each population

#all points
# 
# #mapping cropped 
# CEM_15_utm_all_points <- crop(CEM_15_utm, extent((c(LM_box[1]-100, SD_box[3]+100, SD_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the all points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_all_points, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = fixed_field_data_processed_sf_trans_coordinates)
# 
# #LM
# 
# #mapping cropped 
# CEM_15_utm_LM <- crop(CEM_15_utm, extent((c(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LM, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LM_fixed_field_data_processed)
# 
# #LC
# 
# #mapping cropped 
# CEM_15_utm_LC <- crop(CEM_15_utm, extent((c(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LC, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LC_fixed_field_data_processed)
# 
# #SD
# 
# #mapping cropped 
# CEM_15_utm_SD <- crop(CEM_15_utm, extent((c(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_SD, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = SD_fixed_field_data_processed)


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

#SD
SD_aspect_raster_15_data_pts <- extract(SD_aspect_raster_15, SD_fixed_field_data_processed) #extracting aspect for each point value
SD_slope_raster_15_data_pts <- extract(SD_slope_raster_15, SD_fixed_field_data_processed) #extracting slope for each point value
SD_elevation_raster_15_data_pts <- extract(CEM_15_utm_SD, SD_fixed_field_data_processed) #extracting the elevation for each point value
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed, SD_aspect_raster_15_data_pts) #bind the aspect data for each point to the SD point dataframe
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed_terrain, SD_slope_raster_15_data_pts) #bind the slope data for each point to the SD point dataframe
SD_fixed_field_data_processed_terrain <- cbind(SD_fixed_field_data_processed_terrain, SD_elevation_raster_15_data_pts) #bind the elevation data for each point to the LM point dataframe


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
                                                                        (all_points_aspect_raster_15_data_pts >= 337.5 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 67.5 & all_points_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 112.5 & all_points_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                        (all_points_aspect_raster_15_data_pts >= 157.5 & all_points_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                        (all_points_aspect_raster_15_data_pts >= 202.5 & all_points_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                        (all_points_aspect_raster_15_data_pts >= 247.5 & all_points_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 292.5 & all_points_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees


# North, East, South, West

# the directions are a range of 90 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_4_categorical = case_when((all_points_aspect_raster_15_data_pts >= 0 & all_points_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


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

#LM

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#LC

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed) + # Generate the base plot
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
  geom_histogram(aes(x = LM_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + # Generate the base plot
  geom_histogram(aes(x = LC_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + # Generate the base plot
  geom_histogram(aes(x = SD_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

#barcharts for aspect

# 8 categories of direction

ggplot(LM_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LM_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LC_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = SD_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

# 4 categories of direction

ggplot(LM_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LM_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = LC_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_bar(aes(x = SD_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")



#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)


#### Multiple Linear Regressions ####

#using only the 8 categories

# all points 

#Cook's D
plot(all_points_multiple_lm_SCA)
all_points_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_mlm_SCA_cooks <- cooks.distance(all_points_mlm_SCA) #calculating the cook.s D for each point
plot(all_points_mlm_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- all_points_mlm_SCA_cooks[(all_points_mlm_SCA_cooks > (2 * mean(all_points_mlm_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
all_points_fixed_field_data_processed_terrain_no_NA <- all_points_fixed_field_data_processed_terrain %>%
  filter(is.na(all_points_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(all_points_aspect_raster_15_data_pts_8_categorical) == F)

# SCA

#removing outliers based on which points were deemed influential
all_points_fixed_field_data_processed_terrain_no_NA_No_outliers <- all_points_fixed_field_data_processed_terrain_no_NA[-c(210, 303, 645),] #any with cooks above 0.015

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
all_points_multiple_lm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                 data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(all_points_multiple_lm_SCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing the summary of the model
all_points_multiple_lm_SCA_summary <- summary(all_points_multiple_lm_SCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
all_points_multiple_lm_SCA_vif <- car::vif(all_points_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
all_points_multiple_lm_SCA_VIF_multi_num <- (1 / (1-all_points_multiple_lm_SCA_summary$r.squared))
all_points_multiple_lm_SCA_vif > all_points_multiple_lm_SCA_VIF_multi_num

#checking VIF with reduced models
    #aspect
lm_check_no_aspect <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
lm_check_no_aspect_summary <- summary(lm_check_no_aspect)
lm_check_no_aspect_vif <- car::vif(lm_check_no_aspect) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
lm_check_no_aspect_VIF_multi_num <- (1 / (1-lm_check_no_aspect_summary$r.squared))
lm_check_no_aspect_vif > lm_check_no_aspect_VIF_multi_num

  #slope
lm_check_no_slope <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
lm_check_no_slope_summary <- summary(lm_check_no_slope)
lm_check_no_slope_vif <- car::vif(lm_check_no_slope) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
lm_check_no_slope_VIF_multi_num <- (1 / (1-lm_check_no_slope_summary$r.squared))
lm_check_no_slope_vif > lm_check_no_slope_VIF_multi_num

  #aspect
lm_check_no_elevation <-lm(Canopy_short ~  all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
lm_check_no_elevation_summary <- summary(lm_check_no_elevation)
lm_check_no_elevation_vif <- car::vif(lm_check_no_elevation) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
lm_check_no_elevation_VIF_multi_num <- (1 / (1-lm_check_no_elevation_summary$r.squared))
lm_check_no_elevation_vif > lm_check_no_elevation_VIF_multi_num


#determining our main effects model with two methods: backward's regression and the dredge function 
step(all_points_multiple_lm_SCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(all_points_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
all_points_multiple_lm_SCA_simplified <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_SCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(all_points_multiple_lm_SCA_simplified, all_points_multiple_lm_SCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#the simplified model is the same as the original one

#determing interactions with recursive binary partioning and regression tree
all_points_potential_interactions_SCA <- rpart(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                         all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(all_points_potential_interactions_SCA) # show the tree structure
text(all_points_potential_interactions_SCA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
all_points_multiple_lm_SCA_interacts <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                     all_points_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:all_points_slope_raster_15_data_pts +
                                     I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
summary(all_points_multiple_lm_SCA_interacts)

#slimming down the variables in the interaction model
step(all_points_multiple_lm_SCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(all_points_multiple_lm_SCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model

#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_SCA_interacts_simplified_step <- lm(Canopy_short ~  Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                                             all_points_aspect_raster_15_data_pts_8_categorical + I(all_points_slope_raster_15_data_pts^2) + 
                                                             I(Elevation..m.FIXED^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical,  data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
all_points_multiple_lm_SCA_interacts_simplified_dredge <- lm(Canopy_short ~  all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical +
                                                               I(all_points_slope_raster_15_data_pts^2) + Elevation..m.FIXED + I(Elevation..m.FIXED^2) +
                                                               all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical,
                                                             data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers)
#nested F test to compare simplified dredge and step models
anova(all_points_multiple_lm_SCA_interacts_simplified_step, all_points_multiple_lm_SCA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one


summary(all_points_multiple_lm_SCA_interacts_simplified_dredge) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(all_points_multiple_lm_SCA_interacts_simplified_dredge, all_points_multiple_lm_SCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(all_points_multiple_lm_SCA_interacts_simplified_dredge, all_points_multiple_lm_SCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: all_points_multiple_lm_SCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking linearity with elevation
ggplot(data = all_points_multiple_lm_SCA_simplified, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#checking linearity with slope
ggplot(data = all_points_multiple_lm_SCA_simplified, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#checking linearity by category
all_points_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED
  #elevation
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 

  #slope
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 
  
#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_SCA_simplified, aes(x= all_points_multiple_lm_SCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_SCA_simplified, aes(sample = all_points_multiple_lm_SCA_simplified$residuals))+
  geom_qq()

shapiro.test(all_points_multiple_lm_SCA_simplified$residuals) #shapiro wilk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
all_points_multiple_lm_SCA_simplified_lg <- lm(Canopy_short_lg ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_SCA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)

shapiro.test(all_points_multiple_lm_SCA_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the all_points_multiple_lm_SCA_simplified_lg  Shapiro-Wilk test we need to use a non-parametric model

#checking equal variance
ggplot(data = all_points_multiple_lm_SCA_simplified_lg, aes(x = all_points_multiple_lm_SCA_simplified_lg$fitted.values, y = all_points_multiple_lm_SCA_simplified_lg$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for log(SCA) and Elevation + Slope + Aspect")

#extracting model characteristics and significant
all_points_multiple_lm_SCA_summary <- summary(all_points_multiple_lm_SCA) #sig
all_points_multiple_lm_SCA_simplified_summary <- summary(all_points_multiple_lm_SCA_simplified) #sig
all_points_multiple_lm_SCA_simplified_lg_summary <- summary(all_points_multiple_lm_SCA_simplified_lg) #sig

#Because I could not get transformations and outliers to help the data meet the condition of normalized residuals, we will be using a generalized additive model (nonparametric)
#I am still using the best model: all_points_multiple_lm_SCA_simplified

#additive model
library(mgcv)
all_points_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
             data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA_interact <- gam(Canopy_short ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                             data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed)
gam.check(all_points_add.gam_SCA) 
gam.check(all_points_add.gam_SCA_interact)

#comparing the models' AIC
AIC(all_points_add.gam_SCA, all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed_first_term, 
    all_points_add.gam_SCA.smoothed_second_term, all_points_add.gam_SCA_interact)

#comparing the model's the models summary values
summary(all_points_add.gam_SCA)
summary(all_points_add.gam_SCA.smoothed)
summary(all_points_add.gam_SCA_interact)

#slimming down the variables in the best model
dredge <- dredge(all_points_add.gam_SCA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model

all_points_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#While the dredged model does not include Aspect, I will keep it in because we are curious in how it influences the size variable

#Chosen model: all_points_add.gam_SCA.smoothed

#updating K values
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_SCA.smoothed)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_SCA.smoothed, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")

library(gratia)
draw.gam
draw.gam(all_points_add.gam_SCA.smoothed)

ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA, aes(y =  ))

# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, smooth = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, smooth = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()

# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

library(visreg)
p3 <- visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
library(gridExtra)
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
library(plotly)
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
vis.gam(all_points_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_SCA <- fitted.values(all_points_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                 z=Canopy_short, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#have yet to get this version to work

library(mgcViz)
library(rgl)
all_points_add.gam_SCA <- getViz(all_points_add.gam_SCA)
plotRGL(sm(all_points_add.gam_SCA, 1), fix = c("Canopy_short" = 1), residuals = TRUE)


ggplot() +
  geom_point(data = all_points_fixed_field_data_processed_terrain_no_NA, aes(x = hp, y = mpg)) +prin
  geom_line(data = data.frame(hp = new_data$hp, mpg = predictions$fit), 
            aes(x = hp, y = mpg), color = "blue", size = 1) +
  geom_ribbon(data = data.frame(hp = new_data$hp, fit = predictions$fit, 
                                se = predictions$se.fit), aes(x = hp, 
                                                              ymin = fit - 1.96 * se, 
                                                              ymax = fit + 1.96 * se), alpha = 0.3) +
  labs(title = "Generalized Additive Model (GAM) Fit for mpg vs. hp", 
       x = "Horsepower", y = "Miles per Gallon") +
  theme_minimal()



#The first model has a lower 

# generalized smooth model (using all knots)
all_points_add.gsm_SCA <- gsm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
               data = all_points_fixed_field_data_processed_terrain_no_NA, knots = nrow(all_points_fixed_field_data_processed_terrain_no_NA)) 



summary(all_points_add.gsm_SCA) #summarize the model

#interaction model, difference is the * instead of +
all_points_add.gsm_SCA_interact <- gsm(Canopy_short ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, knots = nrow(all_points_fixed_field_data_processed_terrain_no_NA)) 
all_points_add.gsm_SCA_interact
summary(all_points_add.gsm_SCA_interact) #summarize the model

#comparing the additive and interaction models

# pseudo F test of interaction effect
sse.dif <- all_points_add.gsm_SCA$deviance - all_points_add.gsm_SCA_interact$deviance
df.dif <- all_points_add.gsm_SCA_interact$df - all_points_add.gsm_SCA$df
Fstat <- (sse.dif / df.dif) / all_points_add.gsm_SCA_interact$dispersion
pvalue <- 1 - pf(Fstat, df1 = df.dif, df2 = nrow(Prestige) - all_points_add.gsm_SCA_interact$df)
Fstat
pvalue #not significant, so do not need to consider the interactions

#GSM with gamma and inverse guassian

all_points_add.gsm_SCA_gamma <- gsm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, knots = nrow(all_points_fixed_field_data_processed_terrain_no_NA), family = "Gamma") 

all_points_add.gsm_SCA_inv_guass <- gsm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                    data = all_points_fixed_field_data_processed_terrain_no_NA, knots = nrow(all_points_fixed_field_data_processed_terrain_no_NA), family = "inverse.guassian") 

summary(all_points_add.gsm_SCA_gamma)


#GLM model to find regression for non-normal data

library(DHARMa)
library(aod)

#gamma
all_points_glm_gamma_SCA <- glm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                data = all_points_fixed_field_data_processed_terrain_no_NA, family = "Gamma")


#inverse guassian, for when dependent variable is even more skewed to the right
all_points_glm_inv_guass_SCA <- glm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = "inverse.gaussian")
  
#comparing the results
summary(all_points_multiple_lm_SCA_simplified)
all_points_glm_gamma_SCA_summary <- summary(all_points_glm_gamma_SCA) #results of gamma family, lower aic
all_points_glm_inv_guass_SCA_SCA_summary <- summary(all_points_glm_inv_guass_SCA) #results of inverse guassian

    #the gamma model has a lower AIC value, stronger model

#testing model fit with ANOVA F test and ANOVA Likelihood Ratio Test
anova(all_points_glm_gamma_SCA, test = "F")
anova(all_points_glm_gamma_SCA, test = "LRT")
anova(all_points_glm_gamma_SCA, test = "F")
anova(all_points_glm_gamma_SCA, test = "Chisq")

#dredging to see if there is a better simplified model
dredge <- dredge(all_points_glm_gamma_SCA) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model

#best result: all_points_glm_gamma_SCA

#deviance test, high value means good fit
1 - pchisq(deviance(all_points_glm_gamma_SCA_summary), df.residual(all_points_glm_gamma_SCA_summary))
#good when it is approximately 1, small suggests not fitting weel

# Attempt at using the Pearson chi-square statistic to determine the goodness of fit of the model, high value is good fit

# Compute Pearson chi-square statistic
pearson_chi_sq <- sum(residuals(all_points_glm_gamma_SCA, type = "pearson")^2)

# Get degrees of freedom
df <- df.residual(all_points_glm_gamma_SCA)

# Compute p-value
p_value <- pchisq(pearson_chi_sq, df, lower.tail = F)

# Print results
cat("Pearson Chi-Square Statistic:", pearson_chi_sq, "\n")
cat("Degrees of Freedom:", df, "\n")
cat("P-value:", p_value, "\n") #failed goodness of fit test


# LCA

plot(all_points_multiple_lm_LCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
all_points_multiple_lm_LCA <- lm(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(all_points_multiple_lm_LCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

all_points_multiple_lm_LCA_summary <- summary(all_points_multiple_lm_LCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
all_points_multiple_lm_LCA_vif <- car::vif(all_points_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
all_points_multiple_lm_LCA_VIF_multi_num <- (1 / (1-all_points_multiple_lm_LCA_summary$r.squared))
all_points_multiple_lm_LCA_vif > all_points_multiple_lm_LCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(all_points_multiple_lm_LCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(all_points_multiple_lm_LCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
all_points_multiple_lm_LCA_simplified <- lm(Canopy_long ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_LCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(all_points_multiple_lm_LCA_simplified, all_points_multiple_lm_LCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
all_points_potential_interactions_LCA <- rpart(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                         all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(all_points_potential_interactions_LCA) # show the tree structure
text(all_points_potential_interactions_LCA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree, branches mean that the variables likely have interactions with one another
all_points_multiple_lm_LCA_interacts <- lm(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                     all_points_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:all_points_slope_raster_15_data_pts +
                                     I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical, 
                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_LCA_interacts)

#slimming down the variables in the interaction model
step(all_points_multiple_lm_LCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(all_points_multiple_lm_LCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_fixed_field_data_processed_terrain_no_NA$I(all_points_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_LCA_interacts_simplified_step <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2) + I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts,  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_LCA_interacts_simplified_dredge <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_LCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(all_points_multiple_lm_LCA_interacts_simplified_step, all_points_multiple_lm_LCA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified


#nested F test to compare simplified interactions model made with dredge to full interactions model
anova(all_points_multiple_lm_LCA_interacts_simplified_dredge, all_points_multiple_lm_LCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model made with dredge to simplified model without interactions
anova(all_points_multiple_lm_LCA_interacts_simplified_dredge, all_points_multiple_lm_LCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: all_points_multiple_lm_LCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)
#checking linearity with elevation
ggplot(data = all_points_multiple_lm_LCA_simplified, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#checking linearity by category
all_points_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED
#elevation
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 

#slope
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_LCA_simplified, aes(x= all_points_multiple_lm_LCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_LCA_simplified, aes(sample = all_points_multiple_lm_LCA_simplified$residuals))+
  geom_qq()

shapiro.test(all_points_multiple_lm_LCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro test was singificant, so I will use a transformaed canopy_lung variable
all_points_multiple_lm_LCA_simplified_lg <- lm(Canopy_long_lg ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_LCA_simplified_sqrt <- lm(Canopy_long_sqrt ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)

shapiro.test(all_points_multiple_lm_LCA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking equal variance
ggplot(data = all_points_multiple_lm_LCA_simplified_sqrt, aes(x = all_points_multiple_lm_LCA_simplified_sqrt$fitted.values, y = all_points_multiple_lm_LCA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
all_points_multiple_lm_LCA_summary <- summary(all_points_multiple_lm_LCA)
all_points_multiple_lm_LCA_simplified_summary <- summary(all_points_multiple_lm_LCA_simplified_sqrt) #not significant, the square transformation changes the p-value a lot 
all_points_multiple_lm_LCA_simplified_summary <- summary(all_points_multiple_lm_LCA_simplified) #significant 

# CA

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
all_points_multiple_lm_CA <- lm(Canopy_area ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(all_points_multiple_lm_CA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

all_points_multiple_lm_CA_summary <- summary(all_points_multiple_lm_CA) #storing a summary of the model

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
all_points_multiple_lm_CA_vif <- car::vif(all_points_multiple_lm_CA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
all_points_multiple_lm_CA_VIF_multi_num <- (1 / (1-all_points_multiple_lm_CA_summary$r.squared))
all_points_multiple_lm_CA_vif > all_points_multiple_lm_CA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(all_points_multiple_lm_CA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(all_points_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
all_points_multiple_lm_CA_simplified <- lm(Canopy_area ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(all_points_multiple_lm_CA_simplified, all_points_multiple_lm_CA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
all_points_potential_interactions_CA <- rpart(Canopy_area ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                        all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(all_points_potential_interactions_CA) # show the tree structure
text(all_points_potential_interactions_CA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
all_points_multiple_lm_CA_interacts <- lm(Canopy_area ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                    all_points_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:all_points_slope_raster_15_data_pts +
                                    I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CA_interacts)

#slimming down the variables in the interaction model
step(all_points_multiple_lm_CA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(all_points_multiple_lm_CA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_fixed_field_data_processed_terrain_no_NA$I(all_points_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_CA_interacts_simplified_step <- lm(Canopy_area ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:all_points_slope_raster_15_data_pts,  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = all_points_fixed_field_data_processed_terrain_no_NA)

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(all_points_multiple_lm_CA_interacts_simplified_step, all_points_multiple_lm_CA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified

#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CA_interacts_simplified_dredge) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(all_points_multiple_lm_CA_interacts_simplified_dredge, all_points_multiple_lm_CA_interacts) #results are not significant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(all_points_multiple_lm_CA_interacts_simplified_dredge, all_points_multiple_lm_CA_simplified) #results are significant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: all_points_multiple_lm_CA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)
#checking linearity with elevation
ggplot(data = all_points_multiple_lm_CA_simplified, (aes(x=Elevation..m.FIXED, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Canopy Area")

#checking linearity by category
all_points_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED
#elevation
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=Elevation..m.FIXED, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Canopy Area")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 

#slope
ggplot(data = all_points_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")+
  facet_wrap(~ all_points_aspect_raster_15_data_pts_8_categorical) 


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_CA_simplified, aes(x= all_points_multiple_lm_CA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_CA_simplified, aes(sample = all_points_multiple_lm_CA_simplified$residuals))+
  geom_qq()

shapiro.test(all_points_multiple_lm_CA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
all_points_multiple_lm_CA_simplified_lg <- lm(Canopy_area_lg ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_CA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)

shapiro.test(all_points_multiple_lm_CA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation

#the results of the Shapiro-Wilk test suggest we should use the model where canopy area is square rooted

#checking equal variance
ggplot(data = all_points_multiple_lm_CA_simplified_sqrt, aes(x = all_points_multiple_lm_CA_simplified_sqrt$fitted.values, y = all_points_multiple_lm_CA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
all_points_multiple_lm_CA_summary <- summary(all_points_multiple_lm_CA) #sign
all_points_multiple_lm_CA_simplified_summary <- summary(all_points_multiple_lm_CA_simplified) #sign
all_points_multiple_lm_CA_simplified_sqrt_summary <- summary(all_points_multiple_lm_CA_simplified_sqrt) #sig


# CS
#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
all_points_multiple_lm_CS <- lm(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(all_points_multiple_lm_CS) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing a summary of the model
all_points_multiple_lm_CS_summary <- summary(all_points_multiple_lm_CS)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
all_points_multiple_lm_CS_vif <- car::vif(all_points_multiple_lm_CS) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
all_points_multiple_lm_CS_VIF_multi_num <- (1 / (1-all_points_multiple_lm_CS_summary$r.squared))
all_points_multiple_lm_CS_vif > all_points_multiple_lm_CS_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(all_points_multiple_lm_CS) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(all_points_multiple_lm_CS) #generates the best model and the rank of best models

#both the step and dredge technique produced the same simplified model:
all_points_multiple_lm_CS_simplified <- lm(Crown_spread ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CS_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(all_points_multiple_lm_CS_simplified, all_points_multiple_lm_CS) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#best simplified model without taking into account interactions: all_points_multiple_lm_CS_simplified

#determing interactions with recursive binary partioning and regression tree
all_points_potential_interactions_CS <- rpart(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                        all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(all_points_potential_interactions_CS) # show the tree structure
text(all_points_potential_interactions_CS, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
all_points_multiple_lm_CS_interacts <- lm(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                    all_points_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:all_points_slope_raster_15_data_pts +
                                    I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CS_interacts)

#slimming down the variables in the interaction model
step(all_points_multiple_lm_CS_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(all_points_multiple_lm_CS_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_fixed_field_data_processed_terrain_no_NA$I(all_points_aspect_raster_15_data_pts_8_categorical^2)

#the step and dredge methods produced the sample best model:
#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_CS_interacts_simplified <- lm(Crown_spread ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_CS_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(all_points_multiple_lm_CS_interacts_simplified, all_points_multiple_lm_CS_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(all_points_multiple_lm_CS_interacts_simplified, all_points_multiple_lm_CS_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: all_points_multiple_lm_CS_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_CS_simplified, aes(x= all_points_multiple_lm_CS_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_CS_simplified, aes(sample = all_points_multiple_lm_CS_simplified$residuals))+
  geom_qq()

shapiro.test(all_points_multiple_lm_CS_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
all_points_multiple_lm_CS_simplified_lg <- lm(Crown_spread_lg ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_CS_simplified_sqrt <- lm(Crown_spread_sqrt ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)

shapiro.test(all_points_multiple_lm_CS_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_CS_simplified_sqrt, aes(x= all_points_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for sqrt(Crown Spread) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_CS_simplified_sqrt, aes(sample = all_points_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_multiple_lm_CS_simplified_sqrt, aes(x = all_points_multiple_lm_CS_simplified_sqrt$fitted.values, y = all_points_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
all_points_multiple_lm_SCA_summary <- summary(all_points_multiple_lm_SCA) #sig
all_points_multiple_lm_SCA_simplified_summary <- summary(all_points_multiple_lm_SCA_simplified) #sig
all_points_multiple_lm_SCA_simplified_sqrt_summary <- summary(all_points_multiple_lm_CS_simplified_sqrt) #sig


# DBH_ag

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
all_points_multiple_lm_DBH <- lm(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(all_points_multiple_lm_DBH) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing summary of the model
all_points_multiple_lm_DBH_summary <- summary(all_points_multiple_lm_DBH)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
all_points_multiple_lm_DBH_vif <- car::vif(all_points_multiple_lm_DBH) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
all_points_multiple_lm_DBH_VIF_multi_num <- (1 / (1-all_points_multiple_lm_DBH_summary$r.squared))
all_points_multiple_lm_DBH_vif > all_points_multiple_lm_DBH_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(all_points_multiple_lm_DBH) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(all_points_multiple_lm_DBH) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
all_points_multiple_lm_DBH_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_DBH_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(all_points_multiple_lm_DBH_simplified, all_points_multiple_lm_DBH) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
all_points_potential_interactions_DBH <- rpart(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                         all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(all_points_potential_interactions_DBH) # show the tree structure
text(all_points_potential_interactions_DBH, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
all_points_multiple_lm_DBH_interacts <- lm(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + 
                                     all_points_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:all_points_slope_raster_15_data_pts +
                                     I(all_points_slope_raster_15_data_pts^2) + all_points_slope_raster_15_data_pts:all_points_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:all_points_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_DBH_interacts)

#slimming down the variables in the interaction model
step(all_points_multiple_lm_DBH_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(all_points_multiple_lm_DBH_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_fixed_field_data_processed_terrain_no_NA$I(all_points_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
all_points_multiple_lm_DBH_interacts_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
summary(all_points_multiple_lm_DBH_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(all_points_multiple_lm_DBH_interacts_simplified, all_points_multiple_lm_DBH_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(all_points_multiple_lm_DBH_interacts_simplified, all_points_multiple_lm_DBH_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that the simplified interaction and simplified regular model are the same, so we could use either
# for the now we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: all_points_multiple_lm_DBH_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_DBH_simplified, aes(x= all_points_multiple_lm_DBH_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_DBH_simplified, aes(sample = all_points_multiple_lm_DBH_simplified$residuals))+
  geom_qq()

shapiro.test(all_points_multiple_lm_DBH_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
all_points_multiple_lm_DBH_simplified_lg <- lm(DBH_ag_lg ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_multiple_lm_DBH_simplified_sqrt <- lm(DBH_ag_sqrt ~ Elevation..m.FIXED, data = all_points_fixed_field_data_processed_terrain_no_NA)

shapiro.test(all_points_multiple_lm_DBH_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long varaible that has a logged transformation

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_multiple_lm_DBH_simplified_lg, aes(x= all_points_multiple_lm_DBH_simplified_lg$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for log(DBH) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_multiple_lm_DBH_simplified_lg, aes(sample = all_points_multiple_lm_DBH_simplified_lg$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_multiple_lm_DBH_simplified_lg, aes(x = all_points_multiple_lm_DBH_simplified_lg$fitted.values, y = all_points_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for log(DBH) and Elevation")

#extracting model characteristics and significant
all_points_multiple_lm_DBH_summary <- summary(all_points_multiple_lm_DBH) #not sig
all_points_multiple_lm_DBH_simplified_summary <- summary(all_points_multiple_lm_DBH_simplified) #sig
all_points_multiple_lm_DBH_simplified_lg_summary <- summary(all_points_multiple_lm_DBH_simplified_lg) #sig



# LM

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
LM_fixed_field_data_processed_terrain_no_NA <- LM_fixed_field_data_processed_terrain %>%
  filter(is.na(LM_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)

#Cook's D
LM_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_mlm_SCA_cooks <- cooks.distance(LM_mlm_SCA) #calculating the cook.s D for each point
plot(LM_mlm_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_mlm_SCA_cooks[(LM_mlm_SCA_cooks > (3 * mean(LM_mlm_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_terrain_no_NA_No_outliers <- LM_fixed_field_data_processed_terrain_no_NA[-c(24,26,27),]


# SCA

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_multiple_lm_SCA_summary <- summary(LM_multiple_lm_SCA)

plot(LM_multiple_lm_SCA)

#checking to see which variables might be the most useful
avPlots(LM_multiple_lm_SCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LM_multiple_lm_SCA_vif <- car::vif(LM_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LM_multiple_lm_SCA_VIF_multi_num <- (1 / (1-LM_multiple_lm_SCA_summary$r.squared))
LM_multiple_lm_SCA_vif > LM_multiple_lm_SCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LM_multiple_lm_SCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LM_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LM_multiple_lm_SCA_simplified <- lm(Canopy_short ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
summary(LM_multiple_lm_SCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LM_multiple_lm_SCA_simplified, LM_multiple_lm_SCA) #results are not significant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LM_potential_interactions <- rpart(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                     LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LM_potential_interactions) # show the tree structure
text(LM_potential_interactions, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LM_multiple_lm_SCA_interacts <- lm(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                     LM_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LM_slope_raster_15_data_pts +
                                     I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts:LM_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
summary(LM_multiple_lm_SCA_interacts)

#slimming down the variables in the interaction model
step(LM_multiple_lm_SCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LM_multiple_lm_SCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model


#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_SCA_interacts_simplified <- lm(Canopy_short ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
summary(LM_multiple_lm_SCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LM_multiple_lm_SCA_interacts_simplified, LM_multiple_lm_SCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LM_multiple_lm_SCA_interacts_simplified, LM_multiple_lm_SCA_simplified) #results are signficant, meaning there is compelling evidence to support the larger model than the smaller one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LM_multiple_lm_SCA_interacts_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking linearity with elevation
ggplot(data = LM_multiple_lm_SCA_interacts_simplified, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#checking linearity by category
LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED
#elevation
ggplot(data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")+
  facet_wrap(~ LM_aspect_raster_15_data_pts_8_categorical) 

#slope
ggplot(data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, (aes(x=LM_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")+
  facet_wrap(~ LM_aspect_raster_15_data_pts_8_categorical) 

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_SCA_interacts_simplified, aes(x= LM_multiple_lm_SCA_interacts_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_SCA_interacts_simplified, aes(sample = LM_multiple_lm_SCA_interacts_simplified$residuals))+
  geom_qq()

shapiro.test(LM_multiple_lm_SCA_interacts_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LM_multiple_lm_SCA_simplified_lg <- lm(Canopy_short_lg ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_multiple_lm_SCA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

shapiro.test(LM_multiple_lm_SCA_interacts_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the all_points_multiple_lm_SCA_simplified_lg  Shapiro-Wilk test we need to use a non-parametric test to look at slope, but we could use log to look at the 


#checking equal variance
ggplot(data = LM_multiple_lm_SCA, aes(x = LM_multiple_lm_SCA$fitted.values, y = LM_multiple_lm_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
LM_multiple_lm_SCA_summary <- summary(LM_multiple_lm_SCA)
LM_multiple_lm_SCA_simplified_summary <- summary(LM_multiple_lm_SCA_simplified)

View(LM_fixed_field_data_processed_terrain_no_NA)


# LCA

plot(LM_multiple_lm_LCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_LCA <- lm(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LM_multiple_lm_LCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

LM_multiple_lm_LCA_summary <- summary(LM_multiple_lm_LCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LM_multiple_lm_LCA_vif <- car::vif(LM_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LM_multiple_lm_LCA_VIF_multi_num <- (1 / (1-LM_multiple_lm_LCA_summary$r.squared))
LM_multiple_lm_LCA_vif > LM_multiple_lm_LCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LM_multiple_lm_LCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LM_multiple_lm_LCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LM_multiple_lm_LCA_simplified <- lm(Canopy_long ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_LCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LM_multiple_lm_LCA_simplified, LM_multiple_lm_LCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LM_potential_interactions_LCA <- rpart(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                         LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LM_potential_interactions_LCA) # show the tree structure
text(LM_potential_interactions_LCA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree, branches mean that the variables likely have interactions with one another
LM_multiple_lm_LCA_interacts <- lm(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                     LM_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LM_slope_raster_15_data_pts +
                                     I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts:LM_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_LCA_interacts)

#slimming down the variables in the interaction model
step(LM_multiple_lm_LCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LM_multiple_lm_LCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_fixed_field_data_processed_terrain_no_NA$I(LM_aspect_raster_15_data_pts_8_categorical^2)


#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_LCA_interacts_simplified_step <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2) + I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts,  data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_LCA_interacts_simplified_dredge <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_LCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(LM_multiple_lm_LCA_interacts_simplified_step, LM_multiple_lm_LCA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified


#nested F test to compare simplified interactions model made with dredge to full interactions model
anova(LM_multiple_lm_LCA_interacts_simplified_dredge, LM_multiple_lm_LCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model made with dredge to simplified model without interactions
anova(LM_multiple_lm_LCA_interacts_simplified_dredge, LM_multiple_lm_LCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LM_multiple_lm_LCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_LCA_simplified, aes(x= LM_multiple_lm_LCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_LCA_simplified, aes(sample = LM_multiple_lm_LCA_simplified$residuals))+
  geom_qq()

shapiro.test(LM_multiple_lm_LCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro test was singificant, so I will use a transformaed canopy_lung variable
LM_multiple_lm_LCA_simplified_lg <- lm(Canopy_long_lg ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_LCA_simplified_sqrt <- lm(Canopy_long_sqrt ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LM_multiple_lm_LCA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking equal variance
ggplot(data = LM_multiple_lm_LCA_simplified_sqrt, aes(x = LM_multiple_lm_LCA_simplified_sqrt$fitted.values, y = LM_multiple_lm_LCA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
LM_multiple_lm_LCA_summary <- summary(LM_multiple_lm_LCA)
LM_multiple_lm_LCA_simplified_summary <- summary(LM_multiple_lm_LCA_simplified_sqrt) #not significant, the square transformation changes the p-value a lot 
LM_multiple_lm_LCA_simplified_summary <- summary(LM_multiple_lm_LCA_simplified) #significant 

# CA

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_CA <- lm(Canopy_area ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LM_multiple_lm_CA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

LM_multiple_lm_CA_summary <- summary(LM_multiple_lm_CA) #storing a summary of the model

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LM_multiple_lm_CA_vif <- car::vif(LM_multiple_lm_CA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LM_multiple_lm_CA_VIF_multi_num <- (1 / (1-LM_multiple_lm_CA_summary$r.squared))
LM_multiple_lm_CA_vif > LM_multiple_lm_CA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LM_multiple_lm_CA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LM_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LM_multiple_lm_CA_simplified <- lm(Canopy_area ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LM_multiple_lm_CA_simplified, LM_multiple_lm_CA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LM_potential_interactions_CA <- rpart(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                        LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LM_potential_interactions_CA) # show the tree structure
text(LM_potential_interactions_CA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LM_multiple_lm_CA_interacts <- lm(Canopy_area ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                    LM_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LM_slope_raster_15_data_pts +
                                    I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts:LM_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CA_interacts)

#slimming down the variables in the interaction model
step(LM_multiple_lm_CA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LM_multiple_lm_CA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_fixed_field_data_processed_terrain_no_NA$I(LM_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_CA_interacts_simplified_step <- lm(Canopy_area ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:LM_slope_raster_15_data_pts,  data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LM_fixed_field_data_processed_terrain_no_NA)

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(LM_multiple_lm_CA_interacts_simplified_step, LM_multiple_lm_CA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified

#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CA_interacts_simplified_dredge) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LM_multiple_lm_CA_interacts_simplified_dredge, LM_multiple_lm_CA_interacts) #results are not significant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LM_multiple_lm_CA_interacts_simplified_dredge, LM_multiple_lm_CA_simplified) #results are significant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LM_multiple_lm_CA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_CA_simplified, aes(x= LM_multiple_lm_CA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_CA_simplified, aes(sample = LM_multiple_lm_CA_simplified$residuals))+
  geom_qq()

shapiro.test(LM_multiple_lm_CA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LM_multiple_lm_CA_simplified_lg <- lm(Canopy_area_lg ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_CA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LM_multiple_lm_CA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation

#the results of the Shapiro-Wilk test suggest we should use the model where canopy area is square rooted

#checking equal variance
ggplot(data = LM_multiple_lm_CA_simplified_sqrt, aes(x = LM_multiple_lm_CA_simplified_sqrt$fitted.values, y = LM_multiple_lm_CA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
LM_multiple_lm_CA_summary <- summary(LM_multiple_lm_CA) #sign
LM_multiple_lm_CA_simplified_summary <- summary(LM_multiple_lm_CA_simplified) #sign
LM_multiple_lm_CA_simplified_sqrt_summary <- summary(LM_multiple_lm_CA_simplified_sqrt) #sig


# CS
#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_CS <- lm(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LM_multiple_lm_CS) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing a summary of the model
LM_multiple_lm_CS_summary <- summary(LM_multiple_lm_CS)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LM_multiple_lm_CS_vif <- car::vif(LM_multiple_lm_CS) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LM_multiple_lm_CS_VIF_multi_num <- (1 / (1-LM_multiple_lm_CS_summary$r.squared))
LM_multiple_lm_CS_vif > LM_multiple_lm_CS_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LM_multiple_lm_CS) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LM_multiple_lm_CS) #generates the best model and the rank of best models

#both the step and dredge technique produced the same simplified model:
LM_multiple_lm_CS_simplified <- lm(Crown_spread ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CS_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LM_multiple_lm_CS_simplified, LM_multiple_lm_CS) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#best simplified model without taking into account interactions: LM_multiple_lm_CS_simplified

#determing interactions with recursive binary partioning and regression tree
LM_potential_interactions_CS <- rpart(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                        LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LM_potential_interactions_CS) # show the tree structure
text(LM_potential_interactions_CS, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LM_multiple_lm_CS_interacts <- lm(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                    LM_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LM_slope_raster_15_data_pts +
                                    I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts:LM_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CS_interacts)

#slimming down the variables in the interaction model
step(LM_multiple_lm_CS_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LM_multiple_lm_CS_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_fixed_field_data_processed_terrain_no_NA$I(LM_aspect_raster_15_data_pts_8_categorical^2)

#the step and dredge methods produced the sample best model:
#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_CS_interacts_simplified <- lm(Crown_spread ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_CS_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LM_multiple_lm_CS_interacts_simplified, LM_multiple_lm_CS_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LM_multiple_lm_CS_interacts_simplified, LM_multiple_lm_CS_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LM_multiple_lm_CS_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_CS_simplified, aes(x= LM_multiple_lm_CS_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_CS_simplified, aes(sample = LM_multiple_lm_CS_simplified$residuals))+
  geom_qq()

shapiro.test(LM_multiple_lm_CS_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LM_multiple_lm_CS_simplified_lg <- lm(Crown_spread_lg ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_CS_simplified_sqrt <- lm(Crown_spread_sqrt ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LM_multiple_lm_CS_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_CS_simplified_sqrt, aes(x= LM_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for sqrt(Crown Spread) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_CS_simplified_sqrt, aes(sample = LM_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_multiple_lm_CS_simplified_sqrt, aes(x = LM_multiple_lm_CS_simplified_sqrt$fitted.values, y = LM_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
LM_multiple_lm_SCA_summary <- summary(LM_multiple_lm_SCA) #sig
LM_multiple_lm_SCA_simplified_summary <- summary(LM_multiple_lm_SCA_simplified) #sig
LM_multiple_lm_SCA_simplified_sqrt_summary <- summary(LM_multiple_lm_CS_simplified_sqrt) #sig


# DBH_ag

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_DBH <- lm(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LM_multiple_lm_DBH) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing summary of the model
LM_multiple_lm_DBH_summary <- summary(LM_multiple_lm_DBH)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LM_multiple_lm_DBH_vif <- car::vif(LM_multiple_lm_DBH) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LM_multiple_lm_DBH_VIF_multi_num <- (1 / (1-LM_multiple_lm_DBH_summary$r.squared))
LM_multiple_lm_DBH_vif > LM_multiple_lm_DBH_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LM_multiple_lm_DBH) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LM_multiple_lm_DBH) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LM_multiple_lm_DBH_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_DBH_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LM_multiple_lm_DBH_simplified, LM_multiple_lm_DBH) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LM_potential_interactions_DBH <- rpart(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                         LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LM_potential_interactions_DBH) # show the tree structure
text(LM_potential_interactions_DBH, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LM_multiple_lm_DBH_interacts <- lm(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + 
                                     LM_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LM_slope_raster_15_data_pts +
                                     I(LM_slope_raster_15_data_pts^2) + LM_slope_raster_15_data_pts:LM_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LM_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_DBH_interacts)

#slimming down the variables in the interaction model
step(LM_multiple_lm_DBH_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LM_multiple_lm_DBH_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_fixed_field_data_processed_terrain_no_NA$I(LM_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
LM_multiple_lm_DBH_interacts_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
summary(LM_multiple_lm_DBH_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LM_multiple_lm_DBH_interacts_simplified, LM_multiple_lm_DBH_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LM_multiple_lm_DBH_interacts_simplified, LM_multiple_lm_DBH_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that the simplified interaction and simplified regular model are the same, so we could use either
# for the now we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LM_multiple_lm_DBH_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_DBH_simplified, aes(x= LM_multiple_lm_DBH_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_DBH_simplified, aes(sample = LM_multiple_lm_DBH_simplified$residuals))+
  geom_qq()

shapiro.test(LM_multiple_lm_DBH_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LM_multiple_lm_DBH_simplified_lg <- lm(DBH_ag_lg ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_multiple_lm_DBH_simplified_sqrt <- lm(DBH_ag_sqrt ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LM_multiple_lm_DBH_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long varaible that has a logged transformation

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_multiple_lm_DBH_simplified_lg, aes(x= LM_multiple_lm_DBH_simplified_lg$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for log(DBH) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_multiple_lm_DBH_simplified_lg, aes(sample = LM_multiple_lm_DBH_simplified_lg$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_multiple_lm_DBH_simplified_lg, aes(x = LM_multiple_lm_DBH_simplified_lg$fitted.values, y = LM_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for log(DBH) and Elevation")

#extracting model characteristics and significant
LM_multiple_lm_DBH_summary <- summary(LM_multiple_lm_DBH) #not sig
LM_multiple_lm_DBH_simplified_summary <- summary(LM_multiple_lm_DBH_simplified) #sig
LM_multiple_lm_DBH_simplified_lg_summary <- summary(LM_multiple_lm_DBH_simplified_lg) #sig



# LC

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
LC_fixed_field_data_processed_terrain_no_NA <- LC_fixed_field_data_processed_terrain %>%
  filter(is.na(LC_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)


# SCA

plot(LC_multiple_lm_SCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LC_multiple_lm_SCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing the summary of the model
LM_multiple_lm_SCA_summary <- summary(LC_multiple_lm_SCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LC_multiple_lm_SCA_vif <- car::vif(LC_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LC_multiple_lm_SCA_VIF_multi_num <- (1 / (1-LC_multiple_lm_SCA_summary$r.squared))
LC_multiple_lm_SCA_vif > LC_multiple_lm_SCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LC_multiple_lm_SCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LC_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LC_multiple_lm_SCA_simplified <- lm(Canopy_short ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_SCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LC_multiple_lm_SCA_simplified, LC_multiple_lm_SCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LC_potential_interactions_SCA <- rpart(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                     LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LC_potential_interactions) # show the tree structure
text(LC_potential_interactions, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LC_multiple_lm_SCA_interacts <- lm(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                     LC_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LC_slope_raster_15_data_pts +
                                     I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts:LC_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_SCA_interacts)

#slimming down the variables in the interaction model
step(LC_multiple_lm_SCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LC_multiple_lm_SCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_fixed_field_data_processed_terrain_no_NA$I(LC_aspect_raster_15_data_pts_8_categorical^2)


#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_SCA_interacts_simplified <- lm(Canopy_short ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_SCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LC_multiple_lm_SCA_interacts_simplified, LC_multiple_lm_SCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LC_multiple_lm_SCA_interacts_simplified, LC_multiple_lm_SCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LC_multiple_lm_SCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_SCA_simplified, aes(x= LC_multiple_lm_SCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_SCA_simplified, aes(sample = LC_multiple_lm_SCA_simplified$residuals))+
  geom_qq()

shapiro.test(LC_multiple_lm_SCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro test was singificant, so I will use a transformaed canopy_lung variable
LC_multiple_lm_SCA_simplified_lg <- lm(Canopy_short_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_SCA_simplified_sqrt <- lm(Canopy_short_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LC_multiple_lm_SCA_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we could use the canopy short variale either with log or square root transformation


#checking equal variance
ggplot(data = LC_multiple_lm_SCA_simplified_lg, aes(x = LC_multiple_lm_SCA_simplified_lg$fitted.values, y = LC_multiple_lm_SCA_simplified_lg$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
LC_multiple_lm_SCA_summary <- summary(LC_multiple_lm_SCA)
LC_multiple_lm_SCA_simplified_summary <- summary(LC_multiple_lm_SCA_simplified)
summary(LC_multiple_lm_SCA_simplified_lg)

View(LC_fixed_field_data_processed_terrain_no_NA)


# LCA

plot(LC_multiple_lm_LCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LC_multiple_lm_LCA <- lm(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LC_multiple_lm_LCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

LC_multiple_lm_LCA_summary <- summary(LC_multiple_lm_LCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LC_multiple_lm_LCA_vif <- car::vif(LC_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LC_multiple_lm_LCA_VIF_multi_num <- (1 / (1-LC_multiple_lm_LCA_summary$r.squared))
LC_multiple_lm_LCA_vif > LC_multiple_lm_LCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LC_multiple_lm_LCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LC_multiple_lm_LCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LC_multiple_lm_LCA_simplified <- lm(Canopy_long ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_LCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LC_multiple_lm_LCA_simplified, LC_multiple_lm_LCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LC_potential_interactions_LCA <- rpart(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                         LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LC_potential_interactions_LCA) # show the tree structure
text(LC_potential_interactions_LCA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree, branches mean that the variables likely have interactions with one another
LC_multiple_lm_LCA_interacts <- lm(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                     LC_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LC_slope_raster_15_data_pts +
                                     I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts:LC_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_LCA_interacts)

#slimming down the variables in the interaction model
step(LC_multiple_lm_LCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LC_multiple_lm_LCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_fixed_field_data_processed_terrain_no_NA$I(LC_aspect_raster_15_data_pts_8_categorical^2)


#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_LCA_interacts_simplified_step <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2) + I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts,  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_LCA_interacts_simplified_dredge <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_LCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(LC_multiple_lm_LCA_interacts_simplified_step, LC_multiple_lm_LCA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified


#nested F test to compare simplified interactions model made with dredge to full interactions model
anova(LC_multiple_lm_LCA_interacts_simplified_dredge, LC_multiple_lm_LCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model made with dredge to simplified model without interactions
anova(LC_multiple_lm_LCA_interacts_simplified_dredge, LC_multiple_lm_LCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LC_multiple_lm_LCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_LCA_simplified, aes(x= LC_multiple_lm_LCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_LCA_simplified, aes(sample = LC_multiple_lm_LCA_simplified$residuals))+
  geom_qq()

shapiro.test(LC_multiple_lm_LCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro test was singificant, so I will use a transformaed canopy_lung variable
LC_multiple_lm_LCA_simplified_lg <- lm(Canopy_long_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_LCA_simplified_sqrt <- lm(Canopy_long_sqrt ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LC_multiple_lm_LCA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking equal variance
ggplot(data = LC_multiple_lm_LCA_simplified_sqrt, aes(x = LC_multiple_lm_LCA_simplified_sqrt$fitted.values, y = LC_multiple_lm_LCA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
summary(LC_multiple_lm_LCA)
summary(LC_multiple_lm_LCA_simplified_sqrt) #not significant, the square transformation changes the p-value a lot 
summary(LC_multiple_lm_LCA_simplified) #significant 

# CA

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LC_multiple_lm_CA <- lm(Canopy_area ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LC_multiple_lm_CA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

LC_multiple_lm_CA_summary <- summary(LC_multiple_lm_CA) #storing a summary of the model

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LC_multiple_lm_CA_vif <- car::vif(LC_multiple_lm_CA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LC_multiple_lm_CA_VIF_multi_num <- (1 / (1-LC_multiple_lm_CA_summary$r.squared))
LC_multiple_lm_CA_vif > LC_multiple_lm_CA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LC_multiple_lm_CA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LC_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LC_multiple_lm_CA_simplified <- lm(Canopy_area ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LC_multiple_lm_CA_simplified, LC_multiple_lm_CA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LC_potential_interactions_CA <- rpart(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                        LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LC_potential_interactions_CA) # show the tree structure
text(LC_potential_interactions_CA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LC_multiple_lm_CA_interacts <- lm(Canopy_area ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                    LC_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LC_slope_raster_15_data_pts +
                                    I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts:LC_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CA_interacts)

#slimming down the variables in the interaction model
step(LC_multiple_lm_CA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LC_multiple_lm_CA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_fixed_field_data_processed_terrain_no_NA$I(LC_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_CA_interacts_simplified_step <- lm(Canopy_area ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:LC_slope_raster_15_data_pts,  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LC_fixed_field_data_processed_terrain_no_NA)

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(LC_multiple_lm_CA_interacts_simplified_step, LC_multiple_lm_CA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified

#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CA_interacts_simplified_dredge) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LC_multiple_lm_CA_interacts_simplified_dredge, LC_multiple_lm_CA_interacts) #results are not significant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LC_multiple_lm_CA_interacts_simplified_dredge, LC_multiple_lm_CA_simplified) #results are significant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LC_multiple_lm_CA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_CA_simplified, aes(x= LC_multiple_lm_CA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_CA_simplified, aes(sample = LC_multiple_lm_CA_simplified$residuals))+
  geom_qq()

shapiro.test(LC_multiple_lm_CA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LC_multiple_lm_CA_simplified_lg <- lm(Canopy_area_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_CA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LC_multiple_lm_CA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation

#the results of the Shapiro-Wilk test suggest we should use the model where canopy area is square rooted

#checking equal variance
ggplot(data = LC_multiple_lm_CA_simplified_sqrt, aes(x = LC_multiple_lm_CA_simplified_sqrt$fitted.values, y = LC_multiple_lm_CA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
LC_multiple_lm_CA_summary <- summary(LC_multiple_lm_CA) #sign
LC_multiple_lm_CA_simplified_summary <- summary(LC_multiple_lm_CA_simplified) #sign
LC_multiple_lm_CA_simplified_sqrt_summary <- summary(LC_multiple_lm_CA_simplified_sqrt) #sig


# CS
#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LC_multiple_lm_CS <- lm(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LC_multiple_lm_CS) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing a summary of the model
LC_multiple_lm_CS_summary <- summary(LC_multiple_lm_CS)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LC_multiple_lm_CS_vif <- car::vif(LC_multiple_lm_CS) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LC_multiple_lm_CS_VIF_multi_num <- (1 / (1-LC_multiple_lm_CS_summary$r.squared))
LC_multiple_lm_CS_vif > LC_multiple_lm_CS_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LC_multiple_lm_CS) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LC_multiple_lm_CS) #generates the best model and the rank of best models

#both the step and dredge technique produced the same simplified model:
LC_multiple_lm_CS_simplified <- lm(Crown_spread ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CS_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LC_multiple_lm_CS_simplified, LC_multiple_lm_CS) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#best simplified model without taking into account interactions: LC_multiple_lm_CS_simplified

#determing interactions with recursive binary partioning and regression tree
LC_potential_interactions_CS <- rpart(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                        LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LC_potential_interactions_CS) # show the tree structure
text(LC_potential_interactions_CS, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LC_multiple_lm_CS_interacts <- lm(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                    LC_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LC_slope_raster_15_data_pts +
                                    I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts:LC_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CS_interacts)

#slimming down the variables in the interaction model
step(LC_multiple_lm_CS_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LC_multiple_lm_CS_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_fixed_field_data_processed_terrain_no_NA$I(LC_aspect_raster_15_data_pts_8_categorical^2)

#the step and dredge methods produced the sample best model:
#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_CS_interacts_simplified <- lm(Crown_spread ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_CS_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LC_multiple_lm_CS_interacts_simplified, LC_multiple_lm_CS_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LC_multiple_lm_CS_interacts_simplified, LC_multiple_lm_CS_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LC_multiple_lm_CS_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_CS_simplified, aes(x= LC_multiple_lm_CS_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_CS_simplified, aes(sample = LC_multiple_lm_CS_simplified$residuals))+
  geom_qq()

shapiro.test(LC_multiple_lm_CS_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LC_multiple_lm_CS_simplified_lg <- lm(Crown_spread_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_CS_simplified_sqrt <- lm(Crown_spread_sqrt ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LC_multiple_lm_CS_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_CS_simplified_sqrt, aes(x= LC_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for sqrt(Crown Spread) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_CS_simplified_sqrt, aes(sample = LC_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_multiple_lm_CS_simplified_sqrt, aes(x = LC_multiple_lm_CS_simplified_sqrt$fitted.values, y = LC_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
LC_multiple_lm_SCA_summary <- summary(LC_multiple_lm_SCA) #sig
LC_multiple_lm_SCA_simplified_summary <- summary(LC_multiple_lm_SCA_simplified) #sig
LC_multiple_lm_SCA_simplified_sqrt_summary <- summary(LC_multiple_lm_CS_simplified_sqrt) #sig


# DBH_ag

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LC_multiple_lm_DBH <- lm(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(LC_multiple_lm_DBH) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing summary of the model
LC_multiple_lm_DBH_summary <- summary(LC_multiple_lm_DBH)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
LC_multiple_lm_DBH_vif <- car::vif(LC_multiple_lm_DBH) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
LC_multiple_lm_DBH_VIF_multi_num <- (1 / (1-LC_multiple_lm_DBH_summary$r.squared))
LC_multiple_lm_DBH_vif > LC_multiple_lm_DBH_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(LC_multiple_lm_DBH) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(LC_multiple_lm_DBH) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
LC_multiple_lm_DBH_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_DBH_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(LC_multiple_lm_DBH_simplified, LC_multiple_lm_DBH) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
LC_potential_interactions_DBH <- rpart(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                         LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(LC_potential_interactions_DBH) # show the tree structure
text(LC_potential_interactions_DBH, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
LC_multiple_lm_DBH_interacts <- lm(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + 
                                     LC_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:LC_slope_raster_15_data_pts +
                                     I(LC_slope_raster_15_data_pts^2) + LC_slope_raster_15_data_pts:LC_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:LC_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_DBH_interacts)

#slimming down the variables in the interaction model
step(LC_multiple_lm_DBH_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(LC_multiple_lm_DBH_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_fixed_field_data_processed_terrain_no_NA$I(LC_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
LC_multiple_lm_DBH_interacts_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_multiple_lm_DBH_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(LC_multiple_lm_DBH_interacts_simplified, LC_multiple_lm_DBH_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(LC_multiple_lm_DBH_interacts_simplified, LC_multiple_lm_DBH_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that the simplified interaction and simplified regular model are the same, so we could use either
# for the now we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: LC_multiple_lm_DBH_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_DBH_simplified, aes(x= LC_multiple_lm_DBH_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_DBH_simplified, aes(sample = LC_multiple_lm_DBH_simplified$residuals))+
  geom_qq()

shapiro.test(LC_multiple_lm_DBH_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
LC_multiple_lm_DBH_simplified_lg <- lm(DBH_ag_lg ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)
LC_multiple_lm_DBH_simplified_sqrt <- lm(DBH_ag_sqrt ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed_terrain_no_NA)

shapiro.test(LC_multiple_lm_DBH_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long varaible that has a logged transformation

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_multiple_lm_DBH_simplified_lg, aes(x= LC_multiple_lm_DBH_simplified_lg$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for log(DBH) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_multiple_lm_DBH_simplified_lg, aes(sample = LC_multiple_lm_DBH_simplified_lg$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LC_multiple_lm_DBH_simplified_lg, aes(x = LC_multiple_lm_DBH_simplified_lg$fitted.values, y = LC_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for log(DBH) and Elevation")

#extracting model characteristics and significant
LC_multiple_lm_DBH_summary <- summary(LC_multiple_lm_DBH) #not sig
LC_multiple_lm_DBH_simplified_summary <- summary(LC_multiple_lm_DBH_simplified) #sig
LC_multiple_lm_DBH_simplified_lg_summary <- summary(LC_multiple_lm_DBH_simplified_lg) #sig





# SD


#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
SD_fixed_field_data_processed_terrain_no_NA <- SD_fixed_field_data_processed_terrain %>%
  filter(is.na(SD_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)


# SCA

plot(SD_multiple_lm_SCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
LM_multiple_lm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(SD_multiple_lm_SCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing the summary of the model
LM_multiple_lm_SCA_summary <- summary(SD_multiple_lm_SCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
SD_multiple_lm_SCA_vif <- car::vif(SD_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
SD_multiple_lm_SCA_VIF_multi_num <- (1 / (1-SD_multiple_lm_SCA_summary$r.squared))
SD_multiple_lm_SCA_vif > SD_multiple_lm_SCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(SD_multiple_lm_SCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(SD_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
SD_multiple_lm_SCA_simplified <- lm(Canopy_short ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_SCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(SD_multiple_lm_SCA_simplified, SD_multiple_lm_SCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
SD_potential_interactions_SCA <- rpart(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                         SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(SD_potential_interactions) # show the tree structure
text(SD_potential_interactions, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
SD_multiple_lm_SCA_interacts <- lm(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                     SD_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:SD_slope_raster_15_data_pts +
                                     I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts:SD_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_SCA_interacts)

#slimming down the variables in the interaction model
step(SD_multiple_lm_SCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(SD_multiple_lm_SCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_fixed_field_data_processed_terrain_no_NA$I(SD_aspect_raster_15_data_pts_8_categorical^2)


#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_SCA_interacts_simplified <- lm(Canopy_short ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_SCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(SD_multiple_lm_SCA_interacts_simplified, SD_multiple_lm_SCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(SD_multiple_lm_SCA_interacts_simplified, SD_multiple_lm_SCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: SD_multiple_lm_SCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_SCA_simplified, aes(x= SD_multiple_lm_SCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_SCA_simplified, aes(sample = SD_multiple_lm_SCA_simplified$residuals))+
  geom_qq()

shapiro.test(SD_multiple_lm_SCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#checking equal variance
ggplot(data = SD_multiple_lm_SCA, aes(x = SD_multiple_lm_SCA$fitted.values, y = SD_multiple_lm_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
SD_multiple_lm_SCA_summary <- summary(SD_multiple_lm_SCA)
SD_multiple_lm_SCA_simplified_summary <- summary(SD_multiple_lm_SCA_simplified)

View(SD_fixed_field_data_processed_terrain_no_NA)


# LCA

plot(SD_multiple_lm_LCA)

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
SD_multiple_lm_LCA <- lm(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(SD_multiple_lm_LCA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

SD_multiple_lm_LCA_summary <- summary(SD_multiple_lm_LCA)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
SD_multiple_lm_LCA_vif <- car::vif(SD_multiple_lm_SCA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
SD_multiple_lm_LCA_VIF_multi_num <- (1 / (1-SD_multiple_lm_LCA_summary$r.squared))
SD_multiple_lm_LCA_vif > SD_multiple_lm_LCA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(SD_multiple_lm_LCA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(SD_multiple_lm_LCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
SD_multiple_lm_LCA_simplified <- lm(Canopy_long ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_LCA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(SD_multiple_lm_LCA_simplified, SD_multiple_lm_LCA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
SD_potential_interactions_LCA <- rpart(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                         SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(SD_potential_interactions_LCA) # show the tree structure
text(SD_potential_interactions_LCA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree, branches mean that the variables likely have interactions with one another
SD_multiple_lm_LCA_interacts <- lm(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                     SD_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:SD_slope_raster_15_data_pts +
                                     I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts:SD_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_LCA_interacts)

#slimming down the variables in the interaction model
step(SD_multiple_lm_LCA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(SD_multiple_lm_LCA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_fixed_field_data_processed_terrain_no_NA$I(SD_aspect_raster_15_data_pts_8_categorical^2)


#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_LCA_interacts_simplified_step <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2) + I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts,  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_LCA_interacts_simplified_dredge <- lm(Canopy_long ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_LCA_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(SD_multiple_lm_LCA_interacts_simplified_step, SD_multiple_lm_LCA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified


#nested F test to compare simplified interactions model made with dredge to full interactions model
anova(SD_multiple_lm_LCA_interacts_simplified_dredge, SD_multiple_lm_LCA_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model made with dredge to simplified model without interactions
anova(SD_multiple_lm_LCA_interacts_simplified_dredge, SD_multiple_lm_LCA_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: SD_multiple_lm_LCA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_LCA_simplified, aes(x= SD_multiple_lm_LCA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_LCA_simplified, aes(sample = SD_multiple_lm_LCA_simplified$residuals))+
  geom_qq()

shapiro.test(SD_multiple_lm_LCA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro test was singificant, so I will use a transformaed canopy_lung variable
SD_multiple_lm_LCA_simplified_lg <- lm(Canopy_long_lg ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_LCA_simplified_sqrt <- lm(Canopy_long_sqrt ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)

shapiro.test(SD_multiple_lm_LCA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking equal variance
ggplot(data = SD_multiple_lm_LCA_simplified_sqrt, aes(x = SD_multiple_lm_LCA_simplified_sqrt$fitted.values, y = SD_multiple_lm_LCA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
SD_multiple_lm_LCA_summary <- summary(SD_multiple_lm_LCA)
SD_multiple_lm_LCA_simplified_summary <- summary(SD_multiple_lm_LCA_simplified_sqrt) #not significant, the square transformation changes the p-value a lot 
SD_multiple_lm_LCA_simplified_summary <- summary(SD_multiple_lm_LCA_simplified) #significant 

# CA

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
SD_multiple_lm_CA <- lm(Canopy_area ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(SD_multiple_lm_CA) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

SD_multiple_lm_CA_summary <- summary(SD_multiple_lm_CA) #storing a summary of the model

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
SD_multiple_lm_CA_vif <- car::vif(SD_multiple_lm_CA) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
SD_multiple_lm_CA_VIF_multi_num <- (1 / (1-SD_multiple_lm_CA_summary$r.squared))
SD_multiple_lm_CA_vif > SD_multiple_lm_CA_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(SD_multiple_lm_CA) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(SD_multiple_lm_SCA) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
SD_multiple_lm_CA_simplified <- lm(Canopy_area ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CA_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(SD_multiple_lm_CA_simplified, SD_multiple_lm_CA) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
SD_potential_interactions_CA <- rpart(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                        SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(SD_potential_interactions_CA) # show the tree structure
text(SD_potential_interactions_CA, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
SD_multiple_lm_CA_interacts <- lm(Canopy_area ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                    SD_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:SD_slope_raster_15_data_pts +
                                    I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts:SD_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CA_interacts)

#slimming down the variables in the interaction model
step(SD_multiple_lm_CA_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(SD_multiple_lm_CA_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_fixed_field_data_processed_terrain_no_NA$I(SD_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_CA_interacts_simplified_step <- lm(Canopy_area ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + I(Elevation..m.FIXED^2) + Elevation..m.FIXED:SD_slope_raster_15_data_pts,  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = SD_fixed_field_data_processed_terrain_no_NA)

#nested F test to compare simplified interactions model using step and simplified interactions model using dredge
anova(SD_multiple_lm_CA_interacts_simplified_step, SD_multiple_lm_CA_interacts_simplified_dredge) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one
#because the dredge and step models are not signficiantly different, I will be using the dredge one because it is more simplified

#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_CA_interacts_simplified_dredge <- lm(Canopy_area ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2), data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CA_interacts_simplified_dredge) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(SD_multiple_lm_CA_interacts_simplified_dredge, SD_multiple_lm_CA_interacts) #results are not significant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(SD_multiple_lm_CA_interacts_simplified_dredge, SD_multiple_lm_CA_simplified) #results are significant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: SD_multiple_lm_CA_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_CA_simplified, aes(x= SD_multiple_lm_CA_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_CA_simplified, aes(sample = SD_multiple_lm_CA_simplified$residuals))+
  geom_qq()

shapiro.test(SD_multiple_lm_CA_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
SD_multiple_lm_CA_simplified_lg <- lm(Canopy_area_lg ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_CA_simplified_sqrt <- lm(Canopy_area_sqrt ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)

shapiro.test(SD_multiple_lm_CA_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation

#the results of the Shapiro-Wilk test suggest we should use the model where canopy area is square rooted

#checking equal variance
ggplot(data = SD_multiple_lm_CA_simplified_sqrt, aes(x = SD_multiple_lm_CA_simplified_sqrt$fitted.values, y = SD_multiple_lm_CA_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for sqrt(SCA) and Elevation")


#extracting model characteristics and significant
SD_multiple_lm_CA_summary <- summary(SD_multiple_lm_CA) #sign
SD_multiple_lm_CA_simplified_summary <- summary(SD_multiple_lm_CA_simplified) #sign
SD_multiple_lm_CA_simplified_sqrt_summary <- summary(SD_multiple_lm_CA_simplified_sqrt) #sig


# CS
#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
SD_multiple_lm_CS <- lm(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(SD_multiple_lm_CS) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing a summary of the model
SD_multiple_lm_CS_summary <- summary(SD_multiple_lm_CS)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
SD_multiple_lm_CS_vif <- car::vif(SD_multiple_lm_CS) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
SD_multiple_lm_CS_VIF_multi_num <- (1 / (1-SD_multiple_lm_CS_summary$r.squared))
SD_multiple_lm_CS_vif > SD_multiple_lm_CS_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(SD_multiple_lm_CS) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(SD_multiple_lm_CS) #generates the best model and the rank of best models

#both the step and dredge technique produced the same simplified model:
SD_multiple_lm_CS_simplified <- lm(Crown_spread ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CS_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(SD_multiple_lm_CS_simplified, SD_multiple_lm_CS) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#best simplified model without taking into account interactions: SD_multiple_lm_CS_simplified

#determing interactions with recursive binary partioning and regression tree
SD_potential_interactions_CS <- rpart(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                        SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(SD_potential_interactions_CS) # show the tree structure
text(SD_potential_interactions_CS, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
SD_multiple_lm_CS_interacts <- lm(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                    SD_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:SD_slope_raster_15_data_pts +
                                    I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts:SD_aspect_raster_15_data_pts_8_categorical +
                                    Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                  data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CS_interacts)

#slimming down the variables in the interaction model
step(SD_multiple_lm_CS_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(SD_multiple_lm_CS_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_fixed_field_data_processed_terrain_no_NA$I(SD_aspect_raster_15_data_pts_8_categorical^2)

#the step and dredge methods produced the sample best model:
#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_CS_interacts_simplified <- lm(Crown_spread ~ Elevation..m.FIXED + I(Elevation..m.FIXED^2),  data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_CS_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(SD_multiple_lm_CS_interacts_simplified, SD_multiple_lm_CS_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(SD_multiple_lm_CS_interacts_simplified, SD_multiple_lm_CS_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: SD_multiple_lm_CS_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_CS_simplified, aes(x= SD_multiple_lm_CS_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation + Slope + Aspect")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_CS_simplified, aes(sample = SD_multiple_lm_CS_simplified$residuals))+
  geom_qq()

shapiro.test(SD_multiple_lm_CS_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
SD_multiple_lm_CS_simplified_lg <- lm(Crown_spread_lg ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_CS_simplified_sqrt <- lm(Crown_spread_sqrt ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)

shapiro.test(SD_multiple_lm_CS_simplified_sqrt$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long vaaraible that has a square root transformation


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_CS_simplified_sqrt, aes(x= SD_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for sqrt(Crown Spread) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_CS_simplified_sqrt, aes(sample = SD_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_multiple_lm_CS_simplified_sqrt, aes(x = SD_multiple_lm_CS_simplified_sqrt$fitted.values, y = SD_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation + Slope + Aspect")


#extracting model characteristics and significant
SD_multiple_lm_SCA_summary <- summary(SD_multiple_lm_SCA) #sig
SD_multiple_lm_SCA_simplified_summary <- summary(SD_multiple_lm_SCA_simplified) #sig
SD_multiple_lm_SCA_simplified_sqrt_summary <- summary(SD_multiple_lm_CS_simplified_sqrt) #sig


# DBH_ag

#multiple linear regression base model with all variables, and using the no NA dataset to be able to use the backwards regression
SD_multiple_lm_DBH <- lm(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)

#checking to see which variables might be the most useful
avPlots(SD_multiple_lm_DBH) #added variable plots, looking to see which variables might be most useful in exlaining the size/shape variables 

#storing summary of the model
SD_multiple_lm_DBH_summary <- summary(SD_multiple_lm_DBH)

#checking for any multicollinarity, all of them have great than 1/(1-r^2)  VIF values, meaning there is multicollinarity
SD_multiple_lm_DBH_vif <- car::vif(SD_multiple_lm_DBH) #variance inflation factor, looking for if values is greater than 5 or 10, or if  If the VIF is larger than 1/(1-R2), where R2 is the Multiple R-squared of the regression, then that predictor is more related to the other predictors than it is to the response.
SD_multiple_lm_DBH_VIF_multi_num <- (1 / (1-SD_multiple_lm_DBH_summary$r.squared))
SD_multiple_lm_DBH_vif > SD_multiple_lm_DBH_VIF_multi_num

#determinging our main effects model with two methods: backward's regression and the dredge function 
step(SD_multiple_lm_DBH) #using backwards regression, where last model produced is the best fit

options(na.action = "na.fail") #have to set na.action to na.fail to be able to run dredge
dredge(SD_multiple_lm_DBH) #generates the best model and the rank of best models

#the best simplified multiple linear regression model chosen
SD_multiple_lm_DBH_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_DBH_simplified) #best model, but still only 5% of variability explained

#nested F test comparing the simplified model to the original, If model 1 is really correct, what is the chance that you would randomly obtain data that fits model 2 so much better?
anova(SD_multiple_lm_DBH_simplified, SD_multiple_lm_DBH) #results are not signfiicant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#determing interactions with recursive binary partioning and regression tree
SD_potential_interactions_DBH <- rpart(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                         SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain_no_NA)
par(xpd = TRUE) # allows text to "eXPanD" (spill over outside the plotting area)
plot(SD_potential_interactions_DBH) # show the tree structure
text(SD_potential_interactions_DBH, pretty = 0) # add text labels

#there does appear to be interactions, so we must make an interactions model

#interactions model, based on results of regression tree
SD_multiple_lm_DBH_interacts <- lm(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + 
                                     SD_aspect_raster_15_data_pts_8_categorical + Elevation..m.FIXED:SD_slope_raster_15_data_pts +
                                     I(SD_slope_raster_15_data_pts^2) + SD_slope_raster_15_data_pts:SD_aspect_raster_15_data_pts_8_categorical +
                                     Elevation..m.FIXED:SD_aspect_raster_15_data_pts_8_categorical + I(Elevation..m.FIXED^2), 
                                   data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_DBH_interacts)

#slimming down the variables in the interaction model
step(SD_multiple_lm_DBH_interacts) #using backwards regression, where last model produced is the best fit
dredge <- dredge(SD_multiple_lm_DBH_interacts) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_fixed_field_data_processed_terrain_no_NA$I(SD_aspect_raster_15_data_pts_8_categorical^2)

#including interactions, the best simplified multiple linear regression model chosen
SD_multiple_lm_DBH_interacts_simplified <- lm(DBH_ag ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_multiple_lm_DBH_interacts_simplified) #best model, but still only 5% of variability explained

#nested F test to compare simplified interactions model to full interactions model
anova(SD_multiple_lm_DBH_interacts_simplified, SD_multiple_lm_DBH_interacts) #results are not signficant, meaning there is no compelling evidence to support the larger model and we should stick with the smaller one

#nested F test to compare simplified interactions model to simplified model without interactions
anova(SD_multiple_lm_DBH_interacts_simplified, SD_multiple_lm_DBH_simplified) #results are signficant, meaning there is compelling evidence to support the smaller model than the larger one

# our results indicate that the simplified interaction and simplified regular model are the same, so we could use either
# for the now we should: use the multiple linear regression that is simplified over the full model and the models that include interactions

# Best Model: SD_multiple_lm_DBH_simplified

#the model must satisfy LINES (linearity, independence, normality of residuals, equal variance of residuals, and simple random sample)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_DBH_simplified, aes(x= SD_multiple_lm_DBH_simplified$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

theme(axis.title.x = element_text(size=15),
      axis.title.y = element_text(size=15),
      axis.text.x = element_text(size=12), 
      axis.text.y = element_text(size=12))

#qqnorm plot
ggplot(SD_multiple_lm_DBH_simplified, aes(sample = SD_multiple_lm_DBH_simplified$residuals))+
  geom_qq()

shapiro.test(SD_multiple_lm_DBH_simplified$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed

#shapiro-wilk test is significant, so we will use a model where canopy area is transformed
SD_multiple_lm_DBH_simplified_lg <- lm(DBH_ag_lg ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)
SD_multiple_lm_DBH_simplified_sqrt <- lm(DBH_ag_sqrt ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed_terrain_no_NA)

shapiro.test(SD_multiple_lm_DBH_simplified_lg$residuals) #shapiro welk test for normality, if significant, then the residuals are not likely normally distributed
#based on the Shapiro-Wilk test we should used the canopy long varaible that has a logged transformation

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_multiple_lm_DBH_simplified_lg, aes(x= SD_multiple_lm_DBH_simplified_lg$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for log(DBH) vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_multiple_lm_DBH_simplified_lg, aes(sample = SD_multiple_lm_DBH_simplified_lg$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_multiple_lm_DBH_simplified_lg, aes(x = SD_multiple_lm_DBH_simplified_lg$fitted.values, y = SD_multiple_lm_CS_simplified_sqrt$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for log(DBH) and Elevation")

#extracting model characteristics and significant
SD_multiple_lm_DBH_summary <- summary(SD_multiple_lm_DBH) #not sig
SD_multiple_lm_DBH_simplified_summary <- summary(SD_multiple_lm_DBH_simplified) #sig
SD_multiple_lm_DBH_simplified_lg_summary <- summary(SD_multiple_lm_DBH_simplified_lg) #sig
