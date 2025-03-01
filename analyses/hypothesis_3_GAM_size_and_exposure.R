#After having trouble with multiple linear rgressions because of issues 
#with the normality conditions and some nervousness about linearity. 
#We decided to use Generalized Additive Models (GAMs)


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
library(mgcv) #needed for the gam function for generalized additive models
library(starsExtra) #to use dist_to_nearest
library(MuMIn) #to use the dredge function
library(rpart) #to use the function rpart to check recurissive binary
library(visreg) #package to be able to plot effects of categorical variables
library(gratia) #using the function smooth estimates
library(gridExtra) #way to arrange ggplots in one plot
library(plotly) #3d plotting
devtools::install_github("AckerDWM/gg3D") #3d plotting
library("gg3D") #3d plotting
library(mgcViz) #3d plotting
library(rgl) #3d plotting



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
fixed_field_data_processed_sf_trans_coordinates_download <- write.csv(fixed_field_data_processed_sf_trans_coordinates, "/Users/chewbecca/Morton Arboretum REU 2024/Untitled/QUBR_GenGeoEcoDemoCorr/data/fixed_field_data_processed_sf_trans_coordinates.csv", row.names = F)
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

View(fixed_field_data_processed_sf_trans_coordinates)
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
  DBHle_fill_viridis_c()


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


View(LM_fixed_field_data_processed_terrain)


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
                                                                        (all_points_aspect_raster_15_data_pts >= 337.5 & all_points_aspect_raster_15_data_pts < 359.99999) ~ "N",
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
                                                                        (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


# LM

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_8_categorical = case_when((LM_aspect_raster_15_data_pts > 0 & LM_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LM_aspect_raster_15_data_pts >= 337.5 & LM_aspect_raster_15_data_pts < 359.99999) ~ "N",
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
                                                                (LM_aspect_raster_15_data_pts >= 315 & LM_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LM_aspect_raster_15_data_pts >= 135 & LM_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LM_aspect_raster_15_data_pts >= 225 & LM_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


# LC

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_8_categorical = case_when((LC_aspect_raster_15_data_pts > 0 & LC_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 337.5 & LC_aspect_raster_15_data_pts < 359.99999) ~ "N",
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
                                                                (LC_aspect_raster_15_data_pts >= 315 & LC_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                (LC_aspect_raster_15_data_pts >= 45 & LC_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LC_aspect_raster_15_data_pts >= 135 & LC_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LC_aspect_raster_15_data_pts >= 225 & LC_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315



#SD
# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_8_categorical = case_when((SD_aspect_raster_15_data_pts > 0 & SD_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 337.5 & SD_aspect_raster_15_data_pts < 359.99999) ~ "N",
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
                                                                (SD_aspect_raster_15_data_pts >= 315 & SD_aspect_raster_15_data_pts < 359.99999) ~ "N",
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


#descriptive summary for LM


### Generalized Additive Models ###

#using only the 8 categories

# all points 

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
all_points_fixed_field_data_processed_terrain_no_NA <- all_points_fixed_field_data_processed_terrain %>%
  filter(is.na(all_points_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(all_points_aspect_raster_15_data_pts_8_categorical) == F)

#Cook's D
plot(all_points_multiple_lm_SCA)
all_points_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_mlm_SCA_cooks <- cooks.distance(all_points_mlm_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (2 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential


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

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_SCA, all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed_first_term, 
    all_points_add.gam_SCA.smoothed_second_term, all_points_add.gam_SCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_SCA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_SCA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_SCA.smoothed.quasi, all_points_add.gam_SCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_SCA.smoothed.quasi, all_points_add.gam_SCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_SCA.smoothed.poisson, all_points_add.gam_SCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_SCA)
summary(all_points_add.gam_SCA.smoothed)
summary(all_points_add.gam_SCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_SCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_SCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_SCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  #geom_smooth(se = T) +
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

p3 <- visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
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

# LCA

all_points_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_LCA_interact <- gam(Canopy_long ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_LCA, all_points_add.gam_LCA.smoothed, all_points_add.gam_LCA.smoothed_first_term, 
    all_points_add.gam_LCA.smoothed_second_term, all_points_add.gam_LCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_LCA.smoothed.quasi <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_LCA.smoothed.quasipoisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_LCA.smoothed, all_points_add.gam_LCA.smoothed.quasi, test = "LRT") #smooth vs. quasi
anova(all_points_add.gam_LCA.smoothed.quasi, all_points_add.gam_LCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_LCA.smoothed.quasi, all_points_add.gam_LCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_LCA.smoothed.poisson, all_points_add.gam_LCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a quasipoisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_LCA.smoothed.quasipoisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_LCA)
summary(all_points_add.gam_LCA.smoothed)
summary(all_points_add.gam_LCA.smoothed.poisson)
summary(all_points_add.gam_LCA.smoothed.quasipoisson)

#we do not need to dredge the poisson model, but hear is the 
options(na.action = "na.fail")
dredge <- dredge(all_points_add.gam_LCA.smoothed.quasipoisson, rank = "") #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_LCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_LCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_LCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_LCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_LCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_LCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_LCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_LCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_LCA.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(all_points_add.gam_LCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(all_points_add.gam_LCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_LCA <- fitted.values(all_points_add.gam_LCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                                                                z=Canopy_long, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$Canopy_long ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)



# CA

all_points_add.gam_CA <- gam(Canopy_area ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CA.smoothed_first_term <- gam(Canopy_area ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CA.smoothed_second_term <- gam(Canopy_area ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CA_interact <- gam(Canopy_area ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_CA, all_points_add.gam_CA.smoothed, all_points_add.gam_CA.smoothed_first_term, 
    all_points_add.gam_CA.smoothed_second_term, all_points_add.gam_CA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_CA.smoothed.quasi <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_CA.smoothed.quasipoisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_CA.smoothed, all_points_add.gam_CA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_CA.smoothed.quasi, all_points_add.gam_CA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_CA.smoothed.quasi, all_points_add.gam_CA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_CA.smoothed.poisson, all_points_add.gam_CA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_CA)
summary(all_points_add.gam_CA.smoothed)
summary(all_points_add.gam_CA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_CA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_SCA.smoothed.dredged <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_CA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_CA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_CA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_CA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_CA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_CA.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Canopy Area", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Canopy Area", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(all_points_add.gam_CA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(all_points_add.gam_CA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_CA <- fitted.values(all_points_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                                                                z=Canopy_area, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$Canopy_area ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


# CS

all_points_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_CS_interact <- gam(Crown_spread ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_CS, all_points_add.gam_CS.smoothed, all_points_add.gam_CS.smoothed_first_term, 
    all_points_add.gam_CS.smoothed_second_term, all_points_add.gam_CS_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam.CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_CS.smoothed.quasi <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_CS.smoothed.quasipoisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_CS.smoothed, all_points_add.gam_CS.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_CS.smoothed.quasi, all_points_add.gam_CS.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_CS.smoothed.quasi, all_points_add.gam_CS.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_CS.smoothed.poisson, all_points_add.gam_CS.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CS.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_CS)
summary(all_points_add.gam_CS.smoothed)
summary(all_points_add.gam_CS.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_CS.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_CS.smoothed.dredged <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_CS.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_CS.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_CS.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_CS.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_CS.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_CS.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_CS.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Crown Spread", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Crown Spread", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(all_points_add.gam_CS.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_CS <- fitted.values(all_points_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                                                                z=Crown_spread, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$Crown_spread ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


# DBH_ag

all_points_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_DBH_interact <- gam(DBH_ag ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_DBH, all_points_add.gam_DBH.smoothed, all_points_add.gam_DBH.smoothed_first_term, 
    all_points_add.gam_DBH.smoothed_second_term, all_points_add.gam_DBH_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_DBH.smoothed.quasi <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_DBH.smoothed.quasipoisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_DBH.smoothed, all_points_add.gam_DBH.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_DBH.smoothed.quasi, all_points_add.gam_DBH.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_DBH.smoothed.quasi, all_points_add.gam_DBH.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_DBH.smoothed.poisson, all_points_add.gam_DBH.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_DBH.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_DBH)
summary(all_points_add.gam_DBH.smoothed)
summary(all_points_add.gam_DBH.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_DBH.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_DBH.smoothed.dredged <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_DBH.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_DBH.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_DBH.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_DBH.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_DBH.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_DBH.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on DBH", title = "Smooth Effect of DBH") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on DBH", title = "Smooth Effect of DBH") +
  theme_minimal()

p3 <- visreg(all_points_add.gam_DBH.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(all_points_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_SCA <- fitted.values(all_points_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                                                                z=DBH_ag, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$DBH_ag ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


# LM

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
LM_fixed_field_data_processed_terrain_no_NA <- LM_fixed_field_data_processed_terrain %>%
  filter(is.na(LM_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)

#Cook's D
LM_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_mlm_SCA_cooks <- cooks.distance(LM_mlm_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_terrain_no_NA_No_outliers <- LM_fixed_field_data_processed_terrain_no_NA[-c(24,26,27),]


# SCA

LM_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + LCA_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_SCA_interact <- gam(Canopy_short ~ Elevation..m.FIXED * LM_slope_raster_15_data_pts * LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_SCA, LM_add.gam_SCA.smoothed, LM_add.gam_SCA.smoothed_first_term, 
    LM_add.gam_SCA.smoothed_second_term, LM_add.gam_SCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LCA_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LM_add.gam_SCA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LCA_slope_raster_15_data_pts) + LCA_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasi())
LM_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LCA_slope_raster_15_data_pts) + LCA_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = poisson())
LM_add.gam_SCA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LCA_slope_raster_15_data_pts) + LCA_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LM_add.gam_SCA.smoothed, LM_add.gam_SCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_SCA.smoothed.quasi, LM_add.gam_SCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_SCA.smoothed.quasi, LM_add.gam_SCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LM_add.gam_SCA.smoothed.poisson, LM_add.gam_SCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_SCA)
summary(LM_add.gam_SCA.smoothed)
summary(LM_add.gam_SCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_SCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#Chosen model: LM_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LM_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
k.check(LM_add.gam_SCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LM_add.gam_SCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LM_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LM_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LM_add.gam_SCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LM_add.gam_SCA.smoothed, select = "s(LM_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LM_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LM_add.gam_SCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_short, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LM_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LM_add.gam_SCA <- fitted.values(LM_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LM_fixed_field_data_processed_terrain_no_NA_No_outliers, aes(x=Elevation..m.FIXED, y=LM_slope_raster_15_data_pts, 
                                                                z=Canopy_short, color=LM_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_short ~ 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


# LCA


LM_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_LCA_interact <- gam(Canopy_long ~ Elevation..m.FIXED * LM_slope_raster_15_data_pts * LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_LCA, LM_add.gam_LCA.smoothed, LM_add.gam_LCA.smoothed_first_term, 
    LM_add.gam_LCA.smoothed_second_term, LM_add.gam_LCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LM_add.gam_LCA.smoothed.quasi <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasi())
LM_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = poisson())
LM_add.gam_LCA.smoothed.quasipoisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LM_add.gam_LCA.smoothed, LM_add.gam_LCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_LCA.smoothed.quasi, LM_add.gam_LCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_LCA.smoothed.quasi, LM_add.gam_LCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LM_add.gam_LCA.smoothed.poisson, LM_add.gam_LCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_LCA)
summary(LM_add.gam_LCA.smoothed)
summary(LM_add.gam_LCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_LCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_add.gam_LCA.smoothed.dredged <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#Chosen model: LM_add.gam_LCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LM_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
k.check(LM_add.gam_LCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LM_add.gam_LCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LM_add.gam_LCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LM_add.gam_LCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LM_add.gam_LCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LM_add.gam_LCA.smoothed, select = "s(LM_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LM_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LM_add.gam_LCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_short, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LM_add.gam_LCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LM_add.gam_LCA <- fitted.values(LM_add.gam_LCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LM_fixed_field_data_processed_terrain_no_NA_No_outliers, aes(x=Elevation..m.FIXED, y=LM_slope_raster_15_data_pts, 
                                                                z=Canopy_long, color=LM_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_long ~ 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


# CA
LM_add.gam_CA <- gam(Canopy_area ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CA.smoothed_first_term <- gam(Canopy_area ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CA.smoothed_second_term <- gam(Canopy_area ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CA_interact <- gam(Canopy_area ~ Elevation..m.FIXED * LM_slope_raster_15_data_pts * LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_CA, LM_add.gam_CA.smoothed, LM_add.gam_CA.smoothed_first_term, 
    LM_add.gam_CA.smoothed_second_term, LM_add.gam_CA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LM_add.gam_CA.smoothed.quasi <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasi())
LM_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = poisson())
LM_add.gam_CA.smoothed.quasipoisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LM_add.gam_CA.smoothed, LM_add.gam_CA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_CA.smoothed.quasi, LM_add.gam_CA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_CA.smoothed.quasi, LM_add.gam_CA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LM_add.gam_CA.smoothed.poisson, LM_add.gam_CA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a quasipoisson model is the best fit

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.smoothed.quasipoisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_CA)
summary(LM_add.gam_CA.smoothed)
summary(LM_add.gam_CA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_CA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_add.gam_CA.smoothed.dredged <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#Chosen model: LM_add.gam_CA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LM_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
k.check(LM_add.gam_CA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LM_add.gam_CA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LM_add.gam_CA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LM_add.gam_CA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LM_add.gam_CA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LM_add.gam_CA.smoothed, select = "s(LM_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LM_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LM_add.gam_CA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_area, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LM_add.gam_CA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LM_add.gam_CA <- fitted.values(LM_add.gam_CA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LM_fixed_field_data_processed_terrain_no_NA_No_outliers, aes(x=Elevation..m.FIXED, y=LM_slope_raster_15_data_pts, 
                                                                z=Canopy_area, color=LM_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Canopy_area ~ 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


# CS

LM_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_CS_interact <- gam(Crown_spread ~ Elevation..m.FIXED * LM_slope_raster_15_data_pts * LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_CS, LM_add.gam_CS.smoothed, LM_add.gam_CS.smoothed_first_term, 
    LM_add.gam_CS.smoothed_second_term, LM_add.gam_CS_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LM_add.gam_CS.smoothed.quasi <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasi())
LM_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = poisson())
LM_add.gam_CS.smoothed.quasipoisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LM_add.gam_CS.smoothed, LM_add.gam_CS.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_CS.smoothed.quasi, LM_add.gam_CS.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_CS.smoothed.quasi, LM_add.gam_CS.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LM_add.gam_CS.smoothed.poisson, LM_add.gam_CS.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_CS)
summary(LM_add.gam_CS.smoothed)
summary(LM_add.gam_CS.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_CS.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_add.gam_CS.smoothed.dredged <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#Chosen model: LM_add.gam_CS.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LM_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
k.check(LM_add.gam_CS.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LM_add.gam_CS.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LM_add.gam_CS.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LM_add.gam_CS.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LM_add.gam_CS.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LM_add.gam_CS.smoothed, select = "s(LM_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Crown Spread", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LM_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Crown Spread", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LM_add.gam_CS.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Crown_spread, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LM_add.gam_CS.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LM_add.gam_CS <- fitted.values(LM_add.gam_CS.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LM_fixed_field_data_processed_terrain_no_NA_No_outliers, aes(x=Elevation..m.FIXED, y=LM_slope_raster_15_data_pts, 
                                                                z=Crown_spread, color=LM_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Crown_spread ~ 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)

# DBH_ag

LM_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
LM_add.gam_DBH_interact <- gam(DBH_ag ~ Elevation..m.FIXED * LM_slope_raster_15_data_pts * LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_DBH, LM_add.gam_DBH.smoothed, LM_add.gam_DBH.smoothed_first_term, 
    LM_add.gam_DBH.smoothed_second_term, LM_add.gam_DBH_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LM_add.gam_DBH.smoothed.quasi <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasi())
LM_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = poisson())
LM_add.gam_DBH.smoothed.quasipoisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LM_add.gam_DBH.smoothed, LM_add.gam_DBH.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_DBH.smoothed.quasi, LM_add.gam_DBH.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LM_add.gam_DBH.smoothed.quasi, LM_add.gam_DBH.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LM_add.gam_DBH.smoothed.poisson, LM_add.gam_DBH.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.smoothed.poisson)
gam.check(LM_add.gam_DBH.smoothed)
# based on the issues with normality, we will not use poisson and go with the regular smoothed model


#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_DBH)
summary(LM_add.gam_DBH.smoothed)
summary(LM_add.gam_DBH.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_DBH.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LM_add.gam_DBH.smoothed.dredged <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)

#Chosen model: LM_add.gam_DBH.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LM_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                               data = LM_fixed_field_data_processed_terrain_no_NA_No_outliers)
k.check(LM_add.gam_DBH.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LM_add.gam_DBH.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LM_add.gam_DBH.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LM_add.gam_DBH.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LM_add.gam_DBH.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LM_add.gam_DBH.smoothed, select = "s(LM_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on DBH", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LM_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on DBH", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LM_add.gam_DBH.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$DBH_ag, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LM_add.gam_DBH.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LM_add.gam_DBH <- fitted.values(LM_add.gam_DBH.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LM_fixed_field_data_processed_terrain_no_NA_No_outliers, aes(x=Elevation..m.FIXED, y=LM_slope_raster_15_data_pts, 
                                                                z=DBH_ag, color=LM_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LM_fixed_field_data_processed_terrain_no_NA_No_outliers$DBH_ag ~ 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$Elevation..m.FIXED + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_slope_raster_15_data_pts + 
                LM_fixed_field_data_processed_terrain_no_NA_No_outliers$LM_aspect_raster_15_data_pts_8_categorical)


# LC

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
LC_fixed_field_data_processed_terrain_no_NA <- LC_fixed_field_data_processed_terrain %>%
  filter(is.na(LC_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)


# SCA

LC_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA_interact <- gam(Canopy_short ~ Elevation..m.FIXED * LC_slope_raster_15_data_pts * LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_SCA, LC_add.gam_SCA.smoothed, LC_add.gam_SCA.smoothed_first_term, 
    LC_add.gam_SCA.smoothed_second_term, LC_add.gam_SCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LC_add.gam_SCA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA, family = quasi())
LC_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA, family = poisson())
LC_add.gam_SCA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LC_add.gam_SCA.smoothed, LC_add.gam_SCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_SCA.smoothed.quasi, LC_add.gam_SCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_SCA.smoothed.quasi, LC_add.gam_SCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LC_add.gam_SCA.smoothed.poisson, LC_add.gam_SCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_SCA)
summary(LC_add.gam_SCA.smoothed)
summary(LC_add.gam_SCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LC_add.gam_SCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)

#Chosen model: LC_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LC_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)
k.check(LC_add.gam_SCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LC_add.gam_SCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LC_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LC_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LC_add.gam_SCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LC_add.gam_SCA.smoothed, select = "s(LC_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LC_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LC_add.gam_SCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LC_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LC_add.gam_SCA <- fitted.values(LC_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LC_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=LC_slope_raster_15_data_pts, 
                                                                z=Canopy_short, color=LC_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LC_fixed_field_data_processed_terrain_no_NA$Canopy_short ~ 
                LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


# LCA


LC_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA_interact <- gam(Canopy_long ~ Elevation..m.FIXED * LC_slope_raster_15_data_pts * LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_LCA, LC_add.gam_LCA.smoothed, LC_add.gam_LCA.smoothed_first_term, 
    LC_add.gam_LCA.smoothed_second_term, LC_add.gam_LCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LC_add.gam_LCA.smoothed.quasi <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA, family = quasi())
LC_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA, family = poisson())
LC_add.gam_LCA.smoothed.quasipoisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LC_add.gam_LCA.smoothed, LC_add.gam_LCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_LCA.smoothed.quasi, LC_add.gam_LCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_LCA.smoothed.quasi, LC_add.gam_LCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LC_add.gam_LCA.smoothed.poisson, LC_add.gam_LCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_LCA)
summary(LC_add.gam_LCA.smoothed)
summary(LC_add.gam_LCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LC_add.gam_LCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_add.gam_LCA.smoothed.dredged <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)

#Chosen model: LC_add.gam_LCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LC_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)
k.check(LC_add.gam_LCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LC_add.gam_LCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LC_add.gam_LCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LC_add.gam_LCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LC_add.gam_LCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LC_add.gam_LCA.smoothed, select = "s(LC_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LC_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LC_add.gam_LCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LC_add.gam_LCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LC_add.gam_LCA <- fitted.values(LC_add.gam_LCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LC_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=LC_slope_raster_15_data_pts, 
                                                                z=Canopy_long, color=LC_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LC_fixed_field_data_processed_terrain_no_NA$Canopy_long ~ 
                LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)



# CA

LC_add.gam_CA <- gam(Canopy_area ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA.smoothed_first_term <- gam(Canopy_area ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA.smoothed_second_term <- gam(Canopy_area ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA_interact <- gam(Canopy_area ~ Elevation..m.FIXED * LC_slope_raster_15_data_pts * LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_CA, LC_add.gam_CA.smoothed, LC_add.gam_CA.smoothed_first_term, 
    LC_add.gam_CA.smoothed_second_term, LC_add.gam_CA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LC_add.gam_CA.smoothed.quasi <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA, family = quasi())
LC_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA, family = poisson())
LC_add.gam_CA.smoothed.quasipoisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LC_add.gam_CA.smoothed, LC_add.gam_CA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_CA.smoothed.quasi, LC_add.gam_CA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_CA.smoothed.quasi, LC_add.gam_CA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LC_add.gam_CA.smoothed.poisson, LC_add.gam_CA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_CA)
summary(LC_add.gam_CA.smoothed)
summary(LC_add.gam_CA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LC_add.gam_CA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_add.gam_CA.smoothed.dredged <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)

#Chosen model: LC_add.gam_CA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LC_add.gam_CA.smoothed.poisson <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)
k.check(LC_add.gam_CA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LC_add.gam_CA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LC_add.gam_CA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LC_add.gam_CA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LC_add.gam_CA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LC_add.gam_CA.smoothed, select = "s(LC_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Canopy Area", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LC_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Canopy Area", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LC_add.gam_CA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LC_add.gam_CA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LC_add.gam_CA <- fitted.values(LC_add.gam_CA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LC_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=LC_slope_raster_15_data_pts, 
                                                                z=Canopy_area, color=LC_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LC_fixed_field_data_processed_terrain_no_NA$Canopy_area ~ 
                LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)



# CS

LC_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CS_interact <- gam(Crown_spread ~ Elevation..m.FIXED * LC_slope_raster_15_data_pts * LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_CS, LC_add.gam_CS.smoothed, LC_add.gam_CS.smoothed_first_term, 
    LC_add.gam_CS.smoothed_second_term, LC_add.gam_CS_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LC_add.gam_CS.smoothed.quasi <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA, family = quasi())
LC_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA, family = poisson())
LC_add.gam_CS.smoothed.quasipoisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LC_add.gam_CS.smoothed, LC_add.gam_CS.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_CS.smoothed.quasi, LC_add.gam_CS.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_CS.smoothed.quasi, LC_add.gam_CS.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LC_add.gam_CS.smoothed.poisson, LC_add.gam_CS.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_CS)
summary(LC_add.gam_CS.smoothed)
summary(LC_add.gam_CS.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LC_add.gam_CS.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_add.gam_CS.smoothed.dredged <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)

#Chosen model: LC_add.gam_CS.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LC_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)
k.check(LC_add.gam_CS.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LC_add.gam_CS.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LC_add.gam_CS.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LC_add.gam_CS.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LC_add.gam_CS.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LC_add.gam_CS.smoothed, select = "s(LC_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Crown Spread", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LC_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Crown Spread", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LC_add.gam_CS.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LC_add.gam_CS.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LC_add.gam_CS <- fitted.values(LC_add.gam_CS.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LC_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=LC_slope_raster_15_data_pts, 
                                                                z=Crown_spread, color=LC_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LC_fixed_field_data_processed_terrain_no_NA$Crown_spread ~ 
                LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)



# DBH_ag


LC_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH_interact <- gam(DBH_ag ~ Elevation..m.FIXED * LC_slope_raster_15_data_pts * LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_DBH, LC_add.gam_DBH.smoothed, LC_add.gam_DBH.smoothed_first_term, 
    LC_add.gam_DBH.smoothed_second_term, LC_add.gam_DBH_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
LC_add.gam_DBH.smoothed.quasi <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA, family = quasi())
LC_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA, family = poisson())
LC_add.gam_DBH.smoothed.quasipoisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(LC_add.gam_DBH.smoothed, LC_add.gam_DBH.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_DBH.smoothed.quasi, LC_add.gam_DBH.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(LC_add.gam_DBH.smoothed.quasi, LC_add.gam_DBH.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(LC_add.gam_DBH.smoothed.poisson, LC_add.gam_DBH.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_DBH)
summary(LC_add.gam_DBH.smoothed)
summary(LC_add.gam_DBH.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LC_add.gam_DBH.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
LC_add.gam_DBH.smoothed.dredged <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)

#Chosen model: LC_add.gam_DBH.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
LC_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA)
k.check(LC_add.gam_DBH.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(LC_add.gam_DBH.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(LC_add.gam_DBH.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(LC_add.gam_DBH.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(LC_add.gam_DBH.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(LC_add.gam_DBH.smoothed, select = "s(LC_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on DBH", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = LC_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on DBH", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(LC_add.gam_DBH.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(LC_add.gam_DBH.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_LC_add.gam_DBH <- fitted.values(LC_add.gam_DBH.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(LC_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=LC_slope_raster_15_data_pts, 
                                                                z=DBH_ag, color=LC_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = LC_fixed_field_data_processed_terrain_no_NA$DBH_ag ~ 
                LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts + 
                LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)



# SD


#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
SD_fixed_field_data_processed_terrain_no_NA <- SD_fixed_field_data_processed_terrain %>%
  filter(is.na(SD_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)


# SCA

SD_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA_interact <- gam(Canopy_short ~ Elevation..m.FIXED * SD_slope_raster_15_data_pts * SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_SCA, SD_add.gam_SCA.smoothed, SD_add.gam_SCA.smoothed_first_term, 
    SD_add.gam_SCA.smoothed_second_term, SD_add.gam_SCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
SD_add.gam_SCA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA, family = quasi())
SD_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA, family = poisson())
SD_add.gam_SCA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                    data = SD_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(SD_add.gam_SCA.smoothed, SD_add.gam_SCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_SCA.smoothed.quasi, SD_add.gam_SCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_SCA.smoothed.quasi, SD_add.gam_SCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(SD_add.gam_SCA.smoothed.poisson, SD_add.gam_SCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_SCA)
summary(SD_add.gam_SCA.smoothed)
summary(SD_add.gam_SCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(SD_add.gam_SCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts), 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)

#Chosen model: SD_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
SD_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)
k.check(SD_add.gam_SCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(SD_add.gam_SCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(SD_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(SD_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(SD_add.gam_SCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(SD_add.gam_SCA.smoothed, select = "s(SD_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = SD_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(SD_add.gam_SCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(SD_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_SD_add.gam_SCA <- fitted.values(SD_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(SD_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=SD_slope_raster_15_data_pts, 
                                                                z=Canopy_short, color=SD_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = SD_fixed_field_data_processed_terrain_no_NA$Canopy_short ~ 
                SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)



# LCA

SD_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_LCA_interact <- gam(Canopy_long ~ Elevation..m.FIXED * SD_slope_raster_15_data_pts * SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_LCA, SD_add.gam_LCA.smoothed, SD_add.gam_LCA.smoothed_first_term, 
    SD_add.gam_LCA.smoothed_second_term, SD_add.gam_LCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
SD_add.gam_LCA.smoothed.quasi <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA, family = quasi())
SD_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA, family = poisson())
SD_add.gam_LCA.smoothed.quasipoisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                    data = SD_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(SD_add.gam_LCA.smoothed, SD_add.gam_LCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_LCA.smoothed.quasi, SD_add.gam_LCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_LCA.smoothed.quasi, SD_add.gam_LCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(SD_add.gam_LCA.smoothed.poisson, SD_add.gam_LCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_LCA)
summary(SD_add.gam_LCA.smoothed)
summary(SD_add.gam_LCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(SD_add.gam_LCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_add.gam_LCA.smoothed.dredged <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts), 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)

#Chosen model: SD_add.gam_LCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
SD_add.gam_LCA.smoothed.poisson <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)
k.check(SD_add.gam_LCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(SD_add.gam_LCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(SD_add.gam_LCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(SD_add.gam_LCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(SD_add.gam_LCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(SD_add.gam_LCA.smoothed, select = "s(SD_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = SD_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Long Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(SD_add.gam_LCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(SD_add.gam_LCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_SD_add.gam_LCA <- fitted.values(SD_add.gam_LCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(SD_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=SD_slope_raster_15_data_pts, 
                                                                z=Canopy_long, color=SD_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = SD_fixed_field_data_processed_terrain_no_NA$Canopy_long ~ 
                SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)



# CA
SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts
SD_add.gam_CA <- gam(Canopy_area ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CA.smoothed_first_term <- gam(Canopy_area ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CA.smoothed_second_term <- gam(Canopy_area ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CA_interact <- gam(Canopy_area ~ Elevation..m.FIXED * SD_slope_raster_15_data_pts * SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_CA, SD_add.gam_CA.smoothed, SD_add.gam_CA.smoothed_first_term, 
    SD_add.gam_CA.smoothed_second_term, SD_add.gam_CA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
SD_add.gam_CA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA, family = quasi())
SD_add.gam_CA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA, family = poisson())
SD_add.gam_CA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                    data = SD_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(SD_add.gam_CA.smoothed, SD_add.gam_CA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_CA.smoothed.quasi, SD_add.gam_CA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_CA.smoothed.quasi, SD_add.gam_CA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(SD_add.gam_CA.smoothed.poisson, SD_add.gam_CA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_CA)
summary(SD_add.gam_CA.smoothed)
summary(SD_add.gam_CA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(SD_add.gam_CA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_add.gam_CA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts), 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)

#Chosen model: SD_add.gam_CA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
SD_add.gam_CA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)
k.check(SD_add.gam_CA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(SD_add.gam_CA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(SD_add.gam_CA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(SD_add.gam_CA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(SD_add.gam_CA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(SD_add.gam_CA.smoothed, select = "s(SD_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Canopy Area", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = SD_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Canopy Area", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(SD_add.gam_CA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(SD_add.gam_CA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_SD_add.gam_CA <- fitted.values(SD_add.gam_CA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(SD_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=SD_slope_raster_15_data_pts, 
                                                                z=Canopy_area, color=SD_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = SD_fixed_field_data_processed_terrain_no_NA$Canopy_area ~ 
                SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)



# CS


SD_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CS_interact <- gam(Crown_spread ~ Elevation..m.FIXED * SD_slope_raster_15_data_pts * SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_CS, SD_add.gam_CS.smoothed, SD_add.gam_CS.smoothed_first_term, 
    SD_add.gam_CS.smoothed_second_term, SD_add.gam_CS_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
SD_add.gam_CS.smoothed.quasi <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA, family = quasi())
SD_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA, family = poisson())
SD_add.gam_CS.smoothed.quasipoisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                    data = SD_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(SD_add.gam_CS.smoothed, SD_add.gam_CS.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_CS.smoothed.quasi, SD_add.gam_CS.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_CS.smoothed.quasi, SD_add.gam_CS.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(SD_add.gam_CS.smoothed.poisson, SD_add.gam_CS.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_CS)
summary(SD_add.gam_CS.smoothed)
summary(SD_add.gam_CS.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(SD_add.gam_CS.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_add.gam_CS.smoothed.dredged <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts), 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)

#Chosen model: SD_add.gam_CS.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
SD_add.gam_CS.smoothed.poisson <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)
k.check(SD_add.gam_CS.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(SD_add.gam_CS.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(SD_add.gam_CS.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(SD_add.gam_CS.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(SD_add.gam_CS.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(SD_add.gam_CS.smoothed, select = "s(SD_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = SD_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Canopy Spread", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(SD_add.gam_CS.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Canopy Spread")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(SD_add.gam_CS.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_SD_add.gam_CS <- fitted.values(SD_add.gam_CS.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(SD_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=SD_slope_raster_15_data_pts, 
                                                                z=Crown_spread, color=SD_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = SD_fixed_field_data_processed_terrain_no_NA$Crown_spread ~ 
                SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


# DBH_ag


SD_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_DBH_interact <- gam(DBH_ag ~ Elevation..m.FIXED * SD_slope_raster_15_data_pts * SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_DBH, SD_add.gam_DBH.smoothed, SD_add.gam_DBH.smoothed_first_term, 
    SD_add.gam_DBH.smoothed_second_term, SD_add.gam_DBH_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
SD_add.gam_DBH.smoothed.quasi <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA, family = quasi())
SD_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA, family = poisson())
SD_add.gam_DBH.smoothed.quasipoisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                    data = SD_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(SD_add.gam_DBH.smoothed, SD_add.gam_DBH.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_DBH.smoothed.quasi, SD_add.gam_DBH.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(SD_add.gam_DBH.smoothed.quasi, SD_add.gam_DBH.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(SD_add.gam_DBH.smoothed.poisson, SD_add.gam_DBH.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_DBH)
summary(SD_add.gam_DBH.smoothed)
summary(SD_add.gam_DBH.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(SD_add.gam_DBH.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
SD_add.gam_DBH.smoothed.dredged <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts), 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)

#Chosen model: SD_add.gam_DBH.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
SD_add.gam_DBH.smoothed.poisson <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                               data = SD_fixed_field_data_processed_terrain_no_NA)
k.check(SD_add.gam_DBH.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(SD_add.gam_DBH.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(SD_add.gam_DBH.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(SD_add.gam_DBH.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(SD_add.gam_DBH.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(SD_add.gam_DBH.smoothed, select = "s(SD_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  geom_smooth(se = T) + 
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on DBH", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = SD_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on DBH", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(SD_add.gam_DBH.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(SD_add.gam_DBH.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_SD_add.gam_DBH <- fitted.values(SD_add.gam_DBH.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(SD_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=SD_slope_raster_15_data_pts, 
                                                                z=DBH_ag, color=SD_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = SD_fixed_field_data_processed_terrain_no_NA$DBH_ag ~ 
                SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts + 
                SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

