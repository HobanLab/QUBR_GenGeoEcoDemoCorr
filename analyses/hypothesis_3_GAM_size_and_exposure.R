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

# plot(all_points_multiple_lm_SCA)
# all_points_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
# all_points_mlm_SCA_cooks <- cooks.distance(all_points_mlm_SCA) #calculating the cook.s D for each point
# plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
# influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (2 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
# influential

#SCA


all_points_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
all_points_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_SCA, all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed_first_term, 
    all_points_add.gam_SCA.smoothed_second_term)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_SCA)
summary(all_points_add.gam_SCA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_SCA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(all_points_add.gam_SCA.smoothed.dredge, all_points_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(all_points_add.gam_SCA.smoothed.dredge, all_points_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_SCA.smoothed

#checking K to see if we 
k.check(all_points_add.gam_SCA.smoothed.dredge)
k.check(all_points_add.gam_SCA.smoothed)

#no interaction plots
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
 plot.gam(all_points_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

#looking for interaction
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_SCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_SCA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


# LCA
all_points_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
all_points_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
all_points_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_LCA, all_points_add.gam_LCA.smoothed, all_points_add.gam_LCA.smoothed_first_term, 
    all_points_add.gam_LCA.smoothed_second_term)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_LCA)
summary(all_points_add.gam_LCA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_LCA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 
#the full model is the dredge output

#Chosen model: all_points_add.gam_LCA.smoothed

#checking K to see if we 
k.check(all_points_add.gam_LCA.smoothed)

#plotting the gam results (run in one chunk)
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_LCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
all_points_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_LCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_SCA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot


# CA


all_points_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
all_points_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
all_points_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_CA, all_points_add.gam_CA.smoothed, all_points_add.gam_CA.smoothed_first_term, 
    all_points_add.gam_CA.smoothed_second_term)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_CA)
summary(all_points_add.gam_CA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_CA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(all_points_add.gam_CA.smoothed.dredge, all_points_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(all_points_add.gam_CA.smoothed.dredge, all_points_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_CA.smoothed.dredge
summary(all_points_add.gam_CA.smoothed.dredge)

#checking K to see if we 
k.check(all_points_add.gam_CA.smoothed.dredge)
k.check(all_points_add.gam_CA.smoothed)


par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_CA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
all_points_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_CA.smoothed.inter <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_CA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_CA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# CS


all_points_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                             data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
all_points_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                 data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
all_points_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_CS, all_points_add.gam_CS.smoothed, all_points_add.gam_CS.smoothed_first_term, 
    all_points_add.gam_CS.smoothed_second_term)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_CS)
summary(all_points_add.gam_CS.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_CS.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                              data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(all_points_add.gam_CS.smoothed.dredge, all_points_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(all_points_add.gam_CS.smoothed.dredge, all_points_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_CA.smoothed.dredge
summary(all_points_add.gam_CS.smoothed.dredge)

#checking K to see if we 
k.check(all_points_add.gam_CS.smoothed.dredge)
k.check(all_points_add.gam_CS.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_CA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
all_points_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                            data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_CS.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_CS.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot


# DBH_ag

all_points_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                             data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
all_points_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                 data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
all_points_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_DBH, all_points_add.gam_DBH.smoothed, all_points_add.gam_DBH.smoothed_first_term, 
    all_points_add.gam_DBH.smoothed_second_term)
anova(all_points_add.gam_DBH, all_points_add.gam_DBH.smoothed_first_term, 
    all_points_add.gam_DBH.smoothed_second_term, all_points_add.gam_DBH.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_DBH)
summary(all_points_add.gam_DBH.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_DBH.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_DBH.smoothed.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                              data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(all_points_add.gam_DBH.smoothed.dredge, all_points_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(all_points_add.gam_DBH.smoothed.dredge, all_points_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_DBH.smoothd.dredge
summary(all_points_add.gam_DBH.smoothed.dredge)

#checking K to see if we 
k.check(all_points_add.gam_DBH.smoothed.dredge)
k.check(all_points_add.gam_DBH.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_DBH.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
all_points_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
all_points_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                            data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_DBH.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_DBH.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot



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
#LM_fixed_field_data_processed_terrain_no_NA_No_outliers <- LM_fixed_field_data_processed_terrain_no_NA[-c(24,26,27),]


# SCA

LM_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
LM_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
LM_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_SCA, LM_add.gam_SCA.smoothed_first_term, 
    LM_add.gam_SCA.smoothed_second_term, LM_add.gam_SCA.smoothed)
anova(LM_add.gam_SCA, LM_add.gam_SCA.smoothed_first_term, 
      LM_add.gam_SCA.smoothed_second_term, LM_add.gam_SCA.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_SCA)
summary(LM_add.gam_SCA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_SCA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED),
                                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_SCA.smoothed.dredge

summary(LM_add.gam_SCA.smoothed.dredge)

#checking K to see if we 
k.check(LM_add.gam_SCA.smoothed.dredge)

par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_SCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LM_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_SCA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot



# LCA



LM_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
LM_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
LM_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_LCA, LM_add.gam_LCA.smoothed, LM_add.gam_LCA.smoothed_first_term, 
    LM_add.gam_LCA.smoothed_second_term)
anova(LM_add.gam_LCA, LM_add.gam_LCA.smoothed_second_term, 
      LM_add.gam_LCA.smoothed_first_term, LM_add.gam_LCA.smoothed)


#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_LCA)
summary(LM_add.gam_LCA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_LCA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 
#the full model is the dredge output

#fitting the dredged model
LM_add.gam_SCA.smoothed.dredge <-  gam(Canopy_long ~ s(Elevation..m.FIXED),
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

anova(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_LCA.smoothed)

#Chosen model: LM_add.gam_SCA.smoothed.dredge

#checking K to see if we need to change the K value
k.check(LM_add.gam_SCA.smoothed.dredge)

#plotting the gam results (run in one chunk)
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_LCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
LM_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_LCA.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_LCA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot



# CA

LM_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                             data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
LM_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                      data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
LM_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_CA, LM_add.gam_CA.smoothed, LM_add.gam_CA.smoothed_first_term, 
    LM_add.gam_CA.smoothed_second_term)
anova(LM_add.gam_CA, LM_add.gam_CA.smoothed_first_term, 
    LM_add.gam_CA.smoothed_second_term, LM_add.gam_CA.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_CA)
summary(LM_add.gam_CA.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_CA.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED), 
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(LM_add.gam_CA.smoothed.dredge, LM_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(LM_add.gam_CA.smoothed.dredge, LM_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_CA.smoothed.dredge
summary(LM_add.gam_CA.smoothed.dredge)

#checking K to see if might need to change K (if it is significantly low)
k.check(LM_add.gam_CA.smoothed.dredge)

par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_CA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LM_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_CA.smoothed.inter <- gam(Canopy_area ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_CA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_CA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# CS


LM_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                             data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
LM_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                      data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
LM_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_CS, LM_add.gam_CS.smoothed, LM_add.gam_CS.smoothed_first_term, 
    LM_add.gam_CS.smoothed_second_term)
anova(LM_add.gam_CS, LM_add.gam_CS.smoothed_first_term, 
    LM_add.gam_CS.smoothed_second_term, LM_add.gam_CS.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_CS)
summary(LM_add.gam_CS.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_CS.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED),
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(LM_add.gam_CS.smoothed.dredge, LM_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(LM_add.gam_CS.smoothed.dredge, LM_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_CA.smoothed.dredge
summary(LM_add.gam_CS.smoothed.dredge)

#checking K to see if we 
k.check(LM_add.gam_CS.smoothed.dredge)
k.check(LM_add.gam_CS.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_CA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LM_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                    data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_CS.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_CS.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot


# DBH_ag
LM_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
LM_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
LM_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LM_add.gam_DBH, LM_add.gam_DBH.smoothed, LM_add.gam_DBH.smoothed_first_term, 
    LM_add.gam_DBH.smoothed_second_term)
anova(LM_add.gam_DBH, LM_add.gam_DBH.smoothed_first_term, 
    LM_add.gam_DBH.smoothed_second_term, LM_add.gam_DBH.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LM_add.gam_DBH)
summary(LM_add.gam_DBH.smoothed)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(LM_add.gam_DBH.smoothed) #using the dredge model to narro the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_DBH.smoothd.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Anova F test comparing strength of dredge vs. full model demonstrates dredge performs just as well.
anova(LM_add.gam_DBH.smoothd.dredge, LM_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model 
AIC(LM_add.gam_DBH.smoothd.dredge, LM_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_DBH.smoothed
summary(LM_add.gam_DBH.smoothed)

#checking K to see if we 
k.check(LM_add.gam_DBH.smoothed.dredge)
k.check(LM_add.gam_DBH.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_DBH.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LM_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
LM_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                    data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_DBH.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot



# LC

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
LC_fixed_field_data_processed_terrain_no_NA <- LC_fixed_field_data_processed_terrain %>%
  filter(is.na(LC_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)

# SCA

options(na.action = "na.omit")
LC_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                      data = LC_fixed_field_data_processed_terrain) #na fail makes sure the later dredge does not have to worry about NAs
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain)
LC_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain)
LC_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                           data = LC_fixed_field_data_processed_terrain)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_SCA, LC_add.gam_SCA.smoothed_first_term, 
    LC_add.gam_SCA.smoothed_second_term, LC_add.gam_SCA.smoothed)
anova(LC_add.gam_SCA, LC_add.gam_SCA.smoothed_first_term, 
      LC_add.gam_SCA.smoothed_second_term, LC_add.gam_SCA.smoothed)

#Because dredge was not working with LC, I compared with and without slope and aspect 
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_SCA.smoothed.slope.less <- gam(Canopy_short ~ s(Elevation..m.FIXED)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_SCA.smoothed.elevation.less <- gam(Canopy_short ~ s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                              data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(LC_add.gam_SCA.smoothed.slope.less, LC_add.gam_SCA.smoothed)
anova(LC_add.gam_SCA.smoothed.elevation.less, LC_add.gam_SCA.smoothed)
#We do not really need slope (since the one with only elevation was not significant, so the one with both variables 
#was not as essential as just with elevation)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_SCA)
summary(LC_add.gam_SCA.smoothed)

#Chosen model: LC_add.gam_SCA.smoothed


#checking K to see if we 
k.check(LC_add.gam_SCA.smoothed)

par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_SCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_SCA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_SCA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot


# LCA
options(na.action = "na.omit")
LC_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA) #na fail makes sure the later dredge does not have to worry about NAs
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_LCA, LC_add.gam_LCA.smoothed, LC_add.gam_LCA.smoothed_first_term, 
    LC_add.gam_LCA.smoothed_second_term)
anova(LC_add.gam_LCA, LC_add.gam_LCA.smoothed_first_term, 
    LC_add.gam_LCA.smoothed_second_term, LC_add.gam_LCA.smoothed)


#Because dredge was not working with LC, I compared with and without slope and aspect 
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_LCA.smoothed.slope.less <- gam(Canopy_long ~ s(Elevation..m.FIXED)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_LCA.smoothed.elevation.less <- gam(Canopy_long ~ s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                              data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(LC_add.gam_LCA.smoothed.slope.less, LC_add.gam_LCA.smoothed)
anova(LC_add.gam_LCA.smoothed.elevation.less, LC_add.gam_LCA.smoothed)
#We do not really need slope (since the one with only elevation was not significant, so the one with both variables 
#was not as essential as just with elevation)


#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_LCA)
summary(LC_add.gam_LCA.smoothed)

#Chosen model: LC_add.gam_LCA.smoothed

#checking K to see if we 
k.check(LC_add.gam_LCA.smoothed)

#plotting the gam results (run in one chunk)
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_LCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_LCA.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_LCA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot



# CA

LC_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                             data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
LC_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
LC_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
LC_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_CA, LC_add.gam_CA.smoothed, LC_add.gam_CA.smoothed_first_term, 
    LC_add.gam_CA.smoothed_second_term)
anova(LC_add.gam_CA, LC_add.gam_CA.smoothed_first_term, 
    LC_add.gam_CA.smoothed_second_term, LC_add.gam_CA.smoothed)

summary(LC_add.gam_CA.smoothed)
summary(LC_add.gam_CA.smoothed_first_term)
#we can see that slope does not appear to be as useful

#Because dredge was not working with LC, I compared with and without slope and aspect 
LC_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_CA.smoothed.slope.less <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_CA.smoothed.elevation.less <- gam(log(Canopy_area) ~ s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                              data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(LC_add.gam_CA.smoothed.slope.less, LC_add.gam_CA.smoothed)
anova(LC_add.gam_CA.smoothed.elevation.less, LC_add.gam_CA.smoothed)
#We do not really need slope (since the one with only elevation was not significant, so the one with both variables 
#was not as essential as just with elevation)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_CA)
summary(LC_add.gam_CA.smoothed)
summary(LC_add.gam_CA.smoothed.slope.less)

#Chosen model: LC_add.gam_CA.smoothed


#checking K to see if we 
k.check(LC_add.gam_CA.smoothed)


par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_CA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
LC_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA.smoothed.inter <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                    data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_CA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_CA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# CS


LC_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                             data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
LC_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)

LC_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_CS, LC_add.gam_CS.smoothed, LC_add.gam_CS.smoothed_first_term, 
    LC_add.gam_CS.smoothed_second_term)
anova(LC_add.gam_CS, LC_add.gam_CS.smoothed_first_term, 
    LC_add.gam_CS.smoothed_second_term, LC_add.gam_CS.smoothed)

summary(LC_add.gam_CS.smoothed)
summary(LC_add.gam_CS.smoothed_first_term)
#we can see that elevation and slope does not appear to be as useful

#Because dredge was not working with LC, I compared with and without slope and aspect 
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_CS.smoothed.slope.less <- gam(Crown_spread ~ s(Elevation..m.FIXED)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                         data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_CS.smoothed.elevation.less <- gam(Crown_spread ~ s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(LC_add.gam_CS.smoothed.slope.less, LC_add.gam_CS.smoothed)
anova(LC_add.gam_CS.smoothed.elevation.less, LC_add.gam_CS.smoothed)
#We do not really need slope (since the one with only elevation was not significant, so the one with both variables 
#was not as essential as just with elevation)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_CS)
summary(LC_add.gam_CS.smoothed)

#Chosen model: LC_add.gam_CS.smoothed

#checking K to see if we 
k.check(LC_add.gam_CS.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_CA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                    data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_CS.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_CS.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot



# DBH_ag

LC_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA) #na fail makes sure the later dredge does not have to worry about NAs
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(LC_add.gam_DBH, LC_add.gam_DBH.smoothed, LC_add.gam_DBH.smoothed_first_term, 
    LC_add.gam_DBH.smoothed_second_term)
anova(LC_add.gam_DBH, LC_add.gam_DBH.smoothed_first_term, 
    LC_add.gam_DBH.smoothed_second_term, LC_add.gam_DBH.smoothed)

summary(LC_add.gam_DBH.smoothed)
summary(LC_add.gam_DBH.smoothed_first_term)
#we can see that elevation and slope does not appear to be as useful

#Because dredge was not working with LC, I compared with and without slope and aspect 
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_DBH.smoothed.slope.less <- gam(DBH_ag ~ s(Elevation..m.FIXED)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                         data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
LC_add.gam_DBH.smoothed.elevation.less <- gam(DBH_ag ~ s(LC_slope_raster_15_data_pts)+ LC_aspect_raster_15_data_pts_8_categorical, 
                                             data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(LC_add.gam_DBH.smoothed.slope.less, LC_add.gam_DBH.smoothed)
anova(LC_add.gam_DBH.smoothed.elevation.less, LC_add.gam_DBH.smoothed)
#We do not really need slope or elevation by themselves (since the one with only elevation and one with only slope was not significant)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(LC_add.gam_DBH)
summary(LC_add.gam_DBH.smoothed)

#Chosen model: LC_add.gam_DBH.smoothed
summary(LC_add.gam_DBH.smoothed)

#checking K to see if we 
k.check(LC_add.gam_DBH.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_DBH.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_DBH.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot



# SD


#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
SD_fixed_field_data_processed_terrain_no_NA <- SD_fixed_field_data_processed_terrain %>%
  filter(is.na(SD_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)


# SCA

options(na.action = "na.omit")
SD_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                      data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                          data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                           data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_SCA, SD_add.gam_SCA.smoothed_first_term, 
    SD_add.gam_SCA.smoothed_second_term, SD_add.gam_SCA.smoothed)
anova(SD_add.gam_SCA, SD_add.gam_SCA.smoothed_first_term, 
      SD_add.gam_SCA.smoothed_second_term, SD_add.gam_SCA.smoothed)

#Because dredge was not working with SD, I compared with and without slope and aspect 
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_SCA.smoothed.slope.less <- gam(Canopy_short ~ s(Elevation..m.FIXED)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                          data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_SCA.smoothed.elevation.less <- gam(Canopy_short ~ s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                              data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(SD_add.gam_SCA.smoothed.slope.less, SD_add.gam_SCA.smoothed)
anova(SD_add.gam_SCA.smoothed.elevation.less, SD_add.gam_SCA.smoothed)
#We do not really need elevation (since the one with only slope was not significant, so the one with both variables 
#was not as essential as just with slope)


#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_SCA)
summary(SD_add.gam_SCA.smoothed)

#Chosen model: SD_add.gam_SCA.smoothed

#checking K to see if we 
k.check(SD_add.gam_SCA.smoothed)

par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_SCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_SCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_SCA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot


# LCA

SD_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)


#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_LCA, SD_add.gam_LCA.smoothed, SD_add.gam_LCA.smoothed_first_term, 
    SD_add.gam_LCA.smoothed_second_term)
anova(SD_add.gam_LCA, SD_add.gam_LCA.smoothed_first_term, 
    SD_add.gam_LCA.smoothed_second_term, SD_add.gam_LCA.smoothed)

#Because dredge was not working with SD, I compared with and without slope and aspect 
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_LCA.smoothed.slope.less <- gam(Canopy_long ~ s(Elevation..m.FIXED)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                          data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_LCA.smoothed.elevation.less <- gam(Canopy_long ~ s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                              data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(SD_add.gam_LCA.smoothed.slope.less, SD_add.gam_LCA.smoothed)
anova(SD_add.gam_LCA.smoothed.elevation.less, SD_add.gam_LCA.smoothed)
#We do not really need elevation (since the one with only slope was not significant, so the one with both variables 
#was not as essential as just with slope)


#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_LCA)
summary(SD_add.gam_LCA.smoothed)

#Chosen model: SD_add.gam_LCA.smoothed

#checking K to see if we 
k.check(SD_add.gam_LCA.smoothed)

#plotting the gam results (run in one chunk)
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_LCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_LCA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_LCA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot



# CA


SD_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                             data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
SD_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                 data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_CA, SD_add.gam_CA.smoothed, SD_add.gam_CA.smoothed_first_term, 
    SD_add.gam_CA.smoothed_second_term)

summary(SD_add.gam_CA.smoothed)
#we can see that elevation does not appear to be as useful

#Because dredge was not working with LC, I compared with and without slope and aspect 
SD_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CA.smoothed.slope.less <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                         data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CA.smoothed.elevation.less <- gam(log(Canopy_area) ~ s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(SD_add.gam_CA.smoothed.slope.less, SD_add.gam_CA.smoothed)
anova(SD_add.gam_CA.smoothed.elevation.less, SD_add.gam_CA.smoothed)
#We do not really need elevation (since the one with only slope was not significant, so the one with both variables 
#was not as essential as just with slope)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_CA)
summary(SD_add.gam_CA.smoothed)

#Chosen model: SD_add.gam_CA.smoothed


#checking K to see if we 
k.check(SD_add.gam_CA.smoothed)


par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_CA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
SD_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CA.smoothed.inter <- gam(Canopy_area ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                    data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_CA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_CA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# CS


SD_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                             data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                 data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_CS, SD_add.gam_CS.smoothed, SD_add.gam_CS.smoothed_first_term, 
    SD_add.gam_CS.smoothed_second_term)
AIC(SD_add.gam_CS, SD_add.gam_CS.smoothed_first_term, 
    SD_add.gam_CS.smoothed_second_term, SD_add.gam_CS.smoothed)

summary(SD_add.gam_CS.smoothed)
#we can see that elevation does not appear to be as useful

#Because dredge was not working with LC, I compared with and without slope and aspect 
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CS.smoothed.slope.less <- gam(Crown_spread ~ s(Elevation..m.FIXED)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                         data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_CS.smoothed.elevation.less <- gam(Crown_spread ~ s(SD_slope_raster_15_data_pts)+ SD_aspect_raster_15_data_pts_8_categorical, 
                                             data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(SD_add.gam_CS.smoothed.slope.less, SD_add.gam_CS.smoothed)
anova(SD_add.gam_CS.smoothed.elevation.less, SD_add.gam_CS.smoothed)
#We do not really need elevation (since the one with only slope was not significant, so the one with both variables 
#was not as essential as just with slope)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.smoothed)

#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_CS)
summary(SD_add.gam_CS.smoothed)

#Chosen model: SD_add.gam_CS.smoothed


#checking K to see if we 
k.check(SD_add.gam_CS.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_CA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


#looking for interaction
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                    data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_CS.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_CS.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot



# DBH_ag

SD_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit) #na fail makes sure the later dredge does not have to worry about NAs
SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
SD_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.omit)
#logging canopy area lower the AIC significantly

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(SD_add.gam_DBH, SD_add.gam_DBH.smoothed, SD_add.gam_DBH.smoothed_first_term, 
    SD_add.gam_DBH.smoothed_second_term)
anova(SD_add.gam_DBH, SD_add.gam_DBH.smoothed_first_term, 
    SD_add.gam_DBH.smoothed_second_term, SD_add.gam_DBH.smoothed)


SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_DBH.smoothed.slope.less <- gam(DBH_ag ~ s(Elevation..m.FIXED) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
SD_add.gam_DBH.smoothed.elevation.less <- gam(DBH_ag ~ s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                        data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.omit)
anova(SD_add.gam_DBH.smoothed.slope.less, SD_add.gam_DBH.smoothed)
anova(SD_add.gam_DBH.smoothed.elevation.less, SD_add.gam_DBH.smoothed)
#We do not need elevation,(since the one with only slope was not significant, so the one with both variables 
#was not as essential as just with slope)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.smoothed)
gam.check(SD_add.gam_DBH.smoothed.elevation.less)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest
summary(SD_add.gam_DBH)
summary(SD_add.gam_DBH.smoothed)
summary(SD_add.gam_DBH.smoothed.slope.less)

#Chosen model: SD_add.gam_DBH.smoothed.slope.less
summary(SD_add.gam_DBH.smoothed)


#plotting the model
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_DBH.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

#looking for interaction
SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_DBH.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot


