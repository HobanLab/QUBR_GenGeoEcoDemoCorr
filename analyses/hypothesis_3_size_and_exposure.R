#### Loading libraries and relevant data ####

library(googledrive) #to download files from google drive
library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(perm.t.test) #permutation t test 
library(car) #to create added variable plots and to run levene's test for checking ANOVA conditions
library(stars) # to convert raster into stars
library(gdalUtilities) #to be able to use gdalwarp



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

LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")

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
  mutate(Elevation..m.FIXED = case_when((Elevation..m. == 360) ~ NA,  #change the elevation of 360 which appears to be a miswritten elevation to NA
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

#HERE IS THE IMPORTATION OF THE CROPPED RASTERS
r <- rast() 
crs(r)
#Importing the cropped rasters for LM, LC, and SD and setting crs
CEM_15_utm_LM <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LM.tif"))
terra::crs(CEM_15_utm_LM) <- CRS("+init=epsg:26912")

CEM_15_utm_LC <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LC.tif"))
terra::crs(CEM_15_utm_LC) <- CRS("+init=epsg:26912")

CEM_15_utm_SD <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_SD.tif"))
terra::crs(CEM_15_utm_SD) <- CRS("+init=epsg:26912")

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

#plot the aspects
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

#re-categorizing the aspect data

#setting values of 360 to 0 

#all points
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts = case_when((all_points_aspect_raster_15_data_pts == "360") ~  0,
                                                  (all_points_aspect_raster_15_data_pts != "360")~ all_points_aspect_raster_15_data_pts))
View(all_points_fixed_field_data_processed_terrain)
  
#LM
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts = case_when((LM_aspect_raster_15_data_pts == 360) ~  0,
                                                  (LM_aspect_raster_15_data_pts != 360)~ LM_aspect_raster_15_data_pts))
View(LM_fixed_field_data_processed_terrain)
summary(LM_fixed_field_data_processed_terrain)

#LC

LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts = case_when((LC_aspect_raster_15_data_pts == 360) ~  0, 
                                                  (LC_aspect_raster_15_data_pts != 360) ~ LC_aspect_raster_15_data_pts))

View(LC_fixed_field_data_processed_terrain)
summary(LC_fixed_field_data_processed_terrain)

tmp <- LC_fixed_field_data_processed_terrain %>%
  filter(QUBR_ID %in% c("LC_153", "LC_073", "LC_217")) %>%
  dplyr::select(LC_aspect_raster_15_data_pts)

#SD

SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts = case_when((SD_aspect_raster_15_data_pts == 360) ~  0,
                                                  (SD_aspect_raster_15_data_pts != 360)~ SD_aspect_raster_15_data_pts))
View(SD_fixed_field_data_processed_terrain)

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

View(all_points_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_4_categorical = case_when((all_points_aspect_raster_15_data_pts >= 0 & all_points_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

View(all_points_fixed_field_data_processed_terrain)




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
  
View(LM_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
LM_fixed_field_data_processed_terrain <- LM_fixed_field_data_processed_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_4_categorical = case_when((LM_aspect_raster_15_data_pts >= 0 & LM_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                              (LM_aspect_raster_15_data_pts >= 315 & LM_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                              (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                              (LM_aspect_raster_15_data_pts >= 135 & LM_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                              (LM_aspect_raster_15_data_pts >= 225 & LM_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315
                                                     
View(LM_fixed_field_data_processed_terrain)

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

View(LC_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_4_categorical = case_when((LC_aspect_raster_15_data_pts >= 0 & LC_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 315 & LC_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                (LC_aspect_raster_15_data_pts >= 45 & LC_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LC_aspect_raster_15_data_pts >= 135 & LC_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LC_aspect_raster_15_data_pts >= 225 & LC_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

View(LC_fixed_field_data_processed_terrain)


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

View(SD_fixed_field_data_processed_terrain)

# North, East, South, West

# the directions are a range of 90 degrees 
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_4_categorical = case_when((SD_aspect_raster_15_data_pts >= 0 & SD_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 315 & SD_aspect_raster_15_data_pts < 359.99999) ~ "N",
                                                                (SD_aspect_raster_15_data_pts >= 45 & SD_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (SD_aspect_raster_15_data_pts >= 135 & SD_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (SD_aspect_raster_15_data_pts >= 225 & SD_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

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


### Sizes vs. Elevation ###

# For all trees

fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  drop_na(Elevation..m.FIXED)

#SCA

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_sca_elev, aes(x= all_points_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_sca_elev, aes(sample = all_points_lm_sca_elev$residuals))+
  geom_qq()

shapiro.test(all_points_lm_sca_elev$residuals) #only not signficant for a square root transformation, we could use mann-kendall test for non-parametric data or the square root transformation

#checking equal variance
ggplot(data = all_points_lm_sca_elev, aes(x = all_points_lm_sca_elev$fitted.values, y = all_points_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_sca_elev)

#correlation test
cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_short)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_SCA <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_SCA)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Canopy_short ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#Cook's D
all_points_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates)
all_points_slr_LCA_cooks <- cooks.distance(all_points_slr_LCA) #calculating the cook.s D for each point
plot(all_points_slr_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- all_points_slr_LCA_cooks[(all_points_slr_LCA_cooks > (3 * mean(all_points_slr_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 49, 86, 87, 102, 123, 164, 172, 184, 
                                                                                                                      185, 189, 203, 205, 207, 211, 213, 233, 235, 242, 243, 245, 273, 274,
                                                                                                                      278, 287, 288, 289, 291, 300, 320, 325, 357, 358, 359, 360, 361, 363, 366,
                                                                                                                      380, 415, 429, 453, 467, 487, 544, 620, 632),]
#creating the linear regression
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#without outliers

all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long_lg ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_LM_lm_lca_elev, aes(x= all_points_LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_lca_elev, aes(sample = all_points_lm_lca_elev$residuals))+
  geom_qq()

shapiro.test(all_points_lm_lca_elev$residuals) #shapiro-wilk test is sig despite transformations and removal of outliers, so we have to use Mann-Kendall non-parametric test

#checking equal variance
ggplot(data = all_points_lm_lca_elev, aes(x = all_points_lm_lca_elev$fitted.values, y = all_points_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_lca_elev)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_LCA <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_LCA)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#canopy area

#checking linearity 

#plotting the linear model in ggplot for CA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#Cook's D
all_points_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates)
all_points_slr_CA_cooks <- cooks.distance(all_points_slr_CA) #calculating the cook.s D for each point
plot(all_points_slr_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- all_points_slr_CA_cooks[(all_points_slr_CA_cooks > (3 * mean(all_points_slr_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 87, 153, 164, 172, 185, 194, 203, 205, 207, 211, 213, 237, 241, 242, 243, 287, 
                                                                                                                     288, 291, 300, 303, 320, 325, 357, 358, 359, 360, 361, 366,
                                                                                                                     380, 411, 415, 429, 453, 467, 487, 544, 620, 632),]
                                                                                                                      
#creating the linear regression
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#without outliers

#creating the linear regression
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area_lg ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CA_elev, aes(x= all_points_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CA_elev, aes(sample = all_points_lm_CA_elev$residuals))+
  geom_qq()

shapiro.test(all_points_lm_CA_elev$residuals) #shapiro-wilk test is sig depsite transformations and removal of outliers, so we need to use a non-parametric mann-kendall test

#checking equal variance
ggplot(data = all_points_lm_CA_elev, aes(x = all_points_lm_CA_elev$fitted.values, y = all_points_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_CA_elev)


#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_CA <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_CA)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#Cook's D
all_points_slr_CS <- lm(Crown_spread ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates)
all_points_slr_CS_cooks <- cooks.distance(all_points_slr_CS) #calculating the cook.s D for each point
plot(all_points_slr_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- all_points_slr_CS_cooks[(all_points_slr_CS_cooks > (3 * mean(all_points_slr_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 45, 86, 87, 164, 172, 184, 185, 189, 203, 205, 207, 211, 213, 235, 237, 241, 242, 243, 
                                                                                                                     273, 274, 278, 287, 288, 291, 300, 303, 320, 325, 357, 358, 359, 360, 361, 363, 366,
                                                                                                                     380, 415, 429, 453, 467, 487, 544, 620, 632),]

#creating the linear regression
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#logged CS
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#square rooted CS
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#without outliers

all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)

#logged CS
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread_lg ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)

#square root CS
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CS_elev, aes(x= all_points_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CS_elev, aes(sample = all_points_lm_CS_elev$residuals))+
  geom_qq()

shapiro.test(all_points_lm_CS_elev$residuals) #all were sig, despite transformations and removal of outliers, so need to use non-parametric mann-kendall test

#checking equal variance
ggplot(data = all_points_lm_CS_elev, aes(x = all_points_lm_CS_elev$fitted.values, y = all_points_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_CS_elev)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_CS <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_CS)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#checking linearity 

#plotting the linear model in ggplot for DBH
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#removing outliers

#Cook's D
all_points_slr_DBH <- lm(DBH_ag ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates)
all_points_slr_DBH_cooks <- cooks.distance(all_points_slr_DBH) #calculating the cook.s D for each point
plot(all_points_slr_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- all_points_slr_DBH_cooks[(all_points_slr_DBH_cooks > (3 * mean(all_points_slr_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(16, 49, 64, 77, 87, 93, 122, 149, 157, 164, 172, 201, 203, 205, 207, 213, 238, 239, 241, 242, 243, 
                                                                                                                     245, 266, 273, 287, 288, 290, 296, 300, 303, 320, 325, 335, 351, 355, 356, 358, 360, 361, 366,
                                                                                                                     380, 429, 585, 620, 632),]


#creating the linear regression
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)


#creating linear regressions without outliers
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag_lg ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_DBH_elev, aes(x= all_points_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_DBH_elev, aes(sample = all_points_lm_DBH_elev$residuals))+
  geom_qq()

shapiro.test(all_points_lm_DBH_elev$residuals) #only non-sign result is without outliers and with square rooted dbh

#plotting the linear model in ggplot for DBH
ggplot(data = fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers, (aes(x=Elevation..m.FIXED, y=DBH_ag_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(DBH)")


#checking equal variance
ggplot(data = all_points_lm_DBH_elev, aes(x = all_points_lm_DBH_elev$fitted.values, y = all_points_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_DBH_elev)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_DBH <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_DBH)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

# LM 

#we had to remove the elevation NA to 
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  drop_na(Elevation..m.FIXED)
length(LM_fixed_field_data_processed$Elevation..m.FIXED)

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis (m)")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12))

#Cook's D
LM_slr_SCA <- lm(Canopy_short ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed)
LM_slr_SCA_cooks <- cooks.distance(LM_slr_SCA) #calculating the cook.s D for each point
plot(LM_slr_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_SCA_cooks[(LM_slr_SCA_cooks > (3 * mean(LM_slr_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_sca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 44, 110, 122, 123, 144, 164, 172, 184, 189, 203, 205, 213),]

#creating the linear regression
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed_sca_no_outliers$Canopy_short ~ LM_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#logged sca
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#squared sca
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_elev, aes(x= LM_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_elev, aes(sample = LM_lm_sca_elev$residuals))+
  geom_qq()

shapiro.test(LM_lm_sca_elev$residuals) #significantly not normal, except when outliers are removed

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

#Cook's D
LM_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed)
LM_slr_LCA_cooks <- cooks.distance(LM_slr_LCA) #calculating the cook.s D for each point
plot(LM_slr_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_LCA_cooks[(LM_slr_LCA_cooks > (3 * mean(LM_slr_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_lca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 42, 49, 110, 123, 138, 144, 164, 172, 184, 185, 189, 203, 205, 207, 211, 212, 213),]


#creating the linear regression
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#logged lca
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#square root lca
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#without outliers
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#logged lca
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long_lg ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#square root lca
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long_sqrt ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_lca_elev, aes(x= LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_lca_elev, aes(sample = LM_lm_lca_elev$residuals))+
  geom_qq()

shapiro.test(LM_lm_lca_elev$residuals) #sign not normally distributed with just transformations,but non sign with square root trans and removal of outliers

#checking equal variance
ggplot(data = LM_lm_lca_elev, aes(x = LM_lm_lca_elev$fitted.values, y = LM_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_lca_elev)

#non parametric Mann-Kendall Test
LM_tau_result_LCA <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_LCA)

# Calculate the trend line
LM_trend_line_LCA <- predict(loess(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#canopy area

#checking linearity 

#plotting the linear model in ggplot for CA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#Cook's D
LM_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed)
LM_slr_CA_cooks <- cooks.distance(LM_slr_CA) #calculating the cook.s D for each point
plot(LM_slr_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_CA_cooks[(LM_slr_CA_cooks > (3 * mean(LM_slr_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_ca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 44, 87, 110, 122, 123, 144, 164, 172, 184, 203, 205, 207, 211, 213),]

#creating the linear regression

LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#without outliers

LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area_lg ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area_sqrt ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(Canopy Area_")

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CA_elev, aes(x= LM_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CA_elev, aes(sample = LM_lm_CA_elev$residuals))+
  geom_qq()

shapiro.test(LM_lm_CA_elev$residuals) #with just transformations sig, but with removal of outliers and square root trans it is normal

#checking equal variance
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CA_elev)


#non parametric Mann-Kendall Test
LM_tau_result_CA <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CA)

# Calculate the trend line
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

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

#linear regression with logged transformation of aggregated CS
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated CS
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CS_elev, aes(x= LM_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CS_elev, aes(sample = LM_lm_CS_elev$residuals))+
  geom_qq()

shapiro.test(LM_lm_CS_elev$residuals) #not sign with square root transformation

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(Crown Spread)")

#checking equal variance
ggplot(data = LM_lm_CS_elev, aes(x = LM_lm_CS_elev$fitted.values, y = LM_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CS_elev)


#non parametric Mann-Kendall Test
LM_tau_result_CS <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CS)

# Calculate the trend line
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed$Crown_spread ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#checking linearity 


#creating the linear regression
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("log(DBH)")+
  abline(LM_lm_DBH_elev)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12))

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_DBH_elev, aes(x= LM_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_DBH_elev, aes(sample = LM_lm_DBH_elev$residuals))+
  geom_qq()

shapiro.test(LM_lm_DBH_elev$residuals) #not sign with logged transformation of DBH

#checking equal variance
ggplot(data = LM_lm_DBH_elev, aes(x = LM_lm_DBH_elev$fitted.values, y = LM_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_DBH_elev)

ggplot()+
  geom_sf(data = LM_fixed_field_data_processed, aes(color = Elevation..m.FIXED))


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

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("log(Short Canopy Axis)")


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_sca_elev, aes(x= LC_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_elev, aes(sample = LC_lm_sca_elev$residuals))+
  geom_qq()

shapiro.test(LC_lm_sca_elev$residuals) #not sign with logged transformation of sca

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

shapiro.test(LC_lm_lca_elev$residuals) #not sign with logged transformation

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

shapiro.test(LC_lm_CA_elev$residuals) #not sign with log transformation of CA

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

#plotting the linear model in ggplot for CS
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

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("log(Crown Spread)")

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CS_elev, aes(x= LC_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CS_elev, aes(sample = LC_lm_CS_elev$residuals))+
  geom_qq()

shapiro.test(LC_lm_CS_elev$residuals) #not sign with logged transformation

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

#plotting the linear model in ggplot for DBH
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(DBH)")

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_DBH_elev, aes(x= LC_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_DBH_elev, aes(sample = LC_lm_DBH_elev$residuals))+
  geom_qq()

shapiro.test(LC_lm_DBH_elev$residuals) #not sign with square root transformation

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

SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  drop_na(Elevation..m.FIXED)

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
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("log(Short Canopy Axis)")

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_sca_elev, aes(x= SD_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_elev, aes(sample = SD_lm_sca_elev$residuals))+
  geom_qq()

shapiro.test(SD_lm_sca_elev$residuals) #not sign when logged

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
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_lca_elev, aes(x= SD_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_lca_elev, aes(sample = SD_lm_lca_elev$residuals))+
  geom_qq()

shapiro.test(SD_lm_lca_elev$residuals) # not sign when logged

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

shapiro.test(SD_lm_CA_elev$residuals) #not sign when CA is logged

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

SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CS_elev, aes(x= SD_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CS_elev, aes(sample = SD_lm_CS_elev$residuals))+
  geom_qq()

shapiro.test(SD_lm_CS_elev$residuals) #not sign when logged CS

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

shapiro.test(SD_lm_DBH_elev$residuals) #not sign when logged

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

# For all trees

#SCA

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_sca_elev, aes(x= all_points_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_sca_elev, aes(sample = all_points_lm_sca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_lm_sca_elev, aes(x = all_points_lm_sca_elev$fitted.values, y = all_points_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_sca_elev)

#correlation test
cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression

all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_LM_lm_lca_elev, aes(x= all_points_LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_lca_elev, aes(sample = all_points_lm_lca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_lm_lca_elev, aes(x = all_points_lm_lca_elev$fitted.values, y = all_points_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CA_elev, aes(x= all_points_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CA_elev, aes(sample = all_points_lm_CA_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_lm_CA_elev, aes(x = all_points_lm_CA_elev$fitted.values, y = all_points_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_CA_elev)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CS_elev, aes(x= all_points_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CS_elev, aes(sample = all_points_lm_CS_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_lm_CS_elev, aes(x = all_points_lm_CS_elev$fitted.values, y = all_points_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_CS_elev)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_DBH_elev, aes(x= all_points_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_DBH_elev, aes(sample = all_points_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = all_points_lm_DBH_elev, aes(x = all_points_lm_DBH_elev$fitted.values, y = all_points_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(all_points_lm_DBH_elev)

### linear models comparing slope to size/shape ###

#all points 

#removing NAs preventing us from running tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_slope_raster_15_data_pts)





#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x= all_points_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#creating the linear regression
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#log transformation
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#square root transformation
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_sca_slope, aes(x= all_points_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_sca_slope, aes(sample = all_points_lm_sca_slope$residuals))+
  geom_qq()

shapiro.test(all_points_lm_sca_slope$residuals) #shaprio-welk, never not sig, use the mann-kendall test

#checking equal variance
ggplot(data = all_points_lm_sca_slope, aes(x = all_points_lm_sca_slope$fitted.values, y = all_points_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(all_points_lm_sca_slope)

#correlation test
cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_short)


#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_SCA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_SCA)

# Calculate the trend line
LC_trend_line_LCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_long ~ LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression

all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_lca_slope, aes(x= all_points_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_lca_slope, aes(sample = all_points_lm_lca_slope$residuals))+
  geom_qq()

shapiro.test(all_points_lm_lca_slope$residuals) # shapiro-welk test, all were sig, meaning we will use the Mann-Kendall Test

#checking equal variance
ggplot(data = all_points_lm_lca_slope, aes(x = all_points_lm_lca_slope$fitted.values, y = all_points_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(all_points_lm_lca_slope)


#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_LCA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_LCA)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CA_slope, aes(x= all_points_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CA_slope, aes(sample = all_points_lm_CA_slope$residuals))+
  geom_qq()

shapiro.test(all_points_lm_CA_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test

#checking equal variance
ggplot(data = all_points_lm_CA_slope, aes(x = all_points_lm_CA_slope$fitted.values, y = all_points_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(all_points_lm_CA_slope)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_CA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_CA)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression

all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_CS_slope, aes(x= all_points_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_CS_slope, aes(sample = all_points_lm_CS_slope$residuals))+
  geom_qq()

shapiro.test(all_points_lm_CS_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test


#checking equal variance
ggplot(data = all_points_lm_CS_slope, aes(x = all_points_lm_CS_slope$fitted.values, y = all_points_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(all_points_lm_CS_slope)


#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_CS <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_CS)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#creating the linear regression
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with logged transformation of aggregated DBH
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of aggregated DBH
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(all_points_lm_DBH_slope, aes(x= all_points_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(all_points_lm_DBH_slope, aes(sample = all_points_lm_DBH_slope$residuals))+
  geom_qq()

shapiro.test(all_points_lm_DBH_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test

#checking equal variance
ggplot(data = all_points_lm_DBH_slope, aes(x = all_points_lm_DBH_slope$fitted.values, y = all_points_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(all_points_lm_DBH_slope)

#non parametric Mann-Kendall Test for the version without outliers
all_points_tau_result_DBH <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(all_points_tau_result_DBH)

# Calculate the trend line
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


# LM 



#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x= LM_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#Cook's D
LM_lm_focal_SCA <- lm(DBH_ag ~ sum_SCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_SCA_cooks <- cooks.distance(LM_lm_focal_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_all_focal_trees_no_DBH_outliers <- LM_fixed_field_data_all_focal_trees[-c(3,23,24),]


#creating the linear regression
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_slope, aes(x= LM_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_slope, aes(sample = LM_lm_sca_slope$residuals))+
  geom_qq()

shapiro.test(LM_lm_sca_slope$residuals) #shapiro wilk test 

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
cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, LM_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_lca_slope, aes(x= LM_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_lca_slope, aes(sample = LM_lm_lca_slope$residuals))+
  geom_qq()

shapiro.test(LM_lm_lca_slope$residuals)

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
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CA_slope, aes(x= LM_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CA_slope, aes(sample = LM_lm_CA_slope$residuals))+
  geom_qq()

shapiro.test(LM_lm_CA_slope$residuals)

#checking equal variance
ggplot(data = LM_lm_CA_slope, aes(x = LM_lm_CA_slope$fitted.values, y = LM_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(LM_lm_CA_slope)

#non parametric Mann-Kendall Test
LM_tau_result_CA <- cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, LM_fixed_field_data_processed_terrain$Canopy_area_lg,  method = "kendall")

# Print Kendall's tau and its associated p-value 
print(LM_tau_result_CA)

# Calculate the trend line
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed_terrain$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()



#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression

LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CS_slope, aes(x= LM_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CS_slope, aes(sample = LM_lm_CS_slope$residuals))+
  geom_qq()

shapiro.test(LM_lm_CS_slope$residuals)

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
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#creating the linear regression
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with logged transformation of aggregated DBH
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of aggregated DBH
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_DBH_slope, aes(x= LM_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_DBH_slope, aes(sample = LM_lm_DBH_slope$residuals))+
  geom_qq()

shapiro.test(LM_lm_DBH_slope$residuals)

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
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")


#creating the linear regression
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with logged transformation of short canopy axis
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_sca_slope, aes(x= LC_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_slope, aes(sample = LC_lm_sca_slope$residuals))+
  geom_qq()

shapiro.test(LC_lm_sca_slope$residuals)

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
cor.test(LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts, LC_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Canopy_long_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear transformation with logged long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_lca_slope, aes(x= LC_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_lca_slope, aes(sample = LC_lm_lca_slope$residuals))+
  geom_qq()

shapiro.test(LC_lm_lca_slope$residuals)

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
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CA_slope, aes(x= LC_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CA_slope, aes(sample = LC_lm_CA_slope$residuals))+
  geom_qq()

shapiro.test(LC_lm_CA_slope$residuals)

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
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#creating the linear regression
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear transformation with logged crown spread
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear transformation with square rooted crown spread
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CS_slope, aes(x= LC_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CS_slope, aes(sample = LC_lm_CS_slope$residuals))+
  geom_qq()

shapiro.test(LC_lm_CS_slope$residuals)

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
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope")+
  ylab("DBH")

#creating the linear regression
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with logged transformation of aggregated DBH
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of aggregated DBH
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_DBH_slope, aes(x= LC_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_DBH_slope, aes(sample = LC_lm_DBH_slope$residuals))+
  geom_qq()

shapiro.test(LC_lm_DBH_slope$residuals)

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
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#creating the linear regression
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_sca_slope, aes(x= SD_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_slope, aes(sample = SD_lm_sca_slope$residuals))+
  geom_qq()

shapiro.test(SD_lm_sca_slope$residuals)

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
cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, SD_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#creating the linear regression
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_lca_slope, aes(x= SD_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_lca_slope, aes(sample = SD_lm_lca_slope$residuals))+
  geom_qq()

shapiro.test(SD_lm_lca_slope$residuals)

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
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#creating the linear regression
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of canopy area
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of canopy area
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CA_slope, aes(x= SD_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CA_slope, aes(sample = SD_lm_CA_slope$residuals))+
  geom_qq()

shapiro.test(SD_lm_CA_slope$residuals)

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
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Crown_spread_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("sqrt(Crown Spread)")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

#creating the linear regression
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with logged transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CS_slope, aes(x= SD_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CS_slope, aes(sample = SD_lm_CS_slope$residuals))+
  geom_qq()

shapiro.test(SD_lm_CS_slope$residuals)

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
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#creating the linear regression
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with logged transformation of aggregated DBH
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of aggregated DBH
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_DBH_slope, aes(x= SD_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_DBH_slope, aes(sample = SD_lm_DBH_slope$residuals))+
  geom_qq()

shapiro.test(SD_lm_DBH_slope$residuals)

#checking equal variance
ggplot(data = SD_lm_DBH_slope, aes(x = SD_lm_DBH_slope$fitted.values, y = SD_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test visible in summary of the lm
summary(SD_lm_DBH_slope)

#non parametric Mann-Kendall Test
SD_tau_result_DBH <- cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, SD_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall") #DBH_ag_sqrt

# Print Kendall's tau and its associated p-value
print(SD_tau_result_DBH)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


## Size vs. Aspect ##

# we ran ANOVAs to test difference in size means between cardinal directions

#8 categories for direction


#all points 

#removing NAs preventing us from running tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_aspect_raster_15_data_pts_8_categorical)



#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
all_points_aov_SCA_aspect_8 <- aov(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_SCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_SCA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(all_points_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
all_points_thumb_test_SCA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd)
max(all_points_thumb_test_SCA, na.rm = T) / min(all_points_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
all_points_aov_LCA_aspect_8 <- aov(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_LCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_LCA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_long, 
                                          all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(all_points_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test, sign have to run a non-paramrtric test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
all_points_thumb_test_LCA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd)
max(all_points_thumb_test_LCA, na.rm = T) / min(all_points_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

# canopy area

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
all_points_aov_CA_aspect_8 <- aov(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_CA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_area, 
                                         all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(all_points_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)


#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
all_points_thumb_test_CA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd)
max(all_points_thumb_test_CA, na.rm = T) / min(all_points_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#crown spread

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
all_points_aov_CS_aspect_8 <- aov(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CS_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_CS_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Crown_spread, 
                                         all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(all_points_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
all_points_thumb_test_CS <- tapply(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd)
max(all_points_thumb_test_CS, na.rm = T) / min(all_points_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
all_points_aov_DBH_aspect_8 <- aov(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_DBH_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_DBH_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$DBH_ag, 
                                          all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(all_points_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
all_points_thumb_test_DBH <- tapply(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd)
max(all_points_thumb_test_DBH, na.rm = T) / min(all_points_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# LM

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LM_aov_SCA_aspect_8 <- aov(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect_8)

# checking to see if residuals are normal
hist(LM_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(LM_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LM_thumb_test_SCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd)
max(LM_thumb_test_SCA, na.rm = T) / min(LM_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
LM_aov_LCA_aspect_8 <- aov(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_LCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_LCA_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_long, 
                                          LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(LM_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LM_thumb_test_LCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd)
max(LM_thumb_test_LCA, na.rm = T) / min(LM_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

# canopy area

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
LM_aov_CA_aspect_8 <- aov(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_CA_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_area, 
                                         LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(LM_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)


#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LM_thumb_test_CA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd)
max(LM_thumb_test_CA, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#crown spread

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
LM_aov_CS_aspect_8 <- aov(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CS_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_CS_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Crown_spread, 
                                         LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(LM_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LM_thumb_test_CS <- tapply(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd)
max(LM_thumb_test_CS, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
LM_aov_DBH_aspect_8 <- aov(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_DBH_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_DBH_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$DBH_ag, 
                                          LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(LM_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LM_thumb_test_DBH <- tapply(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd)
max(LM_thumb_test_DBH, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# LC

#short canopy axis
 
#had to remove the NAs to be able to run the function
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(LC_aspect_raster_15_data_pts_8_categorical)


#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LC_aov_SCA_aspect_8 <- aov(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
esummary(LC_aov_SCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(LC_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd)
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
LC_aov_LCA_aspect_8 <- aov(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_LCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_LCA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_long, 
                                          LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(LC_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)


# checking equal variances with levene's test and rule of thumb

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LC_thumb_test_LCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd)
max(LC_thumb_test_LCA, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# canopy area

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
LC_aov_CA_aspect_8 <- aov(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_CA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_area, 
                                         LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_t_test_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(LC_t_test_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)


#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LC_thumb_test_CA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd)
max(LC_thumb_test_CA, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#crown spread

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
LC_aov_CS_aspect_8 <- aov(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CS_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_CS_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Crown_spread, 
                                         LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_t_test_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(LC_t_test_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LC_thumb_test_CS <- tapply(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd)
max(LC_thumb_test_CS, na.rm = T) / min(LC_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
LC_aov_DBH_aspect_8 <- aov(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_DBH_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_DBH_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$DBH_ag, 
                                          LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_t_test_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(LC_t_test_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd)
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# SD

#short canopy axis

#had to remove the NAs to be able to run the function
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(SD_aspect_raster_15_data_pts_8_categorical)


#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
SD_aov_SCA_aspect_8 <- aov(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_SCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_aov_SCA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(SD_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
SD_thumb_test_SCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd)
max(SD_thumb_test_SCA, na.rm = T) / min(SD_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method




#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
SD_aov_LCA_aspect_8 <- aov(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_LCA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_LCA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_long, 
                                          SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(SD_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
SD_thumb_test_LCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd)
max(SD_thumb_test_LCA, na.rm = T) / min(SD_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# canopy area

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
SD_aov_CA_aspect_8 <- aov(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CA_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_CA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_area, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(SD_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
SD_thumb_test_CA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd)
max(SD_thumb_test_CA, na.rm = T) / min(SD_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#crown spread

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
SD_aov_CS_aspect_8 <- aov(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CS_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_CS_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Crown_spread, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(SD_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
SD_thumb_test_CS <- tapply(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd)
max(SD_thumb_test_CS, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")+
  ylim(c(0,1.5))+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

#ANOVA
SD_aov_DBH_aspect_8 <- aov(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_DBH_aspect_8)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_DBH_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$DBH_ag, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(SD_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#rule of thumb test
SD_thumb_test_DBH <- tapply(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd)
max(SD_thumb_test_DBH, na.rm = T) / min(SD_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#4 categories for direction

#all points 

#removing NAs preventing us from running tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_aspect_raster_15_data_pts_4_categorical)

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
all_points_aov_SCA_aspect_4 <- aov(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_SCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_SCA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categoricalv, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(all_points_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
all_points_thumb_test_SCA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd)
max(all_points_thumb_test_SCA_4, na.rm = T) / min(all_points_thumb_test_SCA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
all_points_aov_LCA_aspect_4 <- aov(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_LCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_LCA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_long, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(all_points_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
all_points_thumb_test_LCA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd)
max(all_points_thumb_test_LCA_4, na.rm = T) / min(all_points_thumb_test_LCA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

# canopy area

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
all_points_aov_CA_aspect_4 <- aov(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_CA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_area, 
                                                 all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(all_points_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)


#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
all_points_thumb_test_CA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd)
max(all_points_thumb_test_CA_4, na.rm = T) / min(all_points_thumb_test_CA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#crown spread

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
all_points_aov_CS_aspect_4 <- aov(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CS_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_CS_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Crown_spread, 
                                                 all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(all_points_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
all_points_thumb_test_CS_4 <- tapply(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd)
max(all_points_thumb_test_CS_4, na.rm = T) / min(all_points_thumb_test_CS_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
all_points_aov_DBH_aspect_4 <- aov(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_DBH_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
all_points_t_test_DBH_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$DBH_ag, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(all_points_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#levene's test
leveneTest(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
all_points_thumb_test_DBH_4 <- tapply(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd)
max(all_points_thumb_test_DBH_4, na.rm = T) / min(all_points_thumb_test_DBH_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# LM

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LM_aov_SCA_aspect_4 <- aov(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_SCA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_short, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(LM_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LM_thumb_test_SCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd)
max(LM_thumb_test_SCA, na.rm = T) / min(LM_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_SCA_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_short, 
                                          LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
LM_aov_LCA_aspect_4 <- aov(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LCMaov_LCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_LCA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_long, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_t_test_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(LM_t_test_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LM_thumb_test_LCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd)
max(LM_thumb_test_LCA, na.rm = T) / min(LM_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

# canopy area

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
LM_aov_CA_aspect_4 <- aov(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_CA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_area, 
                                       LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(LM_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LM_thumb_test_CA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd)
max(LM_thumb_test_CA, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

#crown spread

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
LM_aov_CS_aspect_4 <- aov(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CS_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_CS_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Crown_spread, 
                                       LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(LM_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LM_thumb_test_CS <- tapply(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd)
max(LM_thumb_test_CS, na.rm = T) / min(LM_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
LM_aov_DBH_aspect_4 <- aov(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_DBH_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LM_t_test_DBH_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$DBH_ag, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(LM_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LM_fixed_field_data_processed_terrain$LM_t_test_DBH_aspect_4 ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LM_thumb_test_DBH <- tapply(LM_fixed_field_data_processed_terrain$LM_t_test_DBH_aspect_4, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd)
max(LM_thumb_test_DBH, na.rm = T) / min(LM_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method

# LC

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
LC_aov_SCA_aspect_4 <- aov(Canopy_short ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_SCA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(LC_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when normal data
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LC_thumb_test_SCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd)
max(LC_thumb_test_SCA, na.rm = T) / min(LC_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
LC_aov_LCA_aspect_4 <- aov(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_LCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_LCA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_long, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(LC_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LC_thumb_test_LCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd)
max(LC_thumb_test_LCA, na.rm = T) / min(LC_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# canopy area

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
LC_aov_CA_aspect_4 <- aov(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_CA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_area, 
                                       LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(LC_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LC_thumb_test_CA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd)
max(LC_thumb_test_CA, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#crown spread

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
LC_aov_CS_aspect_4 <- aov(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CS_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_CS_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Crown_spread, 
                                       LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(LC_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LC_thumb_test_CS <- tapply(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd)
max(LC_thumb_test_CS, na.rm = T) / min(LC_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
LC_aov_DBH_aspect_4 <- aov(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_DBH_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
LC_t_test_DBH_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$DBH_ag, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(LC_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test, sig so we need to use a non-parametric test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#levene's test
leveneTest(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd)
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# SD

#short canopy axis

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#ANOVA
SD_aov_SCA_aspect_4 <- aov(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_SCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_SCA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect")

qqnorm(SD_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test, sig so need to use a non-parametric test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
SD_thumb_test_SCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd)
max(SD_thumb_test_SCA, na.rm = T) / min(SD_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#long canopy axis

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#ANOVA
SD_aov_LCA_aspect_4 <- aov(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_LCA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_LCA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_long, 
                                          SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(SD_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test, sig so need to use a non-parametric test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
SD_thumb_test_LCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd)
max(SD_thumb_test_LCA, na.rm = T) / min(SD_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


# canopy area

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#ANOVA
SD_aov_CA_aspect_4 <- aov(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CA_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_CA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_area, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect")

qqnorm(SD_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
SD_thumb_test_CA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd)
max(SD_thumb_test_CA, na.rm = T) / min(SD_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#crown spread

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#ANOVA
SD_aov_CS_aspect_4 <- aov(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CS_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_CS_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Crown_spread, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect")

qqnorm(SD_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
SD_thumb_test_CS <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd)
max(SD_thumb_test_CS, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#DBH ag

#boxplot of sizes over the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#ANOVA
SD_aov_DBH_aspect_4 <- aov(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_DBH_aspect_4)

#permutation t.test to see significant differences between categories using a bonferonni adjustment
SD_t_test_DBH_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$DBH_ag, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect")

qqnorm(SD_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test, sig so have to use a non-parametric test

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when 
fligner.test(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#levene's test
leveneTest(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#rule of thumb test
SD_thumb_test_DBH <- tapply(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd)
max(SD_thumb_test_DBH, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#nonparametric tests

#kruskall wallis test
kruskal.test(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method




#Boxplots to look at LC NW, SW, and N tended to have larger means

#grouping NW, N, SW, and the other directions
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  add_column(LC_aspect_raster_15_data_pts_8_regrouped = NA) %>% #creating a new column to regroup the directions
  mutate(LC_aspect_raster_15_data_pts_8_regrouped = case_when(LC_aspect_raster_15_data_pts_8_categorical == "NE" ~ "Other", 
                                                              LC_aspect_raster_15_data_pts_8_categorical == "E" ~  "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "NW" ~ "NW",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "N" ~ "N",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "W" ~ "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "SW" ~ "SW",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "S" ~ "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "SE" ~ "Other",
                                                              ))

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_short))+
  xlab("Direction")+
  ylab("Short Canopy Axis")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_long))+
  xlab("Direction")+
  ylab("Long Canopy Axis")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_area))+
  xlab("Direction")+
  ylab("Canopy Area")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Crown_spread))+
  xlab("Direction")+
  ylab("Crown Spread")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = DBH_ag))+
  xlab("Direction")+
  ylab("DBH")


#Boxplot of SD to look at NE, N, and W tended to have larger means

SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  add_column(SD_aspect_raster_15_data_pts_8_regrouped = NA) %>% #creating a new column to regroup the directions
  mutate(SD_aspect_raster_15_data_pts_8_regrouped = case_when(SD_aspect_raster_15_data_pts_8_categorical == "NE" ~ "NE", 
                                                              SD_aspect_raster_15_data_pts_8_categorical == "E" ~  "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "NW" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "N" ~ "N",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "W" ~ "W",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "SW" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "S" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "SE" ~ "Other",
  ))

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_short))+
  xlab("Direction")+
  ylab("Short Canopy Axis")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_long))+
  xlab("Direction")+
  ylab("Long Canopy Axis")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_area))+
  xlab("Direction")+
  ylab("Canopy Area")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Crown_spread))+
  xlab("Direction")+
  ylab("Crown Spread")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = DBH_ag))+
  xlab("Direction")+
  ylab("DBH")


# Plotting size of shape/size values for points in population with coloring based on aspect

# LC

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_short, color = LC_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_long, color = LC_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_area, color = LC_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Crown_spread, color = LC_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = DBH_ag, color = LC_aspect_raster_15_data_pts_8_regrouped))


#SD

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_short, color = SD_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_long, color = SD_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_area, color = SD_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Crown_spread, color = SD_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = DBH_ag, color = SD_aspect_raster_15_data_pts_8_regrouped))

