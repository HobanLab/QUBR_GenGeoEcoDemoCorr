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
library(mgcv) #to use GAM function 
library(plotly) #to 3d plot variables
library(MuMIn) #to be able to use dredge
library(visreg) # to be able to plot Aspect/categorical variables with GAM

#### Loading and processing relevant data ####

# loading in the tree data (size, elevation, lat/lon, ID, size/shape)

fixed_field_data_processed_source <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#adding a sequential column, "X," to number each tree

fixed_field_data_processed_source <- fixed_field_data_processed_source %>%
  mutate(X = row_number())

# creating the point shapefiles of the tree locations for each population in UTM 12 N

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
fixed_field_data_processed_source_sf <- st_as_sf(fixed_field_data_processed_source, 
                                          coords = c("long", "lat"), crs = 4326)

#creating a transformed point shapefile with UTM 12 N an equal area projection
fixed_field_data_processed_source_sf_transformed_source <- st_transform(fixed_field_data_processed_source_sf, crs = 26912)

#storing point shapefiles for the trees by population

LM_fixed_field_data_processed_source_source_sf_source <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "LM") %>%
  st_as_sf()

LC_fixed_field_data_processed_source_source_sf_source <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "LC") %>%
  st_as_sf()

SD_fixed_field_data_processed_source_source_sf_source <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "SD") %>%
  st_as_sf()

#create dataframe with X and Y UTM coordinates

fixed_field_data_processed_source_sf_trans_coords_source <- st_coordinates(fixed_field_data_processed_source_sf_transformed_source) #creates a dataframe with separate x and y columns from the UTM 12N transformation
fixed_field_data_processed_source_sf_trans_coordinates_source <- fixed_field_data_processed_source_sf_transformed_source %>%
  cbind(fixed_field_data_processed_source_sf_trans_coords_source) #combines the x and y coordinate data frame with the transformed sf dataframe

#### transformations of tree size/shape variables (log and square root) for linear models ####

#creating columns with transformations: logged all of the size/shape variables
fixed_field_data_processed_source_sf_trans_coordinates_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the size/shape variables
fixed_field_data_processed_source_sf_trans_coordinates_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#setting elevation as a numeric value
fixed_field_data_processed_source_sf_trans_coordinates_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

# Creating fixed_field_data_processed_source dataframes for each population with the nearest neighbor columns

LM_fixed_field_data_processed_source_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed_source_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed_source_source <- fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  filter(Locality == "SD")

#### Fixing errors with the elevation data ####

#creating a new column in the whole dataset to get rid of 360 m outlier and turn the values in feet into meter
fixed_field_data_processed_source_sf_trans_coordinates_source <-  fixed_field_data_processed_source_sf_trans_coordinates_source %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. < 700 & Elevation..m. != 360) ~ Elevation..m.,
                                        (Elevation..m. == 360) ~ NA, 
                                        (Elevation..m. > 700) ~ Elevation..m.*0.3048))  #because LM and LC do not have a 360 elevation and SD and LC do have values above 700, this should not effect them

#creating a new elevation column so the values that were mistakenly put in feet are in meters
LM_fixed_field_data_processed_source_source <-  LM_fixed_field_data_processed_source_source %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.)) %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))

#creating a new elevation column so LC, LM, and SD all have this same column, makes it easier for combining the population data frames
LC_fixed_field_data_processed_source_source <-  LC_fixed_field_data_processed_source_source %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))

#plotting all of the tree points by elevation (m) to check the range of values 
ggplot()+
  geom_sf(data = fixed_field_data_processed_source_sf_trans_coordinates_source, aes(color = Elevation..m.FIXED))

#creating a new elevation column so that a 360 m outlier is 460
SD_fixed_field_data_processed_source_source <-  SD_fixed_field_data_processed_source_source %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. == 360) ~ NA, 
                                        (Elevation..m. != 360) ~ Elevation..m.))

#### Loading in ArcGIS river and Baja California Sur shapefile and storing out polygons for each population ####

#Las Matancitas (LM)
river_LM_source <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LM River/LM_Rivers_Final.shp")
river_LM_source  <- river_LM_source$geometry[1]
plot(river_LM_source)

#La Cobriza (LC)
river_LC_source  <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LC River/LC_Rivers_Final.shp")
river_LC_source  <- river_LC_source$geometry[1]
plot(river_LC_source)

#San Dionisio (SD)
river_SD_source <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/SD River/SD_Rivers_Final.shp")
river_SD_source <- river_SD_source$geometry[1]
plot(river_SD_source)

#changing the coordinate reference system of the river polygons to be equal area projection (UTM 12N), uses meters as distance measurement 
#uses meters as distance measurement
river_LM_trans_source <- st_transform(river_LM_source, crs = 26912) 
river_LC_trans_source <- st_transform(river_LC_source, crs = 26912)
river_SD_trans_source <- st_transform(river_SD_source, crs = 26912)

#ensuring the river outlines are shapefiles for the distance measurements
river_LM_trans_source <- st_sf(geometry = river_LM_trans_source)
river_LC_trans_source <- st_sf(geometry = river_LC_trans_source)
river_SD_trans_source <- st_sf(geometry = river_SD_trans_source)

# Importing BCS and River Shapefiles 

#turning the Baja California Sur polygon into a shapefile, to be able to visualize the point locations
BCS_polygon_source <- read_sf("./data/Shapefiles/BCS_Shapefile/bcs_entidad.shp")
BCS_polygon_source <- st_as_sf(BCS_polygon_source) #ensures foreign will be an sf object (in this case a multipolygon)

#### Creating buffers and boundaries around the rivers ####

#creating buffers around the rivers

#LM
river_buffer_LM_source <- st_buffer(river_LM_trans_source, 100) #100 m buffer
ggplot()+  #plotting the river shapefile, the buffer, and the tree points
  geom_sf(data = river_buffer_LM_source)+
  geom_sf(data = river_LM_trans_source)+
  geom_sf(data = LM_fixed_field_data_processed_source_source_sf_source)

#LC
river_buffer_LC_source <- st_buffer(river_LC_trans_source, 100) #230 m buffer
ggplot()+ #plotting the river shapefile, the buffer, and the tree points
  geom_sf(data = river_buffer_LC_source)+
  geom_sf(data = river_LC_trans_source)+
  geom_sf(data = LC_fixed_field_data_processed_source_source_sf_source)

#SD
river_buffer_SD_source <- st_buffer(river_SD_trans_source, 70) #70 m buffer
ggplot()+ #plotting the river shapefile, the buffer, and the tree points
  geom_sf(data = river_buffer_SD_source)+
  geom_sf(data = river_SD_trans_source)+
  geom_sf(data = SD_fixed_field_data_processed_source_source_sf_source)

#creating bounding boxes for each population

#creating a boundary box for LM with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
LM_fixed_field_data_processed_source_source_box <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "LM") %>%
  st_bbox %>%
  st_as_sfc()

#creating a boundary box for LC with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
LC_fixed_field_data_processed_source_source_box <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "LC") %>%
  st_bbox %>%
  st_as_sfc()

#creating a boundary box for SD with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
SD_fixed_field_data_processed_source_source_box <- fixed_field_data_processed_source_sf_transformed_source %>%
  filter(Locality == "SD") %>%
  st_bbox %>%
  st_as_sfc()

#creating bounding boxes for all of the river shapefiles for each population
LM_box <- st_bbox(river_LM_trans_source)
LC_box <- st_bbox(river_LC_trans_source)
SD_box <- st_bbox(river_SD_trans_source)

#### Creating the elevation, aspect, and slope columns in the dataframe ####

# elevation data is from INEGI 15 m and allows us to calculate slope and aspect


#BECAUSE THE ELEVATION RASTERS WERE TOO BIG TO DOWNLOAD DIRECTLY, FROM GOOGLE DRIVE, OR OPEN FROM A ZIP, 
#AFTER LOADING IN THE ORIGINAL DATA, WE CROPPED IT TO FIT OUR POPULATIONS, EXPORTED THOSE FILES, AND THEN DOWNLOADED THOSE
#HERE IS A LINK TO A GOOGLE DRIVE WITH THE INGEI 15 m CONTINUOUS ELEVATION MODEL DATA: 
# https://drive.google.com/drive/folders/17RxjebifsRFFS4iEucDQtMFaqjzRI-Ss?usp=sharing 

#SO WE COMMENTED OUT THE CODE WE USED TO LOAD IN THE ORIGINAL RASTER AND CREATE AND EXPORT THE CROPPED RASTERS FOR EACH POPULATION THAT WE DOWNLOAD LATER

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
#   geom_sf(data = fixed_field_data_processed_source_sf_trans_coordinates_source)
# 
# #LM
# 
# #mapping cropped 
# CEM_15_utm_LM <- crop(CEM_15_utm, extent((c(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LM, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LM_fixed_field_data_processed_source_source)
# 
# #LC
# 
# #mapping cropped 
# CEM_15_utm_LC <- crop(CEM_15_utm, extent((c(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LC, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LC_fixed_field_data_processed_source_source)
# 
# #SD
# 
# #mapping cropped 
# CEM_15_utm_SD <- crop(CEM_15_utm, extent((c(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_SD, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = SD_fixed_field_data_processed_source_source)

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
#   geom_sf(data = fixed_field_data_processed_source_sf_trans_coordinates_source)
# 
# #LM
# 
# #mapping cropped 
# CEM_15_utm_LM <- crop(CEM_15_utm, extent((c(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LM, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LM_fixed_field_data_processed_source_source)
# 
# #LC
# 
# #mapping cropped 
# CEM_15_utm_LC <- crop(CEM_15_utm, extent((c(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_LC, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = LC_fixed_field_data_processed_source_source)
# 
# #SD
# 
# #mapping cropped 
# CEM_15_utm_SD <- crop(CEM_15_utm, extent((c(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))))
# 
# #plotting the LM elevation raster with the LM points
# ggplot()+
#   geom_raster(data= as.data.frame(CEM_15_utm_SD, xy = T), aes(x=x, y=y, fill = CEM_15_utm))+
#   geom_sf(data = SD_fixed_field_data_processed_source_source)


#Importing the cropped rasters for LM, LC, and SD and setting the crs to the same as the points
CEM_15_utm_LM <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LM.tif"))
terra::crs(CEM_15_utm_LM) <- CRS("+init=epsg:26912")

CEM_15_utm_LC <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LC.tif"))
terra::crs(CEM_15_utm_LC) <- CRS("+init=epsg:26912")

CEM_15_utm_SD <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_SD.tif"))
terra::crs(CEM_15_utm_SD) <- CRS("+init=epsg:26912")

#creating the all points raster by merging the LM, LC, and SD rasters
CEM_15_utm_all_points <- raster::merge(CEM_15_utm_LM, CEM_15_utm_LC, CEM_15_utm_SD)

#plotting all of the elevation rasters and tree points
ggplot()+
  geom_raster(data= as.data.frame(CEM_15_utm_all_points, xy = T), aes(x=x, y=y, fill = layer))+
  geom_sf(data = fixed_field_data_processed_source_sf_transformed_source)

## Extracting the slope 

#all points 

#extracting the slope in degrees, using the queens method (neighbor = 8)
all_points_slope_raster_15 <- terra::terrain(CEM_15_utm_all_points, unit = 'degrees', neighbors = 8, 'slope')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(all_points_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = fixed_field_data_processed_source_sf_trans_coordinates_source)+
  scale_fill_viridis_c()


#LM

#extracting the slope in degrees, using the queens method (neighbor = 8)
LM_slope_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'slope')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LM_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()

#LC

#extracting the slope in degrees, using the queens method (neighbor = 8)
LC_slope_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'slope')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = LC_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()

#SD

#extracting the slope in degrees, using the queens method (neighbor = 8)
SD_slope_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'slope')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_slope_raster_15, xy = T), aes(x=x, y=y, fill = slope))+
  geom_sf(data = SD_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()

## Extracting the aspect 

#all points 

#extracting the slope in degrees, using the queens method (neighbor = 8)
all_points_aspect_raster_15 <- terra::terrain(CEM_15_utm_all_points, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(all_points_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = fixed_field_data_processed_source_sf_trans_coordinates_source)+
  scale_fill_viridis_c()


#LM

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LM_aspect_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'aspect')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(LM_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LM_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()

#LC

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LC_aspect_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'aspect')

#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(LC_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = LC_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()

#SD

#extracting the aspect in degrees, using the queens method (neighbor = 8)
SD_aspect_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'aspect')


#plot of the slopes
ggplot()+
  geom_raster(data= as.data.frame(SD_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = SD_fixed_field_data_processed_source_source)+
  scale_fill_viridis_c()


#creating dataframes for each population and the slope and aspect data by extracting the slope and aspect data 
# from each cell for each point and combining the data into a dataframe

#all points

all_points_aspect_raster_15_data_pts <- extract(all_points_aspect_raster_15, fixed_field_data_processed_source_sf_trans_coordinates_source) #extracting aspect for each point value
all_points_slope_raster_15_data_pts <- extract(all_points_slope_raster_15, fixed_field_data_processed_source_sf_trans_coordinates_source) #extracting slope for each point value
all_points_fixed_field_data_processed_source_terrain <- cbind(fixed_field_data_processed_source_sf_trans_coordinates_source, all_points_aspect_raster_15_data_pts) #bind the aspect data for each point to the all point dataframe
all_points_fixed_field_data_processed_source_terrain <- cbind(all_points_fixed_field_data_processed_source_terrain, all_points_slope_raster_15_data_pts) #bind the slope data for each point to the all point dataframe

#LM

LM_aspect_raster_15_data_pts <- extract(LM_aspect_raster_15, LM_fixed_field_data_processed_source_source) #extracting aspect for each point value
LM_slope_raster_15_data_pts <- extract(LM_slope_raster_15, LM_fixed_field_data_processed_source_source) #extracting slope for each point value
LM_elevation_raster_15_data_pts <- extract(CEM_15_utm_LM, LM_fixed_field_data_processed_source_source) #extracting the elevation for each point value
LM_fixed_field_data_processed_source_source_terrain <- cbind(LM_fixed_field_data_processed_source_source, LM_aspect_raster_15_data_pts) #bind the aspect data for each point to the LM point dataframe
LM_fixed_field_data_processed_source_source_terrain <- cbind(LM_fixed_field_data_processed_source_source_terrain, LM_slope_raster_15_data_pts) #bind the slope data for each point to the LM point dataframe
LM_fixed_field_data_processed_source_source_terrain <- cbind(LM_fixed_field_data_processed_source_source_terrain, LM_elevation_raster_15_data_pts) #bind the elevation data for each point to the LM point dataframe

#LC

LC_aspect_raster_15_data_pts <- extract(LC_aspect_raster_15, LC_fixed_field_data_processed_source_source) #extracting aspect for each point value
LC_slope_raster_15_data_pts <- extract(LC_slope_raster_15, LC_fixed_field_data_processed_source_source) #extracting slope for each point value
LC_elevation_raster_15_data_pts <- extract(CEM_15_utm_LC, LC_fixed_field_data_processed_source_source) #extracting the elevation for each point value
LC_fixed_field_data_processed_source_source_terrain <- cbind(LC_fixed_field_data_processed_source_source, LC_aspect_raster_15_data_pts) #bind the aspect data for each point to the LC point dataframe
LC_fixed_field_data_processed_source_source_terrain <- cbind(LC_fixed_field_data_processed_source_source_terrain, LC_slope_raster_15_data_pts) #bind the slope data for each point to the LC point dataframe
LC_fixed_field_data_processed_source_source_terrain <- cbind(LC_fixed_field_data_processed_source_source_terrain, LC_elevation_raster_15_data_pts) #bind the elevation data for each point to the LC point dataframe

#SD

SD_aspect_raster_15_data_pts <- extract(SD_aspect_raster_15, SD_fixed_field_data_processed_source_source) #extracting aspect for each point value
SD_slope_raster_15_data_pts <- extract(SD_slope_raster_15, SD_fixed_field_data_processed_source_source) #extracting slope for each point value
SD_elevation_raster_15_data_pts <- extract(CEM_15_utm_SD, SD_fixed_field_data_processed_source_source) #extracting the elevation for each point value
SD_fixed_field_data_processed_source_source_terrain <- cbind(SD_fixed_field_data_processed_source_source, SD_aspect_raster_15_data_pts) #bind the aspect data for each point to the SD point dataframe
SD_fixed_field_data_processed_source_source_terrain <- cbind(SD_fixed_field_data_processed_source_source_terrain, SD_slope_raster_15_data_pts) #bind the slope data for each point to the SD point dataframe
SD_fixed_field_data_processed_source_source_terrain <- cbind(SD_fixed_field_data_processed_source_source_terrain, SD_elevation_raster_15_data_pts) #bind the elevation data for each point to the SD point dataframe

#re-categorizing the aspect data (using either the 4 or 8 cardinal directions)

#setting values of 360 to 0

#all points
all_points_fixed_field_data_processed_source_terrain <- all_points_fixed_field_data_processed_source_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts = case_when((all_points_aspect_raster_15_data_pts == "360") ~  0,
                                                          (all_points_aspect_raster_15_data_pts != "360")~ all_points_aspect_raster_15_data_pts))

#LM
LM_fixed_field_data_processed_source_source_terrain <- LM_fixed_field_data_processed_source_source_terrain %>%
  mutate(LM_aspect_raster_15_data_pts = case_when((LM_aspect_raster_15_data_pts == "360") ~  0,
                                                  (LM_aspect_raster_15_data_pts != "360")~ LM_aspect_raster_15_data_pts))
#LC

LC_fixed_field_data_processed_source_source_terrain <- LC_fixed_field_data_processed_source_source_terrain %>%
  mutate(LC_aspect_raster_15_data_pts = case_when((LC_aspect_raster_15_data_pts == "360") ~  0,
                                                  (LC_aspect_raster_15_data_pts != "360")~ LC_aspect_raster_15_data_pts))


#SD

SD_fixed_field_data_processed_source_source_terrain <- SD_fixed_field_data_processed_source_source_terrain %>%
  mutate(SD_aspect_raster_15_data_pts = case_when((SD_aspect_raster_15_data_pts == "360") ~  0,
                                                  (SD_aspect_raster_15_data_pts != "360")~ SD_aspect_raster_15_data_pts))


# all points

# Using the 8 cardinal directions: North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are by a range of 45 degrees 
all_points_fixed_field_data_processed_source_terrain <- all_points_fixed_field_data_processed_source_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_8_categorical = case_when((all_points_aspect_raster_15_data_pts > 0 & all_points_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 337.5 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N", #359.99999
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 67.5 & all_points_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 112.5 & all_points_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                        (all_points_aspect_raster_15_data_pts >= 157.5 & all_points_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                        (all_points_aspect_raster_15_data_pts >= 202.5 & all_points_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                        (all_points_aspect_raster_15_data_pts >= 247.5 & all_points_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 292.5 & all_points_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# Using the 4 cardinal directions: North, East, South, West

# the directions are by a range of 90 degrees 
all_points_fixed_field_data_processed_source_terrain <- all_points_fixed_field_data_processed_source_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_4_categorical = case_when((all_points_aspect_raster_15_data_pts >= 0 & all_points_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

# LM

# Using the 8 cardinal directions: North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are by a range of 45 degrees 
LM_fixed_field_data_processed_source_source_terrain <- LM_fixed_field_data_processed_source_source_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_8_categorical = case_when((LM_aspect_raster_15_data_pts > 0 & LM_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LM_aspect_raster_15_data_pts >= 337.5 & LM_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 67.5 & LM_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 112.5 & LM_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (LM_aspect_raster_15_data_pts >= 157.5 & LM_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (LM_aspect_raster_15_data_pts >= 202.5 & LM_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (LM_aspect_raster_15_data_pts >= 247.5 & LM_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (LM_aspect_raster_15_data_pts >= 292.5 & LM_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# Using the 4 cardinal directions: North, East, South, West

# the directions are by a range of 90 degrees 
LM_fixed_field_data_processed_source_source_terrain <- LM_fixed_field_data_processed_source_source_terrain %>%
  mutate(LM_aspect_raster_15_data_pts_4_categorical = case_when((LM_aspect_raster_15_data_pts >= 0 & LM_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LM_aspect_raster_15_data_pts >= 315 & LM_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LM_aspect_raster_15_data_pts >= 22.5 & LM_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LM_aspect_raster_15_data_pts >= 135 & LM_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LM_aspect_raster_15_data_pts >= 225 & LM_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315


# LC

# Using the 8 cardinal directions: North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are by a range of 45 degrees 
LC_fixed_field_data_processed_source_source_terrain <- LC_fixed_field_data_processed_source_source_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_8_categorical = case_when((LC_aspect_raster_15_data_pts > 0 & LC_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 337.5 & LC_aspect_raster_15_data_pts < 359.999999999999999) ~ "N", 
                                                                (LC_aspect_raster_15_data_pts >= 22.5 & LC_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 67.5 & LC_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 112.5 & LC_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (LC_aspect_raster_15_data_pts >= 157.5 & LC_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (LC_aspect_raster_15_data_pts >= 202.5 & LC_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (LC_aspect_raster_15_data_pts >= 247.5 & LC_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (LC_aspect_raster_15_data_pts >= 292.5 & LC_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees


# Using the 4 cardinal directions: North, East, South, West

# the directions are by a range of 90 degrees 
LC_fixed_field_data_processed_source_source_terrain <- LC_fixed_field_data_processed_source_source_terrain %>%
  mutate(LC_aspect_raster_15_data_pts_4_categorical = case_when((LC_aspect_raster_15_data_pts >= 0 & LC_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (LC_aspect_raster_15_data_pts >= 315 & LC_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (LC_aspect_raster_15_data_pts >= 45 & LC_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (LC_aspect_raster_15_data_pts >= 135 & LC_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (LC_aspect_raster_15_data_pts >= 225 & LC_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

#SD
# Using the 8 cardinal directions: North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are by a range of 45 degrees 
SD_fixed_field_data_processed_source_source_terrain <- SD_fixed_field_data_processed_source_source_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_8_categorical = case_when((SD_aspect_raster_15_data_pts > 0 & SD_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 337.5 & SD_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (SD_aspect_raster_15_data_pts >= 22.5 & SD_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 67.5 & SD_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 112.5 & SD_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                (SD_aspect_raster_15_data_pts >= 157.5 & SD_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                (SD_aspect_raster_15_data_pts >= 202.5 & SD_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                (SD_aspect_raster_15_data_pts >= 247.5 & SD_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                (SD_aspect_raster_15_data_pts >= 292.5 & SD_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# Using the 4 cardinal directions: North, East, South, West

# the directions are by a range of 90 degrees 
SD_fixed_field_data_processed_source_source_terrain <- SD_fixed_field_data_processed_source_source_terrain %>%
  mutate(SD_aspect_raster_15_data_pts_4_categorical = case_when((SD_aspect_raster_15_data_pts >= 0 & SD_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                (SD_aspect_raster_15_data_pts >= 315 & SD_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                (SD_aspect_raster_15_data_pts >= 45 & SD_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                (SD_aspect_raster_15_data_pts >= 135 & SD_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                (SD_aspect_raster_15_data_pts >= 225 & SD_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

#### Creating the distance to river columns ####

#LM

#turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
river_LM_trans_points <- st_cast(river_LM_trans_source, "LINESTRING") #turning the polyline of the river into a linestring object
river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #creating a raster out of the river linestring object
plot(river_LM_trans_point_raster) #plotting the river linestring object

#turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
river_LM_buffer_trans_outline <- st_cast(river_buffer_LM_source, "LINESTRING") #turning the polygon of the river buffer into a linestring object
river_buffer_LM_point_raster <- st_rasterize(river_LM_buffer_trans_outline) #creating a raster of river buffer linestring object
plot(river_buffer_LM_point_raster) #plotting the river buffer linestring object

#generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the cell that are not the river buffer linestring raster have a 0 value
dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, river_LM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
plot(dist_near_river_buffer_LM) #plotting the distance to river raster

#LC

#turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
river_LC_trans_points <- st_cast(river_LC_trans_source, "LINESTRING") #turning the polyline of the river into a linestring object
river_LC_trans_point_raster <- st_rasterize(river_LC_trans_points) #creating a raster out of the river linestring object
plot(river_LC_trans_point_raster)

#turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
river_buffer_LC_points <- st_cast(river_buffer_LC_source, "LINESTRING") #turning the polygon of the river buffer into a linestring object
river_buffer_LC_point_raster <- st_rasterize(river_buffer_LC_points) #creating a raster of river buffer linestring object
plot(river_buffer_LC_point_raster) #plotting the river buffer linestring object

#generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
river_buffer_LC_point_raster[is.na(river_buffer_LC_point_raster[])] <- 0  #making sure the cells that are not part of the the river buffer raster have a 0 value
dist_near_river_buffer_LC <- dist_to_nearest(river_buffer_LC_point_raster, river_LC_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
plot(dist_near_river_buffer_LC) #not using inverse distance

#SD

#turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
river_SD_trans_points <- st_cast(river_SD_trans_source, "LINESTRING") #turning the polyline of the river into a linestring object
river_SD_trans_point_raster <- st_rasterize(river_SD_trans_points) #creating a raster out of the river linestring object
plot(river_SD_trans_point_raster)

#turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
river_buffer_SD_points <- st_cast(river_buffer_SD_source, "LINESTRING") #turning the polygon of the river buffer into a linestring object
river_buffer_SD_point_raster <- st_rasterize(river_buffer_SD_points) #creating a raster of river buffer linestring object
plot(river_buffer_SD_point_raster) #plotting the river buffer linestring object

#generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
river_buffer_SD_point_raster[is.na(river_buffer_SD_point_raster[])] <- 0  #making sure the cells that are not part of the the river buffer raster have a 0 value
dist_near_river_buffer_SD <- dist_to_nearest(river_buffer_SD_point_raster, river_SD_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
plot(dist_near_river_buffer_SD) #plotting the distance to river raster


## Making it so the cells in the distance raster within or overlapping with the river raster are assigned 1 

#LM

#Assigning points within and overlapping with the river to be "true"
LM_points_intersects_river <- st_intersects(LM_fixed_field_data_processed_source_source, river_LM_trans_source, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
LM_fixed_field_data_processed_source_source_intersects_river <- cbind(LM_fixed_field_data_processed_source_source, LM_points_intersects_river) #binding the list of true or falses with the point data
#printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
ggplot()+
  geom_sf(data=river_LM_trans_source)+
  geom_sf(data=LM_fixed_field_data_processed_source_source)+
  geom_sf(data=LM_fixed_field_data_processed_source_source_intersects_river, aes(color = LM_points_intersects_river))

#LC 

#Assigning points within and overlapping with the river to be "true"
LC_points_intersects_river <- st_intersects(LC_fixed_field_data_processed_source_source, river_LC_trans_source, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
LC_fixed_field_data_processed_source_source_intersects_river <- cbind(LC_fixed_field_data_processed_source_source, LC_points_intersects_river) #binding the list of true or falses with the point data
#printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
ggplot()+
  geom_sf(data=river_LC_trans_source)+
  geom_sf(data=LC_fixed_field_data_processed_source_source)+
  geom_sf(data=LC_fixed_field_data_processed_source_source_intersects_river, aes(color = LC_points_intersects_river))

#SD

#Assigning points within and overlapping with the river to be "true"
SD_points_intersects_river <- st_intersects(SD_fixed_field_data_processed_source_source, river_SD_trans_source, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
SD_fixed_field_data_processed_source_source_intersects_river <- cbind(SD_fixed_field_data_processed_source_source, SD_points_intersects_river) #binding the list of true or falses with the point data
#printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
ggplot()+
  geom_sf(data=river_SD_trans_source)+
  geom_sf(data=SD_fixed_field_data_processed_source_source)+
  geom_sf(data=SD_fixed_field_data_processed_source_source_intersects_river, aes(color = SD_points_intersects_river))

## Extracting distance to river for each tree using the distance to river raster

#LM
LM_distance_data_pts <- st_extract(dist_near_river_buffer_LM, LM_fixed_field_data_processed_source_source) #extracting distance to river for each tree
LM_fixed_field_data_processed_source_source_distance  <- cbind(LM_fixed_field_data_processed_source_source, LM_distance_data_pts) #binding the distance to river data for each point to the LM point dataframe


#LC
LC_distance_data_pts <- st_extract(dist_near_river_buffer_LC, LC_fixed_field_data_processed_source_source) #extracting distance to river for each tree
LC_fixed_field_data_processed_source_source_distance  <- cbind(LC_fixed_field_data_processed_source_source, LC_distance_data_pts) #binding the distance to river data for each tree to the LC point dataframe


#SD
SD_distance_data_pts <- st_extract(dist_near_river_buffer_SD, SD_fixed_field_data_processed_source_source) #extracting distance to river for each tree
SD_fixed_field_data_processed_source_source_distance  <- cbind(SD_fixed_field_data_processed_source_source, SD_distance_data_pts) #binding the distance to river data for each point to the SD point dataframe

## Assigning all points within/overlapping river to distances of 0

#LM
LM_fixed_field_data_processed_source_source_distance <- LM_fixed_field_data_processed_source_source_distance %>% 
  mutate(d = case_when((LM_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LM_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 

#LC
LC_fixed_field_data_processed_source_source_distance <- LC_fixed_field_data_processed_source_source_distance %>%
  mutate(d = case_when((LC_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LC_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 

#SD
SD_fixed_field_data_processed_source_source_distance <- SD_fixed_field_data_processed_source_source_distance %>%
  mutate(d = case_when((SD_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (SD_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value


#### Combining the elevation, slope, aspect, and distance data into one dataframe ####

#LM
LM_fixed_field_data_processed_source_source_terrain_dist <- cbind(LM_fixed_field_data_processed_source_source_terrain, LM_fixed_field_data_processed_source_source_distance)

#LC
LC_fixed_field_data_processed_source_source_terrain_dist <- cbind(LC_fixed_field_data_processed_source_source_terrain, LC_fixed_field_data_processed_source_source_distance)

#SD
SD_fixed_field_data_processed_source_source_terrain_dist <- cbind(SD_fixed_field_data_processed_source_source_terrain, SD_fixed_field_data_processed_source_source_distance)


#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + # Short Canopy Axis
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + # Long Canopy Axis
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + # Crown Spread
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + # Canopy Area
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")


ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + # DBH
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#LM

ggplot(LM_fixed_field_data_processed_source_source) + # Short Canopy Axis
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed_source_source) + # Long Canopy Axis
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed_source_source) + # Crown Spread
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed_source_source) + # Canopy Area
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_processed_source_source) + # DBH
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#LC

ggplot(LC_fixed_field_data_processed_source_source) + # Short Canopy Axis
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_source_source) + # Long Canopy Axis
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_source_source) + # Crown Spread
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_source_source) + # Canopy Area
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_processed_source_source) + # DBH
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed_source_source) + # Short Canopy Axis
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_source_source) + # Long Canopy Axis
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_source_source) + #  Crown Spread
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_source_source) + # Canopy Area
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_processed_source_source) + # DBH
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")


# Elevation histograms to look at the spread of the data

# all points
ggplot(fixed_field_data_processed_source_sf_trans_coordinates_source) + 
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation")+
  ylab("Frequency")

# LM
ggplot(LM_fixed_field_data_processed_source_source) + 
  geom_histogram(aes(x = Elevation..m.FIXED))+
  xlab("Elevation (m)")+
  ylab("Frequency")

#LC
ggplot(LC_fixed_field_data_processed_source_source) + 
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation (m)")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed_source_source) + 
  geom_histogram(aes(x = Elevation..m.))+
  xlab("Elevation (m)")+
  ylab("Frequency")

#histograms for slope

#LM
ggplot(LM_fixed_field_data_processed_source_source_terrain) + 
  geom_histogram(aes(x = LM_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

#LC
ggplot(LC_fixed_field_data_processed_source_source_terrain) + 
  geom_histogram(aes(x = LC_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed_source_source_terrain) + 
  geom_histogram(aes(x = SD_slope_raster_15_data_pts))+
  xlab("Slope (degrees)")+
  ylab("Frequency")

# barcharts for looking at the spread of aspect

# 8 categories of direction

#LM
ggplot(LM_fixed_field_data_processed_source_source_terrain) + 
  geom_bar(aes(x = LM_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

#LC
ggplot(LC_fixed_field_data_processed_source_source_terrain) + 
  geom_bar(aes(x = LC_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed_source_source_terrain) + 
  geom_bar(aes(x = SD_aspect_raster_15_data_pts_8_categorical))+
  xlab("Direction")+
  ylab("Frequency")

# 4 categories of direction

#LM
ggplot(LM_fixed_field_data_processed_source_source_terrain) + #generate the base plot
  geom_bar(aes(x = LM_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

#LC
ggplot(LC_fixed_field_data_processed_source_source_terrain) + #generate the base plot
  geom_bar(aes(x = LC_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")

#SD
ggplot(SD_fixed_field_data_processed_source_source_terrain) + #generate the base plot
  geom_bar(aes(x = SD_aspect_raster_15_data_pts_4_categorical))+
  xlab("Direction")+
  ylab("Frequency")


#Summary Statistics
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed_source %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)


#### Load in Soil Metric Data ####

# loading in soil textures raster tifs from CONABIO

clay_05 <- raster(paste0("./data/Soil Grid/clay content/clay content 0-5.tif"))
clay_200 <- raster(paste0("./data/Soil Grid/clay content/clay content 100-200.tif"))
silt_05 <- raster(paste0("./data/Soil Grid/silt/silt 0-5.tif"))
silt_200 <-raster(paste0("./data/Soil Grid/silt/silt 100-200.tif"))
sand_05 <- raster(paste0("./data/Soil Grid/sand/sand 0-5.tif"))
sand_200 <- raster(paste0("./data/Soil Grid/sand/sand 100-200.tif"))
ph_05 <- raster(paste0("./data/Soil Grid/pH/ph_0-5.tif")) #0-5 cm ph
ph_200 <- raster(paste0("./data/Soil Grid/pH/ph_100-200.tif")) #100-200 ph
ocd_05 <- raster(paste0("./data/Soil Grid/organic carbon density/OCD_0-5.tif")) #0-5cm organic carbon density
ocd_200 <- raster(paste0("./data/Soil Grid/organic carbon density/OCR_100-200.tif")) #100-200cm organic carbon density
coarse_frag_05 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_0-5.tif")) #0-5 cm coarse fragments
coarse_frag_200 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_100-200.tif")) #100-200 cm coarse fragments
cat_ex_cap_05 <-raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_0-5.tif")) #0-5 cm cation exchange capacity
cat_ex_cap_200 <- raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_100-200.tif")) #100-200 cm cation exchange capacity
bulk_dens_05 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_0-5.tif")) #0-5 cm bulk density
bulk_dens_200 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_100-200.tif")) #100-200 cm bulk density
vol_wat_10kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_0-5.tif"))  #0-5 cm -10 kpa volumn water content
vol_wat_10kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_100-200.tif"))  #100-200 cm -10 kpa volumn water content
vol_wat_33kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -33 kPa /vol_water_0-5.tif")) #0-5 cm -33 kpa volumn water content
vol_wat_33kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -33 kPa /vol_water_100-200.tif")) #100-200 cm -33 kpa volumn water content
vol_wat_1500kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500kPa_0-5.tif"))  #0-5 cm -1500 kpa volumn water content
vol_wat_1500kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500_100-200.tif")) #100-200 cm -1500 kpa volumn water content
nitrogen_05 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 0-5.tif")) #0-5 cm nitrogen content 
nitrogen_200 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 100-200.tif")) #100-200 cm nitrogen content 
Soil_Organic_Carbon_05 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 0-5.tif")) #0-5 cm soil organic carbon
Soil_Organic_Carbon_200 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 100-200.tif")) #100-200 cm soil organic carbon

#project the soil metric rasters to equal area projection (UTM 12N) which uses meters as distance measurement 

clay_05_utm <- projectRaster(clay_05, crs=26912) #converting the 0-5 cm clay raster to utm 12
clay_200_utm <- projectRaster(clay_200, crs=26912) #converting the 90-200 cm clay raster to utm 12
silt_05_utm <- projectRaster(silt_05, crs=26912)
silt_200_utm <- projectRaster(silt_200, crs=26912)
sand_05_utm <- projectRaster(sand_05, crs=26912)
sand_200_utm <- projectRaster(sand_200, crs=26912)
ph_05_utm <- projectRaster(ph_05, crs=26912) 
ph_200_utm <- projectRaster(ph_200, crs=26912) 
ocd_05_utm <- projectRaster(ocd_05, crs=26912)
ocd_200_utm <- projectRaster(ocd_200, crs=26912)
coarse_frag_05_utm <- projectRaster(coarse_frag_05, crs=26912)
coarse_frag_200_utm <- projectRaster(coarse_frag_200, crs=26912)
cat_ex_cap_05_utm <- projectRaster(cat_ex_cap_05, crs=26912)
cat_ex_cap_200_utm <- projectRaster(cat_ex_cap_200, crs=26912)
bulk_dens_05_utm <- projectRaster(bulk_dens_05, crs=26912)
bulk_dens_200_utm <- projectRaster(bulk_dens_200, crs=26912)
vol_wat_10kpa_05_utm <- projectRaster(vol_wat_10kpa_05, crs=26912)
vol_wat_10kpa_200_utm <- projectRaster(vol_wat_10kpa_200, crs=26912)
vol_wat_33kpa_05_utm <- projectRaster(vol_wat_33kpa_05, crs=26912)
vol_wat_33kpa_200_utm <- projectRaster(vol_wat_33kpa_200, crs=26912)
vol_wat_1500kpa_05_utm <- projectRaster(vol_wat_1500kpa_05, crs=26912)
vol_wat_1500kpa_200_utm <- projectRaster(vol_wat_1500kpa_200, crs=26912)
nitrogen_05_utm <- projectRaster(nitrogen_05, crs=26912)
nitrogen_200_utm <- projectRaster(nitrogen_200, crs=26912)
Soil_Organic_Carbon_05_utm <- projectRaster(Soil_Organic_Carbon_05, crs=26912)
Soil_Organic_Carbon_200_utm <- projectRaster(Soil_Organic_Carbon_200, crs=26912)

## cropping the soil metrics to bounding boxes around each population ##

# LM

# using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_LM <- crop(clay_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100)) 
clay_200_LM <- crop(clay_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
silt_05_LM <- crop(silt_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
silt_200_LM <- crop(silt_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
sand_05_LM <- crop(sand_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
sand_200_LM <- crop(sand_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
ph_05_LM <- crop(ph_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
ph_200_LM <- crop(ph_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
ocd_05_LM <- crop(ocd_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
ocd_200_LM <- crop(ocd_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
coarse_frag_05_LM <- crop(coarse_frag_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
coarse_frag_200_LM <- crop(coarse_frag_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
cat_ex_cap_05_LM <- crop(cat_ex_cap_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
cat_ex_cap_200_LM <- crop(cat_ex_cap_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
bulk_dens_05_LM <- crop(bulk_dens_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
bulk_dens_200_LM <- crop(bulk_dens_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_10kpa_05_LM <- crop(vol_wat_10kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_10kpa_200_LM <- crop(vol_wat_10kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_33kpa_05_LM <- crop(vol_wat_33kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_33kpa_200_LM <- crop(vol_wat_33kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_1500kpa_05_LM <- crop(vol_wat_1500kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
vol_wat_1500kpa_200_LM <- crop(vol_wat_1500kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
nitrogen_05_LM <-  crop(nitrogen_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
nitrogen_200_LM <- crop(nitrogen_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
Soil_Organic_Carbon_05_LM <- crop(Soil_Organic_Carbon_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
Soil_Organic_Carbon_200_LM <- crop(Soil_Organic_Carbon_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))


# plotting an example soil metric raster with river shapefile and tree points  
ggplot()+
  geom_raster(data = as.data.frame(Soil_Organic_Carbon_05_LM, xy=T), aes(x=x, y=y, fill = SOC.0.5))+
  geom_sf(data = river_LM_trans_source)+
  geom_sf(data = LM_fixed_field_data_processed_source_source)

#creating stacks of the soil raster layers (using multiple for visualization purposes)
soil_stack_LM_soil_text <- stack(clay_05_LM, clay_200_LM, silt_05_LM, silt_200_LM, sand_05_LM, sand_200_LM) #the stack of all of the soil texture rasters
soil_stack_LM_other <- stack(ph_05_LM, ph_200_LM, ocd_05_LM, ocd_200_LM, coarse_frag_05_LM, coarse_frag_200_LM, #the stack of all of the other soil variables, with different extents than the soil texture rasters
                             cat_ex_cap_05_LM, cat_ex_cap_200_LM, bulk_dens_05_LM, bulk_dens_200_LM, vol_wat_10kpa_05_LM,
                             vol_wat_10kpa_200_LM, vol_wat_33kpa_05_LM, vol_wat_33kpa_200_LM, vol_wat_1500kpa_05_LM, 
                             vol_wat_1500kpa_200_LM) 
soil_stack_LM_extra <- stack(nitrogen_05_LM, nitrogen_200_LM, Soil_Organic_Carbon_05_LM, Soil_Organic_Carbon_200_LM)



#plotting the stacked rasters
plot(soil_stack_LM_soil_text) #version with soil textures
plot(soil_stack_LM_soil_text, zlim = c(100, 710)) #version where the plots have the same scale
plot(soil_stack_LM_other) #version with other variables
plot(soil_stack_LM_other, zlim = c(30, 360)) #version where the plots have the same scale
plot(soil_stack_LM_extra) #version with other variables
plot(soil_stack_LM_extra, zlim = c(30, 360)) #version where the plots have the same scale


#LC

# using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_LC <- crop(clay_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100)) 
clay_200_LC <- crop(clay_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
silt_05_LC <- crop(silt_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
silt_200_LC <- crop(silt_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
sand_05_LC <- crop(sand_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
sand_200_LC <- crop(sand_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
ph_05_LC <- crop(ph_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
ph_200_LC <- crop(ph_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
ocd_05_LC <- crop(ocd_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
ocd_200_LC <- crop(ocd_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
coarse_frag_05_LC <- crop(coarse_frag_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
coarse_frag_200_LC <- crop(coarse_frag_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
cat_ex_cap_05_LC <- crop(cat_ex_cap_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
cat_ex_cap_200_LC <- crop(cat_ex_cap_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
bulk_dens_05_LC <- crop(bulk_dens_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
bulk_dens_200_LC <- crop(bulk_dens_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_10kpa_05_LC <- crop(vol_wat_10kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_10kpa_200_LC <- crop(vol_wat_10kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_33kpa_05_LC <- crop(vol_wat_33kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_33kpa_200_LC <- crop(vol_wat_33kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_1500kpa_05_LC <- crop(vol_wat_1500kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
vol_wat_1500kpa_200_LC <- crop(vol_wat_1500kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
nitrogen_05_LC <-  crop(nitrogen_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
nitrogen_200_LC <- crop(nitrogen_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
Soil_Organic_Carbon_05_LC <- crop(Soil_Organic_Carbon_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
Soil_Organic_Carbon_200_LC <- crop(Soil_Organic_Carbon_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))

#creating stacks of the soil raster layers (using multiple for visualization purposes)
soil_stack_LC_soil_text <- stack(clay_05_LC, clay_200_LC, silt_05_LC, silt_200_LC, sand_05_LC, sand_200_LC) #the stack of all of the soil texture rasters
soil_stack_LC_other <- stack(ph_05_LC, ph_200_LC, ocd_05_LC, ocd_200_LC, coarse_frag_05_LC, coarse_frag_200_LC, #the stack of all of the other soil variables, with different extents than the soil texture rasters
                             cat_ex_cap_05_LC, cat_ex_cap_200_LC, bulk_dens_05_LC, bulk_dens_200_LC, vol_wat_10kpa_05_LC,
                             vol_wat_10kpa_200_LC, vol_wat_33kpa_05_LC, vol_wat_33kpa_200_LC, vol_wat_1500kpa_05_LC, 
                             vol_wat_1500kpa_200_LC) 
soil_stack_LC_extra <- stack(nitrogen_05_LC, nitrogen_200_LC, Soil_Organic_Carbon_05_LC, Soil_Organic_Carbon_200_LC)

#plotting the stacked rasters
plot(soil_stack_LC_soil_text) #version with soil textures
plot(soil_stack_LC_soil_text, zlim = c(100, 710)) #version where the plots have the same scale
plot(soil_stack_LC_other) #version with other variables
plot(soil_stack_LC_other, zlim = c(30, 360)) #version where the plots have the same scale
plot(soil_stack_LC_extra) #version with other variables
plot(soil_stack_LC_extra, zlim = c(30, 180)) #version where the plots have the same scale


#SD

#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_SD <- crop(clay_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100)) 
clay_200_SD <- crop(clay_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
silt_05_SD <- crop(silt_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
silt_200_SD <- crop(silt_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
sand_05_SD <- crop(sand_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
sand_200_SD <- crop(sand_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
ph_05_SD <- crop(ph_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
ph_200_SD <- crop(ph_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
ocd_05_SD <- crop(ocd_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
ocd_200_SD <- crop(ocd_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
coarse_frag_05_SD <- crop(coarse_frag_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
coarse_frag_200_SD <- crop(coarse_frag_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
cat_ex_cap_05_SD <- crop(cat_ex_cap_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
cat_ex_cap_200_SD <- crop(cat_ex_cap_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
bulk_dens_05_SD <- crop(bulk_dens_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
bulk_dens_200_SD <- crop(bulk_dens_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_10kpa_05_SD <- crop(vol_wat_10kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_10kpa_200_SD <- crop(vol_wat_10kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_33kpa_05_SD <- crop(vol_wat_33kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_33kpa_200_SD <- crop(vol_wat_33kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_1500kpa_05_SD <- crop(vol_wat_1500kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
vol_wat_1500kpa_200_SD <- crop(vol_wat_1500kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
nitrogen_05_SD <-  crop(nitrogen_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
nitrogen_200_SD <- crop(nitrogen_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
Soil_Organic_Carbon_05_SD <- crop(Soil_Organic_Carbon_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
Soil_Organic_Carbon_200_SD <- crop(Soil_Organic_Carbon_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))

#creating stacks of the soil raster layers (using multiple for visualization purposes)
soil_stack_SD_soil_text <- stack(clay_05_SD, clay_200_SD, silt_05_SD, silt_200_SD, sand_05_SD, sand_200_SD) #the stack of all of the soil texture rasters
soil_stack_SD_other <- stack(ph_05_SD, ph_200_SD, ocd_05_SD, ocd_200_SD, coarse_frag_05_SD, coarse_frag_200_SD, #the stack of all of the other soil variables, with different extents than the soil texture rasters
                             cat_ex_cap_05_SD, cat_ex_cap_200_SD, bulk_dens_05_SD, bulk_dens_200_SD, vol_wat_10kpa_05_SD,
                             vol_wat_10kpa_200_SD, vol_wat_33kpa_05_SD, vol_wat_33kpa_200_SD, vol_wat_1500kpa_05_SD, 
                             vol_wat_1500kpa_200_SD) 
soil_stack_SD_extra <- stack(nitrogen_05_SD, nitrogen_200_SD,Soil_Organic_Carbon_05_SD,  Soil_Organic_Carbon_200_SD)

#plotting the stacked rasters
plot(soil_stack_SD_soil_text)
plot(soil_stack_SD_soil_text, zlim = c(130, 710)) #version where the plots have the same scale
plot(soil_stack_SD_other)
plot(soil_stack_SD_other, zlim = c(45, 360)) #version where the plots have the same scale
plot(soil_stack_SD_extra)
plot(soil_stack_SD_extra, zlim = c(25, 340)) #version where the plots have the same scale

#creating X sequential columns in LC and SD point data which will make it easier to select random points from each grid later

#creating an x_sequential column that is 1 through the number of LM points, which will make it easier to randomly choose one point
LM_fixed_field_data_processed_source_source <- LM_fixed_field_data_processed_source_source %>%
  mutate(X_sequential = 1:nrow(LM_fixed_field_data_processed_source_source))

#creating an x_sequential column that is 1 through the number of LC points, which will make it easier to randomly choose one point
LC_fixed_field_data_processed_source_source <- LC_fixed_field_data_processed_source_source %>%
  mutate(X_sequential = 1:nrow(LC_fixed_field_data_processed_source_source))

#creating an x_sequential column that is 1 through the number of SD points, which will make it easier to randomly choose one point
SD_fixed_field_data_processed_source_source <- SD_fixed_field_data_processed_source_source %>%
  mutate(X_sequential = 1:nrow(SD_fixed_field_data_processed_source_source))

## Extracting the soil data to the tree points 

#LM
LM_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LM_soil_text, LM_fixed_field_data_processed_source_source) #extracting soil textures for each point value
LM_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LM_other, LM_fixed_field_data_processed_source_source) #extracting the other soil variables for each point value
LM_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LM_extra, LM_fixed_field_data_processed_source_source) #extracting the extra soil variables for each point value
LM_fixed_field_data_processed_source_source_soils <- cbind(LM_fixed_field_data_processed_source_source, LM_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LM point dataframe
LM_fixed_field_data_processed_source_source_soils <- cbind(LM_fixed_field_data_processed_source_source_soils, LM_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LM point dataframe
LM_fixed_field_data_processed_source_source_soils <- cbind(LM_fixed_field_data_processed_source_source_soils, LM_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LM point dataframe

#LC
LC_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LC_soil_text, LC_fixed_field_data_processed_source_source) #extracting soil textures for each point value
LC_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LC_other, LC_fixed_field_data_processed_source_source) #extracting the other soil variables for each point value
LC_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LC_extra, LC_fixed_field_data_processed_source_source) #extracting the extra soil variables for each point value
LC_fixed_field_data_processed_source_source_soils <- cbind(LC_fixed_field_data_processed_source_source, LC_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
LC_fixed_field_data_processed_source_source_soils <- cbind(LC_fixed_field_data_processed_source_source_soils, LC_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
LC_fixed_field_data_processed_source_source_soils <- cbind(LC_fixed_field_data_processed_source_source_soils, LC_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe

#SD
SD_soil_text_raster_250_data_pts <- raster::extract(soil_stack_SD_soil_text, SD_fixed_field_data_processed_source_source) #extracting soil textures for each point value
SD_soil_other_raster_250_data_pts <- raster::extract(soil_stack_SD_other, SD_fixed_field_data_processed_source_source) #extracting the other soil variables for each point value
SD_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_SD_extra, SD_fixed_field_data_processed_source_source) #extracting the extra soil variables for each point value
SD_fixed_field_data_processed_source_source_soils <- cbind(SD_fixed_field_data_processed_source_source, SD_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
SD_fixed_field_data_processed_source_source_soils <- cbind(SD_fixed_field_data_processed_source_source_soils, SD_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
SD_fixed_field_data_processed_source_source_soils <- cbind(SD_fixed_field_data_processed_source_source_soils, SD_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe


#### Creating Sandy and Clay/Loamy Available Water Columns ####

# Making four new soil metric variables for soil available water, which equals the field capacity - permanent wilting point

# LM

LM_fixed_field_data_processed_source_source_soils <- LM_fixed_field_data_processed_source_source_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm

# LC

LC_fixed_field_data_processed_source_source_soils <- LC_fixed_field_data_processed_source_source_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200)  # Clay/Loam Available Water 100-200 cm

# SD

SD_fixed_field_data_processed_source_source_soils <- SD_fixed_field_data_processed_source_source_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm
