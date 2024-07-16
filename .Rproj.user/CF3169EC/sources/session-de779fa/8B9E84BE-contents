#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster)

library(gdalUtils) #to reproject large rasters faster
library(rgdal) #needed to get the soil data
library(XML) #needed to get the soil data

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

#set elevation as a numeric value
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

#creating a new elevation column so the values that were mistakenly 
LM_fixed_field_data_processed <-  LM_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.))
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

#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")


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



##Load in environmental rasters 

#####soil layers... ####
#downloading soil rasters from soil grid
#below from https://git.wur.nl/isric/soilgrids/soilgrids.notebooks/-/blob/master/markdown/webdav_from_R.md... works for ghana.. doesn't work when I try to change the bounding box to anything that isn't on thier website... it then fails by telling me the computed -srcwin has a negative width or height... and computes the same srcwin every time... not sure why.... WORKED WHEN I DONT USE A BBOX AND DOWNLOAD ALL TEH DATA (admittedly huge file but can be deleted once I have my cropped one so :) )
####ph!
voi = "Organic carbon density" # variable of interest
depth = "5-5"
Value = "mean" #mean
voi_layer = paste(voi,depth,quantile, sep="_") # layer of interest 

sg_url= "https://soilgrids.org/"
gdal_translate(paste0(sg_url, voi, "/", voi_layer,'.vrt'),
               "./raster_data/phh2o_igh_r.tif",
               verbose=TRUE)


#loading in soil textures from CONABIO
clay_05 <- raster(paste0("./data/Soil Textur Geotiff geographic coordinates /cly_05cm_mgw/cly_05cm_mgw.tif"))
clay_200 <- raster(paste0("./data/Soil Textur Geotiff geographic coordinates /cly_200cm_mgw/cly_200cm_mgw.tif"))
silt_05 <- raster(paste0("./data/Soil Textur Geotiff geographic coordinates /slt_05cm_pgw/slt_05cm_pgw.tif"))
silt_200 <-raster(paste0("./data/Soil Textur Geotiff geographic coordinates /slt_200cm_pgw/slt_200cm_pgw.tif"))
sand_05 <- raster(paste0("./data/Soil Textur Geotiff geographic coordinates /snd_05cm_mgw/snd_05cm_mgw.tif"))
sand_200 <- raster(paste0("./data/Soil Textur Geotiff geographic coordinates /snd_200cm_mgw/snd_200cm_mgw.tif"))

#project rasters to equal area projection (UTM 12N), uses meters as distance measurement 
clay_05_utm <- projectRaster(clay_05, crs=26912) #converting the 0-5 cm clay raster to utm 12
clay_200_utm <- projectRaster(clay_200, crs=26912) #converting the 90-200 cm clay raster to utm 12
silt_05_utm <- projectRaster(silt_05, crs=26912)
silt_200_utm <- projectRaster(silt_200, crs=26912)
sand_05_utm <- projectRaster(sand_05, crs=26912)
sand_200_utm <- projectRaster(sand_200, crs=26912)


#LM
#examining the layers at different extents
plot(clay_05_utm)
plot(clay_05_LM)
clay_05_LM <- crop(clay_05_utm, extent(LM_box[1]-10000, LM_box[3]+10000, LM_box[2]-10000, LM_box[4]+10000)) 



#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_LM <- crop(clay_05_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4])) 
clay_200_LM <- crop(clay_200_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4]))
silt_05_LM <- crop(silt_05_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4]))
silt_200_LM <- crop(silt_200_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4]))
sand_05_LM <- crop(sand_05_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4]))
sand_200_LM <- crop(sand_200_utm, extent(LM_box[1], LM_box[3], LM_box[2], LM_box[4]))

#attempt of using ggplot to plot clay layer with river shapefile
ggplot()+
  geom_raster(data = as.data.frame(clay_05_LM, xy=T), aes(x=x, y=y, fill = cly_05cm_mgw))+
  geom_sf(data = river_LM_trans)


#creating a stack of the raster layers
soil_stack_LM <- stack(clay_05_LM, clay_200_LM, silt_05_LM, silt_200_LM, sand_05_LM, sand_200_LM)
soil_stack_LM.df <- as.data.frame(getValues(soil_stack_LM))

#attempting to plot the stacks in ggplot
ggplot()+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = cly_05cm_mgw))+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = cly_200cm_mgw))+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = slt_05cm_pgw))+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = slt_200cm_pgw))+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = snd_05cm_mgw))+
  geom_raster(data = as.data.frame(soil_stack_LM, xy=T), aes(x=x, y=y, fill = snd_200cm_mgw))+
  facet_wrap( ~c(cly_05cm_mgw, cly_200cm_mgw, slt_05cm_pgw, slt_200cm_pgw, snd_05cm_mgw, snd_200cm_mgw))
plot(soil_stack_LM, zlim = c(2, 25)) #plotting the soil raster so all of the layers have the same scale
scale_fill_gradientn(colours = terrain.colors(7))


#LC
#using the extent of the box around the rivers to crop the raster for each soil texture layer
#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_LC <- crop(clay_05_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4])) 
clay_200_LC <- crop(clay_200_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4]))
silt_05_LC <- crop(silt_05_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4]))
silt_200_LC <- crop(silt_200_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4]))
sand_05_LC <- crop(sand_05_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4]))
sand_200_LC <- crop(sand_200_utm, extent(LC_box[1], LC_box[3], LC_box[2], LC_box[4]))

#creating a stack of the raster layers
soil_stack_LC <- stack(clay_05_LC, clay_200_LC, silt_05_LC, silt_200_LC, sand_05_LC, sand_200_LC)
plot(soil_stack_LC)


#SD
#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_SD <- crop(clay_05_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4])) 
clay_200_SD <- crop(clay_200_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4]))
silt_05_SD <- crop(silt_05_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4]))
silt_200_SD <- crop(silt_200_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4]))
sand_05_SD <- crop(sand_05_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4]))
sand_200_SD <- crop(sand_200_utm, extent(SD_box[1], SD_box[3], SD_box[2], SD_box[4]))

#creating a stack of the raster layers
soil_stack_SD <- stack(clay_05_SD, clay_200_SD, silt_05_SD, silt_200_SD, sand_05_SD, sand_200_SD)
plot(soil_stack_SD)

