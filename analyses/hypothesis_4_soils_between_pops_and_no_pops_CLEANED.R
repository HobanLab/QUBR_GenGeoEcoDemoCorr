#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

#loading in the data
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

#add average nearest neighbor for each individual column
fixed_field_data_processed_NN_UTM <- fixed_field_data_processed_sf_trans_coordinates %>%  #creates a dataframe with the ANN of the closest 5 individual trees for each individual
  mutate(dist1 = nndist(X = X.1, Y= Y, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
  mutate(dist2 = nndist(X = X.1, Y= Y, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
  mutate(dist3 = nndist(X = X.1, Y= Y, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
  mutate(dist4 = nndist(X = X.1, Y= Y, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
  mutate(dist5 = nndist(X = X.1, Y= Y, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
  rowwise()%>% #so that in the next part we take the averages across rows
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5)))  %>% #creates a column of the average distances (1-5) of each individual
  dplyr::select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances


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



## Load in environmental rasters ##


#loading in soil textures from CONABIO, theses are too larger, about 1 km^2 I believe
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
nitrogen_05 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 0-5.tif"))
nitrogen_200 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 100-200.tif"))
Soil_Organic_Carbon_05 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 0-5.tif"))
Soil_Organic_Carbon_200 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 100-200.tif"))


#project rasters to equal area projection (UTM 12N), uses meters as distance measurement 
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
sandy_avail_water_0.5_utm <- vol_wat_33kpa_05_utm - vol_wat_1500kpa_05_utm
sandy_avail_water_100.200_utm <- vol_wat_33kpa_200_utm - vol_wat_1500kpa_200_utm
clay_loam_avail_water_0.5_utm <- vol_wat_10kpa_05_utm - vol_wat_1500kpa_05_utm
clay_loam_avail_water_100.200_utm <- vol_wat_10kpa_200_utm - vol_wat_1500kpa_200_utm


#LM
#examining the layers at different extents

#using the extent of the box around the rivers to crop the raster for each soil texture layer
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


#attempt of using ggplot to plot clay layer with river shapefile
ggplot()+
  geom_raster(data = as.data.frame(Soil_Organic_Carbon_05_LM, xy=T), aes(x=x, y=y, fill = SOC.0.5))+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed)


ggplot()+
  geom_raster(data = as.data.frame(ph_200_LM, xy=T), aes(x=x, y=y, fill = ph_100.200))+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed)

#creating a stack of the raster layers
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
#using the extent of the box around the rivers to crop the raster for each soil texture layer
#using the extent of the box around the rivers to crop the raster for each soil texture layer
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


#creating a stack of the raster layers 
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


#creating a stack of the raster layers
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


#FOR A PRESENTATION
library(gridExtra)
nitrogen0.5 <- ggplot() +
  geom_raster(data= as.data.frame(soil_stack_SD_extra, xy = T), aes(x=x, y=y, fill = nitrogen.0.5))+
  geom_sf(data = SD_fixed_field_data_processed, fill = NA, col = "white") +
  labs(title = "Nitrogen at 0-5 cm",
       fill = "Nitrogen (cg/kg)",
       x = "",
       y = "") +
  scale_fill_viridis_c(limits = c(50, 330))

nitrogen100.200 <- ggplot() +
  geom_raster(data= as.data.frame(soil_stack_SD_extra, xy = T), aes(x=x, y=y, fill = nitrogen.100.200))+
  geom_sf(data = SD_fixed_field_data_processed, fill = NA, col = "black") +
  labs(title = "Nitrogen at 100-200 cm", 
       fill = "Nitrogen (cg/kg)",
       x = "",
       y = "") +
  scale_fill_viridis_c(limits = c(50, 330))

grid.arrange(nitrogen0.5, nitrogen100.200, nrow = 2)


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

#Extracting the soil data to the tree points 

#LM
LM_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LM_soil_text, LM_fixed_field_data_processed) #extracting soil textures for each point value
LM_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LM_other, LM_fixed_field_data_processed) #extracting the other soil variables for each point value
LM_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LM_extra, LM_fixed_field_data_processed) #extracting the extra soil variables for each point value
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed, LM_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LM point dataframe
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LM point dataframe
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LM point dataframe


#LC
LC_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LC_soil_text, LC_fixed_field_data_processed) #extracting soil textures for each point value
LC_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LC_other, LC_fixed_field_data_processed) #extracting the other soil variables for each point value
LC_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LC_extra, LC_fixed_field_data_processed) #extracting the extra soil variables for each point value
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed, LC_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe


#SD
SD_soil_text_raster_250_data_pts <- raster::extract(soil_stack_SD_soil_text, SD_fixed_field_data_processed) #extracting soil textures for each point value
SD_soil_other_raster_250_data_pts <- raster::extract(soil_stack_SD_other, SD_fixed_field_data_processed) #extracting the other soil variables for each point value
SD_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_SD_extra, SD_fixed_field_data_processed) #extracting the extra soil variables for each point value
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed, SD_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe


#### Creating Sandy and Clay/Loamy Available Water Columns ####

# LM

LM_fixed_field_data_processed_soils <- LM_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(sandy_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) %>%
  mutate(clay_loam_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(clay_loam_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) 


# LC

LC_fixed_field_data_processed_soils <- LC_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(sandy_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) %>%
  mutate(clay_loam_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(clay_loam_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200)

# SD

SD_fixed_field_data_processed_soils <- SD_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(sandy_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) %>%
  mutate(clay_loam_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(clay_loam_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200)

#combining all of the population soil data into one dataframe 
fixed_field_data_processed_soils <- rbind(LM_fixed_field_data_processed_soils, LC_fixed_field_data_processed_soils)
fixed_field_data_processed_soils <- rbind(fixed_field_data_processed_soils, SD_fixed_field_data_processed_soils)



#### PART 3: Comparing average soil values from inside populations to outside populations ####

#downloading the locations of the 20 known populations
all_pop_locations.df <- read.csv(file = "./data/Known QUBR populations.xlsx - More accurate GPD coords for pops (12_2024).csv")
View(all_pop_locations.df)


#transforming the data into shapefiles with either WGS84 
all_pop_locations.df_sf <- st_as_sf(all_pop_locations.df, 
                                          coords = c("Classic.Longitude", "Classic.Latitude"), crs = 4326)

#transforming the shapefile of trees from WGS84 into equal area projection UTM 12N
all_pop_locations.df_sf_transformed <- st_transform(all_pop_locations.df_sf, crs = 26912) # this in UTM 12 N an equal area projection

#create dataframe with X and Y UTM coordinates
all_pop_locations.df_sf_trans_coords <- st_coordinates(all_pop_locations.df_sf_transformed) #creates a dataframe with seperate x and y columns from the UTM 12N transformation
all_pop_locations.df_sf_trans_coordinates <- all_pop_locations.df_sf_transformed %>%
  cbind(all_pop_locations.df_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

# importing the baja sur shapefile in

#turn the BCS polygon into a shapefile and visualize its outline
BCS_polygon <- read_sf("./data/Shapefiles/BCS_Shapefile/bcs_entidad.shp")
BCS_polygon <- st_as_sf(BCS_polygon)
BCS_polygon_UTM <- st_transform(BCS_polygon, crs = 26912) # this in UTM 12 N an equal area projection
BCS_polygon_UTM <- st_as_sf(BCS_polygon_UTM)


#cropping the BCS polygon to just be the southern region of where the 20 known populations are with a 10 km radiu

#BCS_polygon_box <- st_bbox(BCS_polygon_UTM)


# Get bbox of points
bbox_points <- st_bbox(all_pop_locations.df_sf_trans_coordinates)

# Convert bbox to polygon
bbox_poly <- st_as_sfc(bbox_points)

# Buffer polygon by 7,000 meters (10 km)
bbox_poly_buffered <- st_buffer(bbox_poly, dist = 7000)

# Make sure CRS matches your big polygon
bbox_poly_buffered <- st_transform(bbox_poly_buffered, st_crs(BCS_polygon_UTM))

# Crop (intersect) the big polygon with buffered bbox polygon
BCS_polygon_box_sf_cropped <- st_intersection(BCS_polygon_UTM, bbox_poly_buffered)

# Check bbox of cropped polygon (should be larger than original bbox of points)
print(st_bbox(BCS_polygon_box_sf_cropped))

#plotting the original extent and polygon and the cropped extent
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=bbox_poly_buffered)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)

#plotting just the cropped area
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)



#using the extent of the box around the rivers to crop the raster for each soil texture layer
#clay 05
clay_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(clay_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
clay_05_clipped <- crop(clay_05_utm, clay_05_bbox_poly_buffered) # Crop the raster to the polygon
clay_05_all_pop <- mask(clay_05_clipped, clay_05_bbox_poly_buffered) # Mask the raster to the polygon
#clay 200
clay_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(clay_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
clay_200_clipped <- crop(clay_200_utm, clay_200_bbox_poly_buffered) # Crop the raster to the polygon
clay_200_all_pop <- mask(clay_200_clipped, clay_200_bbox_poly_buffered) # Mask the raster to the polygon
#silt 05
silt_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(silt_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
silt_05_clipped <- crop(silt_05_utm, silt_05_bbox_poly_buffered) # Crop the raster to the polygon
silt_05_all_pop <- mask(silt_05_clipped, silt_05_bbox_poly_buffered) # Mask the raster to the polygon
#silt 200
silt_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(silt_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
silt_200_clipped <- crop(silt_200_utm, silt_200_bbox_poly_buffered) # Crop the raster to the polygon
silt_200_all_pop <- mask(silt_200_clipped, silt_200_bbox_poly_buffered) # Mask the raster to the polygon
#sand 05
sand_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(sand_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
sand_05_clipped <- crop(sand_05_utm, sand_05_bbox_poly_buffered) # Crop the raster to the polygon
sand_05_all_pop <- mask(sand_05_clipped, sand_05_bbox_poly_buffered) # Mask the raster to the polygon
#sand 200
sand_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(sand_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
sand_200_clipped <- crop(sand_200_utm, sand_200_bbox_poly_buffered) # Crop the raster to the polygon
sand_200_all_pop <- mask(sand_200_clipped, sand_200_bbox_poly_buffered) # Mask the raster to the polygon
#ph 05
ph_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ph_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ph_05_clipped <- crop(ph_05_utm, ph_05_bbox_poly_buffered) # Crop the raster to the polygon
ph_05_all_pop <- mask(ph_05_clipped, ph_05_bbox_poly_buffered) # Mask the raster to the polygon
#ph 200
ph_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ph_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ph_200_clipped <- crop(ph_200_utm, ph_200_bbox_poly_buffered) # Crop the raster to the polygon
ph_200_all_pop <- mask(ph_200_clipped, ph_200_bbox_poly_buffered) # Mask the raster to the polygon
#ph 05
ph_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ph_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ph_05_clipped <- crop(ph_05_utm, ph_05_bbox_poly_buffered) # Crop the raster to the polygon
ph_05_all_pop <- mask(ph_05_clipped, ph_05_bbox_poly_buffered) # Mask the raster to the polygon
#ph 200
ph_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ph_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ph_200_clipped <- crop(ph_200_utm, ph_200_bbox_poly_buffered) # Crop the raster to the polygon
ph_200_all_pop <- mask(ph_200_clipped, ph_200_bbox_poly_buffered) # Mask the raster to the polygon
#ocd 05
ocd_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ocd_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ocd_05_clipped <- crop(ocd_05_utm, ocd_05_bbox_poly_buffered) # Crop the raster to the polygon
ocd_05_all_pop <- mask(ocd_05_clipped, ocd_05_bbox_poly_buffered) # Mask the raster to the polygon
#ocd 200
ocd_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(ocd_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
ocd_200_clipped <- crop(ocd_200_utm, ocd_200_bbox_poly_buffered) # Crop the raster to the polygon
ocd_200_all_pop <- mask(ocd_200_clipped, ocd_200_bbox_poly_buffered) # Mask the raster to the polygon
#coarse_frag 05
coarse_frag_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(coarse_frag_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
coarse_frag_05_clipped <- crop(coarse_frag_05_utm, coarse_frag_05_bbox_poly_buffered) # Crop the raster to the polygon
coarse_frag_05_all_pop <- mask(coarse_frag_05_clipped, coarse_frag_05_bbox_poly_buffered) # Mask the raster to the polygon
#coarse_frag 200
coarse_frag_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(coarse_frag_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
coarse_frag_200_clipped <- crop(coarse_frag_200_utm, coarse_frag_200_bbox_poly_buffered) # Crop the raster to the polygon
coarse_frag_200_all_pop <- mask(coarse_frag_200_clipped, coarse_frag_200_bbox_poly_buffered) # Mask the raster to the polygon
#cat_ex_cap 05
cat_ex_cap_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(cat_ex_cap_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
cat_ex_cap_05_clipped <- crop(cat_ex_cap_05_utm, cat_ex_cap_05_bbox_poly_buffered) # Crop the raster to the polygon
cat_ex_cap_05_all_pop <- mask(cat_ex_cap_05_clipped, cat_ex_cap_05_bbox_poly_buffered) # Mask the raster to the polygon
#cat_ex_cap 200
cat_ex_cap_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(cat_ex_cap_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
cat_ex_cap_200_clipped <- crop(cat_ex_cap_200_utm, cat_ex_cap_200_bbox_poly_buffered) # Crop the raster to the polygon
cat_ex_cap_200_all_pop <- mask(cat_ex_cap_200_clipped, cat_ex_cap_200_bbox_poly_buffered) # Mask the raster to the polygon
#bulk_dens 05
bulk_dens_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(bulk_dens_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
bulk_dens_05_clipped <- crop(bulk_dens_05_utm, bulk_dens_05_bbox_poly_buffered) # Crop the raster to the polygon
bulk_dens_05_all_pop <- mask(bulk_dens_05_clipped, bulk_dens_05_bbox_poly_buffered) # Mask the raster to the polygon
#bulk_dens 200
bulk_dens_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(bulk_dens_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
bulk_dens_200_clipped <- crop(bulk_dens_200_utm, bulk_dens_200_bbox_poly_buffered) # Crop the raster to the polygon
bulk_dens_200_all_pop <- mask(bulk_dens_200_clipped, bulk_dens_200_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_10kpa 05
vol_wat_10kpa_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_10kpa_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_10kpa_05_clipped <- crop(vol_wat_10kpa_05_utm, vol_wat_10kpa_05_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_10kpa_05_all_pop <- mask(vol_wat_10kpa_05_clipped, vol_wat_10kpa_05_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_10kpa 200
vol_wat_10kpa_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_10kpa_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_10kpa_200_clipped <- crop(vol_wat_10kpa_200_utm, vol_wat_10kpa_200_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_10kpa_200_all_pop <- mask(vol_wat_10kpa_200_clipped, vol_wat_10kpa_200_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_33kpa 05
vol_wat_33kpa_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_33kpa_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_33kpa_05_clipped <- crop(vol_wat_33kpa_05_utm, vol_wat_33kpa_05_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_33kpa_05_all_pop <- mask(vol_wat_33kpa_05_clipped, vol_wat_33kpa_05_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_33kpa 200
vol_wat_33kpa_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_33kpa_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_33kpa_200_clipped <- crop(vol_wat_33kpa_200_utm, vol_wat_33kpa_200_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_33kpa_200_all_pop <- mask(vol_wat_33kpa_200_clipped, vol_wat_33kpa_200_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_1500kpa 05
vol_wat_1500kpa_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_1500kpa_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_1500kpa_05_clipped <- crop(vol_wat_1500kpa_05_utm, vol_wat_1500kpa_05_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_1500kpa_05_all_pop <- mask(vol_wat_1500kpa_05_clipped, vol_wat_1500kpa_05_bbox_poly_buffered) # Mask the raster to the polygon
#vol_wat_1500kpa 200
vol_wat_1500kpa_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(vol_wat_1500kpa_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
vol_wat_1500kpa_200_clipped <- crop(vol_wat_1500kpa_200_utm, vol_wat_1500kpa_200_bbox_poly_buffered) # Crop the raster to the polygon
vol_wat_1500kpa_200_all_pop <- mask(vol_wat_1500kpa_200_clipped, vol_wat_1500kpa_200_bbox_poly_buffered) # Mask the raster to the polygon
#nitrogen 05
nitrogen_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(nitrogen_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
nitrogen_05_clipped <- crop(nitrogen_05_utm, nitrogen_05_bbox_poly_buffered) # Crop the raster to the polygon
nitrogen_05_all_pop <- mask(nitrogen_05_clipped, nitrogen_05_bbox_poly_buffered) # Mask the raster to the polygon
#nitrogen 200
nitrogen_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(nitrogen_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
nitrogen_200_clipped <- crop(nitrogen_200_utm, nitrogen_200_bbox_poly_buffered) # Crop the raster to the polygon
nitrogen_200_all_pop <- mask(nitrogen_200_clipped, nitrogen_200_bbox_poly_buffered) # Mask the raster to the polygon
#Soil_Organic_Carbon 05
Soil_Organic_Carbon_05_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(Soil_Organic_Carbon_05_utm)) # Make sure your sfc polygon is in the same CRS as the raster
Soil_Organic_Carbon_05_clipped <- crop(Soil_Organic_Carbon_05_utm, Soil_Organic_Carbon_05_bbox_poly_buffered) # Crop the raster to the polygon
Soil_Organic_Carbon_05_all_pop <- mask(Soil_Organic_Carbon_05_clipped, Soil_Organic_Carbon_05_bbox_poly_buffered) # Mask the raster to the polygon
#Soil_Organic_Carbon 200
Soil_Organic_Carbon_200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(Soil_Organic_Carbon_200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
Soil_Organic_Carbon_200_clipped <- crop(Soil_Organic_Carbon_200_utm, Soil_Organic_Carbon_200_bbox_poly_buffered) # Crop the raster to the polygon
Soil_Organic_Carbon_200_all_pop <- mask(Soil_Organic_Carbon_200_clipped, Soil_Organic_Carbon_200_bbox_poly_buffered) # Mask the raster to the polygon
#sandy_avail_water_0.5
sandy_avail_water_0.5_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(sandy_avail_water_0.5_utm)) # Make sure your sfc polygon is in the same CRS as the raster
sandy_avail_water_0.5_clipped <- crop(sandy_avail_water_0.5_utm, sandy_avail_water_0.5_bbox_poly_buffered) # Crop the raster to the polygon
sandy_avail_water_0.5_all_pop <- mask(sandy_avail_water_0.5_clipped, sandy_avail_water_0.5_bbox_poly_buffered) # Mask the raster to the polygon
#sandy_avail_water_100.200
sandy_avail_water_100.200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(sandy_avail_water_100.200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
sandy_avail_water_100.200_clipped <- crop(sandy_avail_water_100.200_utm, sandy_avail_water_100.200_bbox_poly_buffered) # Crop the raster to the polygon
sandy_avail_water_100.200_all_pop <- mask(sandy_avail_water_100.200_clipped, sandy_avail_water_100.200_bbox_poly_buffered) # Mask the raster to the polygon
#clay_loam_avail_water_0.5
clay_loam_avail_water_0.5_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(clay_loam_avail_water_0.5_utm)) # Make sure your sfc polygon is in the same CRS as the raster
clay_loam_avail_water_0.5_clipped <- crop(clay_loam_avail_water_0.5_utm, clay_loam_avail_water_0.5_bbox_poly_buffered) # Crop the raster to the polygon
clay_loam_avail_water_0.5_all_pop <- mask(clay_loam_avail_water_0.5_clipped, clay_loam_avail_water_0.5_bbox_poly_buffered) # Mask the raster to the polygon
#clay_loam_avail_water_100.200
clay_loam_avail_water_100.200_bbox_poly_buffered <- st_transform(BCS_polygon_box_sf_cropped, crs(clay_loam_avail_water_100.200_utm)) # Make sure your sfc polygon is in the same CRS as the raster
clay_loam_avail_water_100.200_clipped <- crop(clay_loam_avail_water_100.200_utm, clay_loam_avail_water_100.200_bbox_poly_buffered) # Crop the raster to the polygon
clay_loam_avail_water_100.200_all_pop <- mask(clay_loam_avail_water_100.200_clipped, clay_loam_avail_water_100.200_bbox_poly_buffered) # Mask the raster to the polygon



#confirming I properly cropped the rasters by plotting the clay rasters with the cropped polygon around it
ggplot()+
  geom_raster(data = as.data.frame(clay_05_utm, xy=T), aes(x=x, y=y, fill = clay.content.0.5)) +
  scale_fill_gradientn(colours=c("yellow","red"), name = "clay_05_utm")+
  # Add new fill scale
  ggnewscale::new_scale_fill() +
  geom_raster(data = as.data.frame(clay_05_all_pop, xy=T), aes(x=x, y=y, fill = clay.content.0.5)) +
  scale_fill_gradientn(colours=c("lightblue","darkblue"), name = "clay_05_all_pop") +
  geom_sf(data = BCS_polygon_box_sf_cropped, fill = NA, color = "green")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)

ggplot() +
  #geom_sf(data=BCS_polygon_UTM)+
  #geom_sf(data=bbox_poly_buffered)+
  geom_raster(data= as.data.frame(vol_wat_33kpa_05_all_pop, xy = T), aes(x=x, y=y, fill = vol_water_0.5))+
  geom_sf(data = BCS_polygon_box_sf_cropped, fill = NA, color = "green")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates) +
  geom_sf(data=random_20, color ="red") + 
  labs(title = "Volume of Water Content at 0-5 cm",
       fill = "Water Content (10-2 cm3 cm-3 )*10",
       x = "",
       y = "") +
  scale_fill_viridis_c(limits = c(50, 330)) +
  theme(legend.title =  element_text(size = 10))


#creating a stack of the raster layers for the original rasters

soil_stack_clay <- stack(clay_05_all_pop, clay_200_all_pop) #stacked clay
soil_stack_silt <- stack(silt_05_all_pop, silt_200_all_pop) #stacked silt
soil_stack_sand <- stack(sand_05_all_pop, sand_200_all_pop) #stacked sand
soil_stack_ph <- stack(ph_05_all_pop, ph_200_all_pop) #stacked ph
soil_stack_ocd <- stack(ocd_05_all_pop, ocd_200_all_pop)  #stacked ocd
soil_stack_coarse_frag <- stack(coarse_frag_05_all_pop, coarse_frag_200_all_pop)  #stacked coarse fragment
soil_stack_cat_ex <- stack(cat_ex_cap_05_all_pop, cat_ex_cap_200_all_pop)  #stacked cation exchange capacity
soil_stack_bulk_dens <- stack(bulk_dens_05_all_pop, bulk_dens_200_all_pop)  #stacked bulk density
soil_stack_vol_wat_10kpa <- stack(vol_wat_10kpa_05_all_pop, vol_wat_10kpa_200_all_pop)  #stacked volume water content at 10 kpa
soil_stack_vol_wat_33kpa <- stack(vol_wat_33kpa_05_all_pop, vol_wat_33kpa_200_all_pop)  #stacked volume water content at 33 kpa
soil_stack_vol_wat_1500kpa <- stack(vol_wat_1500kpa_05_all_pop, vol_wat_1500kpa_200_all_pop) #stacked volume water content at 1500 kpa
soil_stack_nitrogen <- stack(nitrogen_05_all_pop, nitrogen_200_all_pop)  #stacked volume water content at 10 kpa
soil_stack_soc <- stack(Soil_Organic_Carbon_05_all_pop, Soil_Organic_Carbon_200_all_pop)  #stacked volume water content at 10 kpa
soil_stack_sandy_water <- stack(sandy_avail_water_0.5_all_pop, sandy_avail_water_100.200_all_pop) #stacked available water in sandy soils
soil_stack_clay_loam_water <- stack(clay_loam_avail_water_0.5_all_pop, clay_loam_avail_water_100.200_all_pop) #stacked available water in clay and loam soils

#plotting the stacked rasters, example with clay
plot(soil_stack_clay) #version with soil textures
plot(soil_stack_clay, zlim = c(0, 350)) #version where the plots have the same scale


#extracting the soil data for each point of the known 20 points 
all_known_pop_soil_clay <- raster::extract(soil_stack_clay, all_pop_locations.df_sf_trans_coordinates) #extracting soil textures for each point value
all_known_pop_soil_silt <- raster::extract(soil_stack_silt, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_sand <- raster::extract(soil_stack_sand, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_ph <- raster::extract(soil_stack_ph, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_ocd <- raster::extract(soil_stack_ocd, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_coarse_frag <- raster::extract(soil_stack_coarse_frag, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_cat_ex <- raster::extract(soil_stack_cat_ex, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_bulk_dens <- raster::extract(soil_stack_bulk_dens, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_vol_wat_10kpa <- raster::extract(soil_stack_vol_wat_10kpa, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_vol_wat_33kpa <- raster::extract(soil_stack_vol_wat_33kpa, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_vol_wat_1500kpa <- raster::extract(soil_stack_vol_wat_1500kpa, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_nitrogen <- raster::extract(soil_stack_nitrogen, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_soc <- raster::extract(soil_stack_soc, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_sandy_water <- raster::extract(soil_stack_sandy_water, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_clay_loam_water <- raster::extract(soil_stack_clay_loam_water, all_pop_locations.df_sf_trans_coordinates)

#bind the soil textures data for each point to the all point dataframe so the soil values are available for each population point
all_known_pop_soils <- cbind(all_pop_locations.df_sf_trans_coordinates, all_known_pop_soil_clay) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_silt) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_sand) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_ph) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_ocd)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_cat_ex)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_bulk_dens)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_vol_wat_10kpa)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_vol_wat_33kpa)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_vol_wat_1500kpa)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_nitrogen)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_soc)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_sandy_water)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_clay_loam_water)

#correcting column names for clarity
all_known_pop_soils <- all_known_pop_soils %>%
  mutate(sandy_avail_water_0.5 = layer.1) %>%
  mutate(sandy_avail_water_100.200 = layer.2) %>% 
  mutate(clay_loam_avail_water_0.5 = layer.1.1) %>%
  mutate(clay_loam_avail_water_100.200 = layer.2.1) 

#creating the list of soil metrics to iterate over
Soil.metrics <- c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                  "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                  "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                  "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                  "Volume of water content -1500 kpa 100-200", 
                  "Nitrogen 0-5", "Nitrogen 100-200", 
                  "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                  "Sand Available Water 0-5", "Sand Available Water 100-200",
                  "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200")


random_pop_soils <- function(){
  

  #creating empty list to collect p_values
  random_soil_p_values <- c()  #for the p values of the randomly generated population means compared to our known soil mean
  known_soil_means <- c()  #for the p values of the known population means
  
  
  plot_list <- list()   #to store plots
  
  for (i in 1:length(Soil.metrics)){
    
    #assigning the population based on the current tree
    if (Soil.metrics[i] == "Clay 0-5"){ 
      soil_stack = soil_stack_clay
      soil_metric = all_known_pop_soils$clay.content.0.5
    } else if (Soil.metrics[i] == "Clay 100-200"){
      soil_stack = soil_stack_clay
      soil_metric = all_known_pop_soils$clay.content.100.200
    } else if (Soil.metrics[i] == "Silt 0-5"){
      soil_stack = soil_stack_silt
      soil_metric = all_known_pop_soils$silt.0.5
    } else if (Soil.metrics[i] == "Silt 100-200"){
      soil_stack = soil_stack_silt
      soil_metric = all_known_pop_soils$silt.100.200
    } else if (Soil.metrics[i] == "Sand 0-5"){
      soil_stack = soil_stack_sand
      soil_metric = all_known_pop_soils$sand.0.5
    } else if (Soil.metrics[i] == "Sand 100-200"){
      soil_stack = soil_stack_sand
      soil_metri = all_known_pop_soils$sand.100.200
    } else if (Soil.metrics[i] == "Ph 0-5"){
      soil_stack = soil_stack_ph
      soil_metric = all_known_pop_soils$ph_0.5
    } else if (Soil.metrics[i] == "Ph 100-200"){
      soil_stack = soil_stack_ph
      soil_metric = all_known_pop_soils$ph_100.200
    } else if (Soil.metrics[i] == "Volume of water content -10 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_10kpa
      soil_metric = all_known_pop_soils$vol_water_.10_0.5
    } else if (Soil.metrics[i] == "Volume of water content -10 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_10kpa
      soil_metric = all_known_pop_soils$vol_water_.10_100.200
    } else if (Soil.metrics[i] == "Volume of water content -33 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_33kpa
      soil_metric = all_known_pop_soils$vol_water_0.5
    } else if (Soil.metrics[i] == "Volume of water content -33 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_33kpa
      soil_metric = all_known_pop_soils$vol_water_100.200
    } else if (Soil.metrics[i] == "Volume of water content -1500 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_1500kpa
      soil_metric = all_known_pop_soils$vol_water_.1500kPa_0.5
    } else if (Soil.metrics[i] == "Volume of water content -1500 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_1500kpa
      soil_metric = all_known_pop_soils$vol_water_.1500_100.200
    } else if (Soil.metrics[i] == "Nitrogen 0-5"){
      soil_stack = soil_stack_nitrogen
      soil_metric = all_known_pop_soils$nitrogen.0.5
    } else if (Soil.metrics[i] == "Nitrogen 100-200"){
      soil_stack = soil_stack_nitrogen
      soil_metric = all_known_pop_soils$nitrogen.100.200
    } else if (Soil.metrics[i] == "Soil Organic Carbon 0-5"){
      soil_stack = soil_stack_soc
      soil_metric = all_known_pop_soils$SOC.0.5
    } else if (Soil.metrics[i] == "Soil Organic Carbon 100-200"){
      soil_stack = soil_stack_soc
      soil_metric = all_known_pop_soils$SOC.100.200
    } else if (Soil.metrics[i] == "Sand Available Water 0-5"){
      soil_stack = soil_stack_sandy_water
      soil_metric = all_known_pop_soils$sandy_avail_water_0.5
    } else if (Soil.metrics[i] == "Sand Available Water 100-200"){
      soil_stack = soil_stack_sandy_water
      soil_metric = all_known_pop_soils$sandy_avail_water_100.200
    } else if (Soil.metrics[i] == "Clay/Loam Available Water 0-5"){
      soil_stack = soil_stack_clay_loam_water
      soil_metric = all_known_pop_soils$clay_loam_avail_water_0.5
    } else if (Soil.metrics[i] == "Clay/Loam Available Water 100-200"){
      soil_stack = soil_stack_clay_loam_water
      soil_metric = all_known_pop_soils$clay_loam_avail_water_100.200
    } 
    
    #creating empty list to collect means
    random_soil_means <- c()  #for the means of the randomly generated population means
     
    set.seed(20)
    for (y in 1:1000){ #for 1000 permutations
      
      random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
      random_20 <- random_20 %>%
        st_as_sf()
      random_20_pop_soil <- raster::extract(soil_stack, random_20) #extracting the soil metrics for the random points
      
      #storing the mean of the soil metric depending on if it is the 0-5 or 100-200 cm version
      if (i %% 2 == 1){ #if the iteration we are on is odd, then we use the 0-5 cm variable
        random_mean <- mean(random_20_pop_soil[,1]) #storing the mean of the 0-5 value
      } else {  #if the iteration we are on is odd, then we use the 100-200 cm variable
        random_mean <- mean(random_20_pop_soil[,2]) #storing the mean of the 100-200 value
      }
      
      random_soil_means <- c(random_soil_means, random_mean) #adding the 0-5 or 100-200 cm mean to the list of means
      
    }
    
    #plotting the randomly selected points on the Baja polygon
    random_points_BCS <- ggplot()+
      geom_sf(data=BCS_polygon_UTM)+
      geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
      geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
      geom_sf(data=random_20, color ="blue")
    
    #plotting the randomly selected points just on the cropped polygon
    random_points_BCS_crop <- ggplot()+
      geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
      geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
      geom_sf(data=random_20, color ="blue")
  
    #storing the real means
    all_known_mean <- mean(soil_metric)
    
    #adding the known population mean soil metric to the list
    known_soil_means <- c(known_soil_means, all_known_mean)
    
    #plotting the histogram of the randomly distributed p-values and our real slope
    plot_out_histogram <- ggplot()+
      geom_histogram(aes(x=random_soil_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
      geom_vline(xintercept=all_known_mean, col = "red")+ #line of our real slope
      xlab(paste0("Mean ", Soil.metrics[i], " of Random Populations vs.Known Populations (n=20)"))+
      theme_classic()
    
    # store the histogram in list with a descriptive name
    plot_name_histogram <- paste(Soil.metrics[i], "Histogram",
                       sep = "_")
    plot_list[[plot_name_histogram]] <- plot_out_histogram
    
    random_soil_means <- na.omit(random_soil_means) #removing NAs
    
    # if using greater than hypothesis
    
    #calculating pseudo p-value for 
    # total = 0  #set empty value
    # for (k in 1:length(random_soil_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
    #   if (random_soil_means[k] < all_known_mean){
    #     total = total + 1
    #   }
    # } #add number of values of in the random set of means values that are less than our mean ANN
    # random_p.value <- 1 - (total / length(random_soil_means)) #the proportion of random ANNs that are greater than our ANN
    
    # using the significantly different alternative hypothesis 
    
    p_value_greater_than <- sum(random_soil_means >= all_known_mean)/length(random_soil_means)   # proportion of simulated slopes higher than our real slope
    p_value_less_than <- sum(random_soil_means <= all_known_mean)/length(random_soil_means)   # proportion of simulated slopes lower than our real slope
    random_p.value <- min(1, 2 * min(p_value_greater_than, p_value_less_than)) # take the smaller tail (the "more extreme" one), then double it
    
    
    #adding the p value to total list of p-values for all soil metrics
    random_soil_p_values <- c(random_soil_p_values, random_p.value)
    
    #print(paste("Updating:", i))
    
  }
  
  return(list(known_soil_means = known_soil_means,
         random_soil_p_values = random_soil_p_values, 
         random_points_BCS = random_points_BCS, 
         random_points_BCS_crop = random_points_BCS_crop,
         plot_list = plot_list))

}

#running and storing the function and its results
random_pop_soils_function <- random_pop_soils()

#Example of extracting one of the histograms comparing the slopes for our original soil vs. size metrics to the shuffled ones
random_pop_soils_function$plot_list$`Clay 0-5_Histogram`
plot <- random_pop_soils_function$plot_list$`Clay 0-5_Histogram`

#if you want to see all of the plots at once run: 
#random_pop_soils_function$plot_list

# Bonferroni correcting for multiple testing
p_bonf_corrected <- p.adjust(random_pop_soils_function$random_soil_p_values, method = "bonferroni")
p_bonf_corrected
  
#making a dataframe from the function output
random_pop.df <- data.frame("Soil.metrics" = Soil.metrics, 
                            "P_values" = random_pop_soils_function$random_soil_p_values, 
                            "P_values_bonf_corrected" = p_bonf_corrected,
                            "Significance" = c(rep(NA, 22)))

#creating the significance column for the p-values
random_pop.df <- random_pop.df %>%
  mutate(Significance = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                  p_bonf_corrected >= 0.05 ~ "N"))


#Heat Map 


#labeled p-values
ggplot(aes(x = fct_reorder(Soil.metrics, P_values), y = Significance, fill = P_values), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Association Between Soil Metrics and Population Locations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_values < 0.05, round(P_values, 4), NA)), col = "white") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))

