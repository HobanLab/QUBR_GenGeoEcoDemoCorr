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


### Comparing the soil metrics between populations ###

#LM

#creating a grid over the soil cells
LM_tree_grid_cropped <- st_make_grid(soil_stack_LM_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_LM_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = LM_tree_grid_cropped, fill = NA)

#selecting a point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, LM_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
LM_list_grids_and_trees <- lapply(LM_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})

#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
LM_list_grids_and_point_trees_df <- as.data.frame(unlist(LM_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LM_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
LM_list_grids_and_trees_fixed <- LM_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
LM_fixed_field_data_processed_trees_soils <- LM_fixed_field_data_processed_soils %>%
  filter(X %in% LM_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data= LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_fixed_field_data_processed_trees_soils, color = "red")


#LC

#creating a grid over the soil cells
LC_tree_grid_cropped <- st_make_grid(soil_stack_LC_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_LC_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = LC_tree_grid_cropped, fill = NA)

####NOTE FROM ASH:I DONT THINK THIS CODE WORKS RIGHT NOW CORRECTLY BC THE ST_CONTAINS NEEDS TO BE THE CROPPED RASTER AND THEN ALL OF THE TREE DATA FOR THE ROW CALLS TO WORK CORRECTLY#####
#selecting a point from each grid cell with trees within them
LC_list_grids_and_points <- st_contains(LC_tree_grid_cropped, LC_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
LC_list_grids_and_trees <- lapply(LC_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})


#creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
LC_list_grids_and_point_trees_df <- as.data.frame(unlist(LC_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LC_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
LC_list_grids_and_trees_fixed <- LC_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
  mutate(data_row = LC_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them


#filtering out point data to be just the focal points
LC_fixed_field_data_processed_trees_soils <- LC_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% LC_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data= LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_fixed_field_data_processed_trees_soils, color = "red")


#SD

#creating a grid over the soil cells
SD_tree_grid_cropped <- st_make_grid(soil_stack_SD_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_SD_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = SD_tree_grid_cropped, fill = NA)


#selecting a point from each grid cell with trees within them
SD_list_grids_and_points <- st_contains(SD_tree_grid_cropped, SD_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
SD_list_grids_and_trees <- lapply(SD_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})


#creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
SD_list_grids_and_point_trees_df <- as.data.frame(unlist(SD_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(SD_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
SD_list_grids_and_trees_fixed <- SD_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
  mutate(data_row = SD_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them


#filtering out point data to be just the focal points
SD_fixed_field_data_processed_trees_soils <- SD_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% SD_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = SD_tree_grid_cropped)+
  geom_sf(data= SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_fixed_field_data_processed_trees_soils, color = "red")


ggplot()+
  geom_raster(data = as.data.frame(clay_05_SD, xy = T), aes(x=x, y=y, fill = clay.content.0.5)) +
 # geom_sf(data = st_boundary(SD_tree_grid_cropped)) + 
  labs(fill = "Clay Content (g/kg) (0-5 cm)") +
  geom_sf(data= SD_fixed_field_data_processed_sf)+
 # geom_sf(data = SD_fixed_field_data_processed_trees_soils, color = "red") +
  scale_fill_viridis_c()


#combining the LM, LC, and SD tree randomly chosen tree point data into one dataframe

fixed_field_data_processed_trees_soils <- rbind(LM_fixed_field_data_processed_trees_soils, LC_fixed_field_data_processed_trees_soils) #combining the LM and LC soil and randomly chosen tree data
fixed_field_data_processed_trees_soils <- rbind(fixed_field_data_processed_trees_soils, SD_fixed_field_data_processed_trees_soils) #combining the SD tree point data to the LM and LC soil and randomly chosen tree point data


#creating a locality as factor column to be able to use Tamhane's T2 Test later
fixed_field_data_processed_trees_soils$Locality_Factor <- as.factor(fixed_field_data_processed_trees_soils$Locality)

#### Comparing the soil vs. size values ####


#FOR THIS WE ARE GOING TO SHUFFLE THE SHAPE/SIZE VALUES BETWEEN ALL OF THE POINTS IN 
#EACH PERMUTATION AND CALCULATE THE SLOPE OF THE SHAPE AND SOIL VALUE, WE THEN COMPARE THE RANDOMIZED SLOPE 
#TO OUR REAL SLOPE AND EXTRACT THE P-VALUE

slopes_simulations <- function () {
  
  Populations <- c("LM", "LC", "SD")
  Size.metrics <- c("SCA", "LCA", "CA", "CS", "DBH")
  Soil.metrics <- c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                     "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                     "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                    "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                    "Volume of water content -1500 kpa 100-200", 
                    "Nitrogen 0-5", "Nitrogen 100-200", 
                    "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                    "Sand Available Water 0-5", "Sand Available Water 100-200",
                    "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200")
  slopes_array <- array(NA, dim = c(22, 5, 3), 
                        list(Soil.metrics, Size.metrics, Populations))
  p_values_array <- array(NA, dim = c(22, 5, 3), 
                        list(Soil.metrics, Size.metrics, Populations))

  
  for (i in 1:length(Populations)) { #iterating over the total tree dataframe
   
     #assigning the population based on the current tree
     if (Populations[i] == "LM"){ 
       Population = "LM"
       dataframe_soils = LM_fixed_field_data_processed_soils
     } else if (Populations[i] == "LC"){
       Population = "LC"
       dataframe_soils = LC_fixed_field_data_processed_soils
     } else if (Populations[i] == "SD"){
       Population = "SD"
       dataframe_soils = SD_fixed_field_data_processed_soils
     }
    dataframe_soils
    #creating a dataframe of just the soil characteristics
    fixed_field_data_processed_soils.condensed <-  st_drop_geometry(dataframe_soils)
    fixed_field_data_processed_soils.condensed <- fixed_field_data_processed_soils.condensed[, c(33:40, 49:62)]
    #creating a dataframe of just the size variables
    fixed_field_data_processed_size_variables <- st_drop_geometry(dataframe_soils)
    fixed_field_data_processed_size_variables <- fixed_field_data_processed_size_variables[, c(22:23, 29, 24, 20)]
    
    for (j in 1:(ncol(fixed_field_data_processed_soils.condensed))){ #iterating over the number of soil variables
      
      # soil_characteristics <- as.factor(soil_characteristics_list[i])
      
      for (k in 1:(ncol(fixed_field_data_processed_size_variables))) {#iterating over the number of soil variables
        #extracting slopes from comparing soil values with randomized shape/size values with linear regressions
        
        
        slopes <- c() #creating empty list to collect slope values
        
        set.seed(21)
        for (y in 1:1000){ #for 1000 permutations
          fixed_field_data_processed_soils_shuffled <- transform(dataframe_soils, variable.shuffled = sample(fixed_field_data_processed_size_variables[,k])) #create a data frame with a shuffled
          lm <- lm(fixed_field_data_processed_soils_shuffled$variable.shuffled ~ fixed_field_data_processed_soils.condensed[,j]) #LM_fixed_field_data_processed_soils.condensed[,i]
          lm_sum <- summary(lm) #extracting the linear regression information
          slopes <- c(slopes, lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
        }
        
        #extracting the slope of our points
        Size_Variable = colnames(fixed_field_data_processed_soils.condensed)[j] #as.factor(Size_Variable)
        lm_real <- lm(fixed_field_data_processed_size_variables[,k] ~ fixed_field_data_processed_soils.condensed[,j]) #creating the linear regression
        lm_real_sum <- summary(lm_real) #extract the summary 
        lm_real_slope <- lm_real_sum$coefficients[2] #storing the slope
        
        #plotting the histogram of the randomly distributed p-values and our real slope
        ggplot()+
          geom_histogram(aes(x=slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
          geom_vline(xintercept=lm_real_slope, col = "red") + #line of our real slope
          xlab(paste(Population, "Slopes of Shuffled", Size_Variable, "vs. our", Size_Variable)) +
          theme_classic()
        
        
        #calculating pseudo p-value for 
        total = 0  #set empty value
        for (p in 1:length(slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
          if (slopes[p] > lm_real_slope){
            total = total + 1
          }
        } #add number of values of in the random set of ANN values that are less than our mean ANN
        p_value <- (total / length(slopes)) #the proportion of random ANNs that are less than our ANN
        
        slopes_array[j, k, i] = lm_real_slope #assigning the each index in the array with the new slope value
        p_values_array[j, k, i] = p_value #assigning the each index in the array with the new p value
        
        print(paste("Updating:", j, k, i)) #soil metric, size metric, population
        print(paste("Real slope:", lm_real_slope))
        print(paste("P-value:", p_value))
        
      }
    }
    
  }
  
  return(list(slopes_array = slopes_array, p_values_array = p_values_array))
  
}

#running the simulation function
slope_simultations <- slopes_simulations()

#turning the simulation results into a data frame
slopes_df <- as.data.frame.table(slope_simultations$slopes_array, 
                                 responseName = "slope") #creating a dataframe from the slope array
pvals_df <- as.data.frame.table(slope_simultations$p_values_array, 
                                responseName = "p_value") #creating a dataframe from the p value array

size.pop.slopes.df <- merge(slopes_df, pvals_df, by = c("Var3", "Var2", "Var1")) #merging the two dataframes into one
names(size.pop.slopes.df) <- c("Population", "Size.Variable", "Soil.Metric", "Slope", "P.value") #re-naming the columns to be more appropriate

#creating a column for whether the p values are significant or not
size.pop.slopes.df <- size.pop.slopes.df %>%
  mutate(Significance = case_when(P.value < 0.05 ~ "Y",
                                  P.value >= 0.05 ~ "N"))

#ensuring the categorical variables have the desired order of levels
size.pop.slopes.df$Population <- factor(size.pop.slopes.df$Population, levels = c("LM", "LC", "SD"))
size.pop.slopes.df$Size.Variable <- factor(size.pop.slopes.df$Size.Variable, levels = c("SCA", "LCA", "CA", "CS", "DBH"))
size.pop.slopes.df$Soil.Metric <- factor(size.pop.slopes.df$Soil.Metric, levels = c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                                                                                    "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                                                                                    "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                                                                                    "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                                                                                    "Volume of water content -1500 kpa 100-200", 
                                                                                    "Nitrogen 0-5", "Nitrogen 100-200", 
                                                                                    "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                                                                                    "Sand Available Water 0-5", "Sand Available Water 100-200",
                                                                                    "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200"))

#rearranging the rows based on the specified levels for the categorical variables
size.pop.slopes.df <- size.pop.slopes.df %>%
  arrange(Population, Size.Variable, Soil.Metric)
  

#labeled p-values
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P Values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, P_Value, NA))) +
  coord_flip()


ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P Values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, P_Value, NA))) +
  coord_flip()

# Making a dataframe from all of these slopes and p-values

p_values <- c(LM_sca_clay_0.5_p_value, LM_sca_clay_100_200_p_value, LM_sca_silt_0.5_p_value, LM_sca_silt_100_200_p_value, LM_sca_sand_0.5_p_value, LM_sca_sand_100_200_p_value,
              LM_sca_ph_0.5_p_value, LM_sca_ph_100_200_p_value, LM_sca_soc_0.5_p_value, LM_sca_soc_100_200_p_value, LM_sca_vol_10_0.5_p_value, LM_sca_vol_10_100_200_p_value, 
              LM_sca_vol_1500_0.5_p_value, LM_sca_vol_1500_100_200_p_value, LM_sca_nitrogen_0.5_p_value, LM_sca_nitrogen_100_200_p_value, LM_lca_clay_0.5_p_value, LM_lca_clay_100_200_p_value,
              LM_sca_clay_loam_avail_water_0.5_p_value, LM_sca_clay_loam_avail_water_100.200_p_value,
              LM_sca_sandy_avail_water_0.5_p_value, LM_sca_sandy_avail_water_100.200_p_value,
              LM_lca_silt_0.5_p_value, LM_lca_silt_100_200_p_value, LM_lca_sand_0.5_p_value, LM_lca_sand_100_200_p_value, LM_lca_ph_0.5_p_value, LM_lca_ph_100_200_p_value, 
              LM_lca_soc_0.5_p_value, LM_lca_soc_100_200_p_value, LM_lca_vol_10_0.5_p_value, LM_lca_vol_10_100_200_p_value, LM_lca_vol_1500_0.5_p_value, LM_lca_vol_1500_100_200_p_value,
              LM_lca_nitrogen_0.5_p_value, LM_lca_nitrogen_100_200_p_value, 
              LM_lca_clay_loam_avail_water_0.5_p_value, LM_lca_clay_loam_avail_water_100.200_p_value,
              LM_lca_sandy_avail_water_0.5_p_value, LM_lca_sandy_avail_water_100.200_p_value,
              LM_ca_clay_0.5_p_value, LM_ca_clay_100_200_p_value, LM_ca_silt_0.5_p_value, LM_ca_silt_100_200_p_value,
              LM_ca_sand_0.5_p_value, LM_ca_sand_100_200_p_value, LM_ca_ph_0.5_p_value, LM_ca_ph_100_200_p_value, LM_ca_soc_0.5_p_value, LM_ca_soc_100_200_p_value, LM_ca_vol_10_0.5_p_value,
              LM_ca_vol_10_100_200_p_value, LM_ca_vol_1500_0.5_p_value, LM_ca_vol_1500_100_200_p_value, LM_ca_nitrogen_0.5_p_value, LM_ca_nitrogen_100_200_p_value,
              LM_cs_clay_0.5_p_value, LM_cs_clay_100_200_p_value, LM_cs_silt_0.5_p_value, LM_cs_silt_100_200_p_value, LM_cs_sand_0.5_p_value, LM_cs_sand_100_200_p_value, 
              LM_cs_ph_0.5_p_value, LM_cs_ph_100_200_p_value, LM_cs_soc_0.5_p_value, LM_cs_soc_100_200_p_value, LM_cs_vol_10_0.5_p_value, LM_cs_vol_10_100_200_p_value,
              LM_cs_vol_1500_0.5_p_value, LM_cs_vol_1500_100_200_p_value, LM_cs_nitrogen_0.5_p_value, LM_cs_nitrogen_100_200_p_value, LM_dbh_clay_0.5_p_value, 
              LM_cs_clay_loam_avail_water_0.5_p_value, LM_cs_clay_loam_avail_water_100.200_p_value,
              LM_cs_sandy_avail_water_0.5_p_value, LM_cs_sandy_avail_water_100.200_p_value,
              LM_dbh_clay_100_200_p_value, LM_dbh_silt_0.5_p_value, LM_dbh_silt_100_200_p_value, LM_dbh_sand_0.5_p_value, LM_dbh_sand_100_200_p_value, LM_dbh_ph_0.5_p_value,
              LM_dbh_ph_100_200_p_value, LM_dbh_soc_0.5_p_value, LM_dbh_soc_100_200_p_value, LM_dbh_vol_10_0.5_p_value, LM_dbh_vol_10_100_200_p_value, LM_dbh_vol_1500_0.5_p_value, 
              LM_dbh_vol_1500_100_200_p_value, LM_dbh_nitrogen_0.5_p_value, LM_dbh_nitrogen_100_200_p_value,
              LM_dbh_clay_loam_avail_water_0.5_p_value, LM_dbh_clay_loam_avail_water_100.200_p_value,
              LM_dbh_sandy_avail_water_0.5_p_value, LM_dbh_sandy_avail_water_100.200_p_value,
              LC_sca_clay_0.5_p_value,LC_sca_clay_100_200_p_value, LC_sca_silt_0.5_p_value, LC_sca_silt_100_200_p_value, LC_sca_sand_0.5_p_value, LC_sca_sand_100_200_p_value,
              LC_sca_ph_0.5_p_value, LC_sca_ph_100_200_p_value, LC_sca_soc_0.5_p_value, LC_sca_soc_100_200_p_value, LC_sca_vol_10_0.5_p_value, LC_sca_vol_10_100_200_p_value, 
              LC_sca_vol_1500_0.5_p_value, LC_sca_vol_1500_100_200_p_value, LC_sca_nitrogen_0.5_p_value, LC_sca_nitrogen_100_200_p_value, LC_lca_clay_0.5_p_value, LC_lca_clay_100_200_p_value,
              LC_sca_clay_loam_avail_water_0.5_p_value, LC_sca_clay_loam_avail_water_100.200_p_value,
              LC_sca_sandy_avail_water_0.5_p_value, LC_sca_sandy_avail_water_100.200_p_value,
              LC_lca_silt_0.5_p_value, LC_lca_silt_100_200_p_value, LC_lca_sand_0.5_p_value, LC_lca_sand_100_200_p_value, LC_lca_ph_0.5_p_value, LC_lca_ph_100_200_p_value, 
              LC_lca_soc_0.5_p_value, LC_lca_soc_100_200_p_value, LC_lca_vol_10_0.5_p_value, LC_lca_vol_10_100_200_p_value, LC_lca_vol_1500_0.5_p_value, LC_lca_vol_1500_100_200_p_value,
              LC_lca_nitrogen_0.5_p_value, LC_lca_nitrogen_100_200_p_value, LC_ca_clay_0.5_p_value, LC_ca_clay_100_200_p_value, LC_ca_silt_0.5_p_value, LC_ca_silt_100_200_p_value,
              LC_lca_clay_loam_avail_water_0.5_p_value, LC_lca_clay_loam_avail_water_100.200_p_value,
              LC_lca_sandy_avail_water_0.5_p_value, LC_lca_sandy_avail_water_100.200_p_value,
              LC_ca_sand_0.5_p_value, LC_ca_sand_100_200_p_value, LC_ca_ph_0.5_p_value, LC_ca_ph_100_200_p_value, LC_ca_soc_0.5_p_value, LC_ca_soc_100_200_p_value, LC_ca_vol_10_0.5_p_value,
              LC_ca_vol_10_100_200_p_value, LC_ca_vol_1500_0.5_p_value, LC_ca_vol_1500_100_200_p_value, LC_ca_nitrogen_0.5_p_value, LC_ca_nitrogen_100_200_p_value,
              LC_ca_clay_loam_avail_water_0.5_p_value, LC_ca_clay_loam_avail_water_100.200_p_value,
              LC_ca_sandy_avail_water_0.5_p_value, LC_ca_sandy_avail_water_100.200_p_value,
               LC_cs_clay_0.5_p_value, LC_cs_clay_100_200_p_value, LC_cs_silt_0.5_p_value, LC_cs_silt_100_200_p_value, LC_cs_sand_0.5_p_value, LC_cs_sand_100_200_p_value, 
              LC_cs_ph_0.5_p_value, LC_cs_ph_100_200_p_value, LC_cs_soc_0.5_p_value, LC_cs_soc_100_200_p_value, LC_cs_vol_10_0.5_p_value, LC_cs_vol_10_100_200_p_value,
              LC_cs_vol_1500_0.5_p_value, LC_cs_vol_1500_100_200_p_value, LC_cs_nitrogen_0.5_p_value, LC_cs_nitrogen_100_200_p_value, LC_dbh_clay_0.5_p_value, 
              LC_cs_clay_loam_avail_water_0.5_p_value, LC_cs_clay_loam_avail_water_100.200_p_value,
              LC_cs_sandy_avail_water_0.5_p_value, LC_cs_sandy_avail_water_100.200_p_value,
               LC_dbh_clay_100_200_p_value, LC_dbh_silt_0.5_p_value, LC_dbh_silt_100_200_p_value, LC_dbh_sand_0.5_p_value, LC_dbh_sand_100_200_p_value, LC_dbh_ph_0.5_p_value,
              LC_dbh_ph_100_200_p_value, LC_dbh_soc_0.5_p_value, LC_dbh_soc_100_200_p_value, LC_dbh_vol_10_0.5_p_value, LC_dbh_vol_10_100_200_p_value, LC_dbh_vol_1500_0.5_p_value, 
              LC_dbh_vol_1500_100_200_p_value, LC_dbh_nitrogen_0.5_p_value, LC_dbh_nitrogen_100_200_p_value,
              LC_dbh_clay_loam_avail_water_0.5_p_value, LC_dbh_clay_loam_avail_water_100.200_p_value,
              LC_dbh_sandy_avail_water_0.5_p_value, LC_dbh_sandy_avail_water_100.200_p_value,
               SD_sca_clay_0.5_p_value,SD_sca_clay_100_200_p_value, SD_sca_silt_0.5_p_value, SD_sca_silt_100_200_p_value, SD_sca_sand_0.5_p_value, SD_sca_sand_100_200_p_value,
              SD_sca_ph_0.5_p_value, SD_sca_ph_100_200_p_value, SD_sca_soc_0.5_p_value, SD_sca_soc_100_200_p_value, SD_sca_vol_10_0.5_p_value, SD_sca_vol_10_100_200_p_value, 
              SD_sca_vol_1500_0.5_p_value, SD_sca_vol_1500_100_200_p_value, SD_sca_nitrogen_0.5_p_value, SD_sca_nitrogen_100_200_p_value, SD_lca_clay_0.5_p_value, SD_lca_clay_100_200_p_value,
              SD_sca_clay_loam_avail_water_0.5_p_value, SD_sca_clay_loam_avail_water_100.200_p_value,
              SD_sca_sandy_avail_water_0.5_p_value, SD_sca_sandy_avail_water_100.200_p_value,
              SD_lca_silt_0.5_p_value, SD_lca_silt_100_200_p_value, SD_lca_sand_0.5_p_value, SD_lca_sand_100_200_p_value, SD_lca_ph_0.5_p_value, SD_lca_ph_100_200_p_value, 
              SD_lca_soc_0.5_p_value, SD_lca_soc_100_200_p_value, SD_lca_vol_10_0.5_p_value, SD_lca_vol_10_100_200_p_value, SD_lca_vol_1500_0.5_p_value, SD_lca_vol_1500_100_200_p_value,
              SD_lca_nitrogen_0.5_p_value, SD_lca_nitrogen_100_200_p_value, SD_ca_clay_0.5_p_value, SD_ca_clay_100_200_p_value, SD_ca_silt_0.5_p_value, SD_ca_silt_100_200_p_value,
              SD_lca_clay_loam_avail_water_0.5_p_value, SD_lca_clay_loam_avail_water_100.200_p_value,
              SD_lca_sandy_avail_water_0.5_p_value, SD_lca_sandy_avail_water_100.200_p_value,
              SD_ca_sand_0.5_p_value, SD_ca_sand_100_200_p_value, SD_ca_ph_0.5_p_value, SD_ca_ph_100_200_p_value, SD_ca_soc_0.5_p_value, SD_ca_soc_100_200_p_value, SD_ca_vol_10_0.5_p_value,
              SD_ca_vol_10_100_200_p_value, SD_ca_vol_1500_0.5_p_value, SD_ca_vol_1500_100_200_p_value, SD_ca_nitrogen_0.5_p_value, SD_ca_nitrogen_100_200_p_value,
              SD_ca_clay_loam_avail_water_0.5_p_value, SD_ca_clay_loam_avail_water_100.200_p_value,
              SD_ca_sandy_avail_water_0.5_p_value, SD_ca_sandy_avail_water_100.200_p_value,
              SD_cs_clay_0.5_p_value, SD_cs_clay_100_200_p_value, SD_cs_silt_0.5_p_value, SD_cs_silt_100_200_p_value, SD_cs_sand_0.5_p_value, SD_cs_sand_100_200_p_value, 
              SD_cs_ph_0.5_p_value, SD_cs_ph_100_200_p_value, SD_cs_soc_0.5_p_value, SD_cs_soc_100_200_p_value, SD_cs_vol_10_0.5_p_value, SD_cs_vol_10_100_200_p_value,
              SD_cs_vol_1500_0.5_p_value, SD_cs_vol_1500_100_200_p_value, SD_cs_nitrogen_0.5_p_value, SD_cs_nitrogen_100_200_p_value, SD_dbh_clay_0.5_p_value, 
              SD_cs_clay_loam_avail_water_0.5_p_value, SD_cs_clay_loam_avail_water_100.200_p_value,
              SD_cs_sandy_avail_water_0.5_p_value, SD_cs_sandy_avail_water_100.200_p_value,
              SD_dbh_clay_100_200_p_value, SD_dbh_silt_0.5_p_value, SD_dbh_silt_100_200_p_value, SD_dbh_sand_0.5_p_value, SD_dbh_sand_100_200_p_value, SD_dbh_ph_0.5_p_value,
              SD_dbh_ph_100_200_p_value, SD_dbh_soc_0.5_p_value, SD_dbh_soc_100_200_p_value, SD_dbh_vol_10_0.5_p_value, SD_dbh_vol_10_100_200_p_value, SD_dbh_vol_1500_0.5_p_value, 
              SD_dbh_vol_1500_100_200_p_value, SD_dbh_nitrogen_0.5_p_value, SD_dbh_nitrogen_100_200_p_value,
              SD_dbh_clay_loam_avail_water_0.5_p_value, SD_dbh_clay_loam_avail_water_100.200_p_value,
              SD_dbh_sandy_avail_water_0.5_p_value, SD_dbh_sandy_avail_water_100.200_p_value)
              

slopes <- c(LM_sca_clay_0_5_lm_real_slope, LM_sca_clay_100_200_lm_real_slope, LM_sca_silt_0_5_lm_real_slope, LM_sca_silt_100_200_lm_real_slope, LM_sca_sand_0_5_lm_real_slope, LM_sca_sand_100_200_lm_real_slope,
              LM_sca_ph_0_5_lm_real_slope, LM_sca_ph_100_200_lm_real_slope, LM_sca_soc_0_5_lm_real_slope, LM_sca_soc_100_200_lm_real_slope, LM_sca_vol_10_0_5_lm_real_slope, LM_sca_vol_10_100_200_lm_real_slope, 
              LM_sca_vol_1500_0_5_lm_real_slope, LM_sca_vol_1500_100_200_lm_real_slope, LM_sca_nitrogen_0_5_lm_real_slope, LM_sca_nitrogen_100_200_lm_real_slope, LM_lca_clay_0_5_lm_real_slope, LM_lca_clay_100_200_lm_real_slope,
             LM_sca_clay_loam_avail_water_0.5_lm_real_slope, LM_sca_clay_loam_avail_water_100.200_lm_real_slope,
             LM_sca_sandy_avail_water_0.5_lm_real_slope, LM_sca_sandy_avail_water_100.200_lm_real_slope,
             LM_lca_silt_0_5_lm_real_slope, LM_lca_silt_100_200_lm_real_slope, LM_lca_sand_0_5_lm_real_slope, LM_lca_sand_100_200_lm_real_slope, LM_lca_ph_0_5_lm_real_slope, LM_lca_ph_100_200_lm_real_slope, 
              LM_lca_soc_0_5_lm_real_slope, LM_lca_soc_100_200_lm_real_slope, LM_lca_vol_10_0_5_lm_real_slope, LM_lca_vol_10_100_200_lm_real_slope, LM_lca_vol_1500_0_5_lm_real_slope, LM_lca_vol_1500_100_200_lm_real_slope,
              LM_lca_nitrogen_0_5_lm_real_slope, LM_lca_nitrogen_100_200_lm_real_slope, LM_ca_clay_0_5_lm_real_slope, LM_ca_clay_100_200_lm_real_slope, LM_ca_silt_0_5_lm_real_slope, LM_ca_silt_100_200_lm_real_slope,
            LM_lca_clay_loam_avail_water_0.5_lm_real_slope, LM_lca_clay_loam_avail_water_100.200_lm_real_slope,
            LM_lca_sandy_avail_water_0.5_lm_real_slope, LM_lca_sandy_avail_water_100.200_lm_real_slope, 
             LM_ca_sand_0_5_lm_real_slope, LM_ca_sand_100_200_lm_real_slope, LM_ca_ph_0_5_lm_real_slope, LM_ca_ph_100_200_lm_real_slope, LM_ca_soc_0_5_lm_real_slope, LM_ca_soc_100_200_lm_real_slope, LM_ca_vol_10_0_5_lm_real_slope,
              LM_ca_vol_10_100_200_lm_real_slope, LM_ca_vol_1500_0_5_lm_real_slope, LM_ca_vol_1500_100_200_lm_real_slope, LM_ca_nitrogen_0_5_lm_real_slope, LM_ca_nitrogen_100_200_lm_real_slope,
            LM_ca_clay_loam_avail_water_0.5_lm_real_slope, LM_ca_clay_loam_avail_water_100.200_lm_real_slope,
            LM_ca_sandy_avail_water_0.5_lm_real_slope, LM_ca_sandy_avail_water_100.200_lm_real_slope,  
             LM_cs_clay_0_5_lm_real_slope, LM_cs_clay_100_200_lm_real_slope, LM_cs_silt_0_5_lm_real_slope, LM_cs_silt_100_200_lm_real_slope, LM_cs_sand_0_5_lm_real_slope, LM_cs_sand_100_200_lm_real_slope, 
              LM_cs_ph_0_5_lm_real_slope, LM_cs_ph_100_200_lm_real_slope, LM_cs_soc_0_5_lm_real_slope, LM_cs_soc_100_200_lm_real_slope, LM_cs_vol_10_0_5_lm_real_slope, LM_cs_vol_10_100_200_lm_real_slope,
              LM_cs_vol_1500_0_5_lm_real_slope, LM_cs_vol_1500_100_200_lm_real_slope, LM_cs_nitrogen_0_5_lm_real_slope, LM_cs_nitrogen_100_200_lm_real_slope, LM_dbh_clay_0_5_lm_real_slope, 
            LM_cs_clay_loam_avail_water_0.5_lm_real_slope, LM_cs_clay_loam_avail_water_100.200_lm_real_slope,
            LM_cs_sandy_avail_water_0.5_lm_real_slope, LM_cs_sandy_avail_water_100.200_lm_real_slope,    
            LM_dbh_clay_100_200_lm_real_slope, LM_dbh_silt_0_5_lm_real_slope, LM_dbh_silt_100_200_lm_real_slope, LM_dbh_sand_0_5_lm_real_slope, LM_dbh_sand_100_200_lm_real_slope, LM_dbh_ph_0_5_lm_real_slope,
              LM_dbh_ph_100_200_lm_real_slope, LM_dbh_soc_0_5_lm_real_slope, LM_dbh_soc_100_200_lm_real_slope, LM_dbh_vol_10_0_5_lm_real_slope, LM_dbh_vol_10_100_200_lm_real_slope, LM_dbh_vol_1500_0_5_lm_real_slope, 
              LM_dbh_vol_1500_100_200_lm_real_slope, LM_dbh_nitrogen_0_5_lm_real_slope, LM_dbh_nitrogen_100_200_lm_real_slope,
            LM_dbh_clay_loam_avail_water_0.5_lm_real_slope, LM_dbh_clay_loam_avail_water_100.200_lm_real_slope,
            LM_dbh_sandy_avail_water_0.5_lm_real_slope, LM_dbh_sandy_avail_water_100.200_lm_real_slope,   
            LC_sca_clay_0_5_LC_real_slope,LC_sca_clay_100_200_LC_real_slope, LC_sca_silt_0_5_LC_real_slope, LC_sca_silt_100_200_LC_real_slope, LC_sca_sand_0_5_LC_real_slope, LC_sca_sand_100_200_LC_real_slope,
              LC_sca_ph_0_5_LC_real_slope, LC_sca_ph_100_200_LC_real_slope, LC_sca_soc_0_5_lm_real_slope, LC_sca_soc_100_200_lm_real_slope, LC_sca_vol_10_0_5_lm_real_slope, LC_sca_vol_10_100_200_lm_real_slope, 
              LC_sca_vol_1500_0_5_lm_real_slope, LC_sca_vol_1500_100_200_lm_real_slope, LC_sca_nitrogen_0_5_lm_real_slope, LC_sca_nitrogen_100_200_lm_real_slope, LC_lca_clay_0_5_lm_real_slope, LC_lca_clay_100_200_lm_real_slope,
            LC_sca_clay_loam_avail_water_0.5_lm_real_slope, LC_sca_clay_loam_avail_water_100.200_lm_real_slope,
            LC_sca_sandy_avail_water_0.5_lm_real_slope, LC_sca_sandy_avail_water_100.200_lm_real_slope,  
            LC_lca_silt_0_5_lm_real_slope, LC_lca_silt_100_200_lm_real_slope, LC_lca_sand_0_5_lm_real_slope, LC_lca_sand_100_200_lm_real_slope, LC_lca_ph_0_5_lm_real_slope, LC_lca_ph_100_200_lm_real_slope, 
              LC_lca_soc_0_5_lm_real_slope, LC_lca_soc_100_200_lm_real_slope, LC_lca_vol_10_0_5_lm_real_slope, LC_lca_vol_10_100_200_lm_real_slope, LC_lca_vol_1500_0_5_lm_real_slope, LC_lca_vol_1500_100_200_lm_real_slope,
              LC_lca_nitrogen_0_5_lm_real_slope, LC_lca_nitrogen_100_200_lm_real_slope, LC_ca_clay_0_5_lm_real_slope, LC_ca_clay_100_200_lm_real_slope, LC_ca_silt_0_5_lm_real_slope, LC_ca_silt_100_200_lm_real_slope,
            LC_lca_clay_loam_avail_water_0.5_lm_real_slope, LC_lca_clay_loam_avail_water_100.200_lm_real_slope,
            LC_lca_sandy_avail_water_0.5_lm_real_slope, LC_lca_sandy_avail_water_100.200_lm_real_slope,    
            LC_ca_sand_0_5_lm_real_slope, LC_ca_sand_100_200_lm_real_slope, LC_ca_ph_0_5_lm_real_slope, LC_ca_ph_100_200_lm_real_slope, LC_ca_soc_0_5_lm_real_slope, LC_ca_soc_100_200_lm_real_slope, LC_ca_vol_10_0_5_lm_real_slope,
              LC_ca_vol_10_100_200_lm_real_slope, LC_ca_vol_1500_0_5_lm_real_slope, LC_ca_vol_1500_100_200_lm_real_slope, LC_ca_nitrogen_0_5_lm_real_slope, LC_ca_nitrogen_100_200_lm_real_slope,
            LC_ca_clay_loam_avail_water_0.5_lm_real_slope, LC_ca_clay_loam_avail_water_100.200_lm_real_slope,
            LC_ca_sandy_avail_water_0.5_lm_real_slope, LC_ca_sandy_avail_water_100.200_lm_real_slope,      
            LC_cs_clay_0_5_lm_real_slope, LC_cs_clay_100_200_lm_real_slope, LC_cs_silt_0_5_lm_real_slope, LC_cs_silt_100_200_lm_real_slope, LC_cs_sand_0_5_lm_real_slope, LC_cs_sand_100_200_lm_real_slope, 
              LC_cs_ph_0_5_lm_real_slope, LC_cs_ph_100_200_lm_real_slope, LC_cs_soc_0_5_lm_real_slope, LC_cs_soc_100_200_lm_real_slope, LC_cs_vol_10_0_5_lm_real_slope, LC_cs_vol_10_100_200_lm_real_slope,
              LC_cs_vol_1500_0_5_lm_real_slope, LC_cs_vol_1500_100_200_lm_real_slope, LC_cs_nitrogen_0_5_lm_real_slope, LC_cs_nitrogen_100_200_lm_real_slope, LC_dbh_clay_0_5_lm_real_slope, 
            LC_cs_clay_loam_avail_water_0.5_lm_real_slope, LC_cs_clay_loam_avail_water_100.200_lm_real_slope,
            LC_cs_sandy_avail_water_0.5_lm_real_slope, LC_cs_sandy_avail_water_100.200_lm_real_slope,  
             LC_dbh_clay_100_200_lm_real_slope, LC_dbh_silt_0_5_lm_real_slope, LC_dbh_silt_100_200_lm_real_slope, LC_dbh_sand_0_5_lm_real_slope, LC_dbh_sand_100_200_lm_real_slope, LC_dbh_ph_0_5_lm_real_slope,
              LC_dbh_ph_100_200_lm_real_slope, LC_dbh_soc_0_5_lm_real_slope, LC_dbh_soc_100_200_lm_real_slope, LC_dbh_vol_10_0_5_lm_real_slope, LC_dbh_vol_10_100_200_lm_real_slope, LC_dbh_vol_1500_0_5_lm_real_slope, 
              LC_dbh_vol_1500_100_200_lm_real_slope, LC_dbh_nitrogen_0_5_lm_real_slope, LC_dbh_nitrogen_100_200_lm_real_slope,
            LC_dbh_clay_loam_avail_water_0.5_lm_real_slope, LC_dbh_clay_loam_avail_water_100.200_lm_real_slope,
            LC_dbh_sandy_avail_water_0.5_lm_real_slope, LC_dbh_sandy_avail_water_100.200_lm_real_slope,   
            SD_sca_clay_0_5_SD_real_slope,SD_sca_clay_100_200_SD_real_slope, SD_sca_silt_0_5_SD_real_slope, SD_sca_silt_100_200_SD_real_slope, SD_sca_sand_0_5_SD_real_slope, SD_sca_sand_100_200_SD_real_slope,
              SD_sca_ph_0_5_SD_real_slope, SD_sca_ph_100_200_SD_real_slope, SD_sca_soc_0_5_lm_real_slope, SD_sca_soc_100_200_lm_real_slope, SD_sca_vol_10_0_5_lm_real_slope, SD_sca_vol_10_100_200_lm_real_slope, 
              SD_sca_vol_1500_0_5_lm_real_slope, SD_sca_vol_1500_100_200_lm_real_slope, SD_sca_nitrogen_0_5_lm_real_slope, SD_sca_nitrogen_100_200_lm_real_slope, SD_lca_clay_0_5_lm_real_slope, SD_lca_clay_100_200_lm_real_slope,
            SD_sca_clay_loam_avail_water_0.5_lm_real_slope, SD_sca_clay_loam_avail_water_100.200_lm_real_slope,
            SD_sca_sandy_avail_water_0.5_lm_real_slope, SD_sca_sandy_avail_water_100.200_lm_real_slope,    
            SD_lca_silt_0_5_lm_real_slope, SD_lca_silt_100_200_lm_real_slope, SD_lca_sand_0_5_lm_real_slope, SD_lca_sand_100_200_lm_real_slope, SD_lca_ph_0_5_lm_real_slope, SD_lca_ph_100_200_lm_real_slope, 
              SD_lca_soc_0_5_lm_real_slope, SD_lca_soc_100_200_lm_real_slope, SD_lca_vol_10_0_5_lm_real_slope, SD_lca_vol_10_100_200_lm_real_slope, SD_lca_vol_1500_0_5_lm_real_slope, SD_lca_vol_1500_100_200_lm_real_slope,
              SD_lca_nitrogen_0_5_lm_real_slope, SD_lca_nitrogen_100_200_lm_real_slope, SD_ca_clay_0_5_lm_real_slope, SD_ca_clay_100_200_lm_real_slope, SD_ca_silt_0_5_lm_real_slope, SD_ca_silt_100_200_lm_real_slope,
            SD_lca_clay_loam_avail_water_0.5_lm_real_slope, SD_lca_clay_loam_avail_water_100.200_lm_real_slope,
            SD_lca_sandy_avail_water_0.5_lm_real_slope, SD_lca_sandy_avail_water_100.200_lm_real_slope,   
            SD_ca_sand_0_5_lm_real_slope, SD_ca_sand_100_200_lm_real_slope, SD_ca_ph_0_5_lm_real_slope, SD_ca_ph_100_200_lm_real_slope, SD_ca_soc_0_5_lm_real_slope, SD_ca_soc_100_200_lm_real_slope, SD_ca_vol_10_0_5_lm_real_slope,
              SD_ca_vol_10_100_200_lm_real_slope, SD_ca_vol_1500_0_5_lm_real_slope, SD_ca_vol_1500_100_200_lm_real_slope, SD_ca_nitrogen_0_5_lm_real_slope, SD_ca_nitrogen_100_200_lm_real_slope,
            SD_ca_clay_loam_avail_water_0.5_lm_real_slope, SD_ca_clay_loam_avail_water_100.200_lm_real_slope,
            SD_ca_sandy_avail_water_0.5_lm_real_slope, SD_ca_sandy_avail_water_100.200_lm_real_slope,   
            SD_cs_clay_0_5_lm_real_slope, SD_cs_clay_100_200_lm_real_slope, SD_cs_silt_0_5_lm_real_slope, SD_cs_silt_100_200_lm_real_slope, SD_cs_sand_0_5_lm_real_slope, SD_cs_sand_100_200_lm_real_slope, 
              SD_cs_ph_0_5_lm_real_slope, SD_cs_ph_100_200_lm_real_slope, SD_cs_soc_0_5_lm_real_slope, SD_cs_soc_100_200_lm_real_slope, SD_cs_vol_10_0_5_lm_real_slope, SD_cs_vol_10_100_200_lm_real_slope,
              SD_cs_vol_1500_0_5_lm_real_slope, SD_cs_vol_1500_100_200_lm_real_slope, SD_cs_nitrogen_0_5_lm_real_slope, SD_cs_nitrogen_100_200_lm_real_slope, SD_dbh_clay_0_5_lm_real_slope, 
            SD_cs_clay_loam_avail_water_0.5_lm_real_slope, SD_cs_clay_loam_avail_water_100.200_lm_real_slope,
            SD_cs_sandy_avail_water_0.5_lm_real_slope, SD_cs_sandy_avail_water_100.200_lm_real_slope,     
            SD_dbh_clay_100_200_lm_real_slope, SD_dbh_silt_0_5_lm_real_slope, SD_dbh_silt_100_200_lm_real_slope, SD_dbh_sand_0_5_lm_real_slope, SD_dbh_sand_100_200_lm_real_slope, SD_dbh_ph_0_5_lm_real_slope,
              SD_dbh_ph_100_200_lm_real_slope, SD_dbh_soc_0_5_lm_real_slope, SD_dbh_soc_100_200_lm_real_slope, SD_dbh_vol_10_0_5_lm_real_slope, SD_dbh_vol_10_100_200_lm_real_slope, SD_dbh_vol_1500_0_5_lm_real_slope, 
              SD_dbh_vol_1500_100_200_lm_real_slope, SD_dbh_nitrogen_0_5_lm_real_slope, SD_dbh_nitrogen_100_200_lm_real_slope,
            SD_dbh_clay_loam_avail_water_0.5_lm_real_slope, SD_dbh_clay_loam_avail_water_100.200_lm_real_slope,
            SD_dbh_sandy_avail_water_0.5_lm_real_slope, SD_dbh_sandy_avail_water_100.200_lm_real_slope)


size.pop.slopes.df <- data.frame("Population" = c(rep('LM', 80), rep('LC', 80), rep('SD', 80)),
                                 "Variable" = rep(c(rep("SCA", 16), rep("LCA", 16), rep("CA", 16), rep("CS", 16),
                                                  rep("DBH", 16)), 3),
                                 "Shape.Size" = rep(c("Clay 0-5 cm", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                                                      "Ph 0-5", "Ph 100-200", "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200", "Volume of water content -10 kpa 0-5",
                                                      "Volume of water content -10 kpa 100-200", "Volume of water content -1500 kpa 0-5", "Volume of water content -1500 kpa 100-200", 
                                                      "Nitrogen 0-5", "Nitrogen 100-200", "sandy_avail_water_0.5", "sandy_avail_water_100.200",
                                                      "clay_loam_avail_water_0.5", "clay_loam_avail_water_100.200")),
                                 "Slope" = slopes, 
                                 "P_Value" = p_values,
                                 "Significance" = ifelse(p_values < 0.05, "Y", "N"))
View(size.pop.slopes.df)


#Summarizing results across population, size/shape, and soil variables


#summarizing the data frame 
summary(size.pop.slopes.df)

#number of yes or no by population 
table(size.pop.slopes.df$Population, size.pop.slopes.df$Significance)

# number of yes or nos for each size variables 

#LM
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="SD"])

# number of yes or nos for each soil variables 

#LM
table(size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="SD"])


## number of yes or nos for each size variables and each soil variable 

#LM 
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="LM"])

#LC 
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="LC"])

#SD 
table(size.pop.slopes.df$Shape.Size[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Variable[size.pop.slopes.df$Population=="SD"])


#plots to visualize the results

#plot of population, yes/no, and the interaciton between size/shape and soil values
ggplot(size.pop.slopes.df, aes(x = Significance, fill = interaction(Shape.Size,Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and  size/shape 
ggplot(size.pop.slopes.df, aes(x = Significance, fill = interaction(Shape.Size)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and soil values
ggplot(size.pop.slopes.df, aes(x = Significance, fill = interaction(Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()


### Heat Map
size.pop.slopes.df$Population <- as.factor(size.pop.slopes.df$Population)
levels(size.pop.slopes.df$Population)

options(digits=3)

#across all populations
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Values",
       title = "All Sites Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) +
  coord_flip()

#across LM
size.pop.slopes.df.LM <- size.pop.slopes.df[size.pop.slopes.df$Population == "LM", ] #isolating the LM p-values/slopes

#labeled p-values
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, P_Value, NA))) +
  coord_flip() 

#labeled slopes
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, round(Slope, 3), NA)), color = "black") + 
  coord_flip()  



#across LC
size.pop.slopes.df.LC <- size.pop.slopes.df[size.pop.slopes.df$Population == "LC", ] #isolating the LM p-values/slopes

#labeled p-values
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(P_Value < 0.05, P_Value, NA))) +
  coord_flip()

#labeled slopes
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(P_Value < 0.05, round(Slope, 3), NA))) +
  coord_flip()



#across SD
size.pop.slopes.df.SD <- size.pop.slopes.df[size.pop.slopes.df$Population == "SD", ] #isolating the LM p-values/slopes

#labeled p-values
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P Values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, P_Value, NA))) +
  coord_flip()

#labeled slopes
ggplot(aes(x = Shape.Size, y = Variable, fill = ifelse(P_Value < 0.05, P_Value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.05, round(Slope, 3), NA))) +
  coord_flip()


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
 # geom_sf(data=bbox_poly_buffered)+
  geom_raster(data= as.data.frame(clay_05_all_pop, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = BCS_polygon_box_sf_cropped, fill = NA, color = "green")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates) +
  geom_sf(data=random_20, color ="red") + 
  labs(title = "Clay at 0-5 cm",
       fill = "Clay Content (g/kg)",
       x = "",
       y = "") +
  scale_fill_viridis_c(limits = c(50, 330))
vol_wat_33kpa_05_all_pop$vol_water_0.5
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

View(all_known_pop_soils)

#clay

#extracting means from randomly selected 16 points 

random_clay_0.5_means <- c() #creating empty list to collect means
random_clay_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations

  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
  random_20 <- random_20 %>%
    st_as_sf()
  random_20_pop_soil_clay <- raster::extract(soil_stack_clay, random_20) #extracting the soil metrics for the random points
  
  random_clay_0.5_mean <- mean(random_20_pop_soil_clay[,1]) #storing the mean of the 0-5 value
  random_clay_100.200_mean <- mean(random_20_pop_soil_clay[,2]) #storing the mean of the 100-200 value
  
  random_clay_0.5_means <- c(random_clay_0.5_means, random_clay_0.5_mean) #adding the 0-5 mean to the list of means
  random_clay_100.200_means <- c(random_clay_100.200_means, random_clay_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_clay_0.5_mean <- mean(all_known_pop_soils$clay.content.0.5)
all_known_clay_100.200_mean <- mean(all_known_pop_soils$clay.content.100.200)

#for Clay 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean Clay 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_clay_0.5_means <- na.omit(random_clay_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_clay_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_0.5_means[i] < all_known_clay_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
clay_0.5_random_p.value <- (total / length(random_clay_0.5_means)) #the proportion of random ANNs that are less than our ANN, our p-value

1- (total / length(random_clay_0.5_means)) #the proportion of random ANNs that are greater than our ANN


#for Clay 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean Clay 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_clay_100.200_means <- na.omit(random_clay_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_clay_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_100.200_means[i] > all_known_clay_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
clay_100.200_random_p.value <- (total / length(random_clay_100.200_means)) #the proportion of random ANNs that are less than our ANN
clay_100.200_random_p.value



#silt

#extracting means from randomly selected 16 points 

random_silt_0.5_means <- c() #creating empty list to collect means
random_silt_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_silt <- raster::extract(soil_stack_silt, random_20) #extracting the soil metrics for the random points
  
  random_silt_0.5_mean <- mean(random_20_pop_soil_silt[,1]) #storing the mean of the 0-5 value
  random_silt_100.200_mean <- mean(random_20_pop_soil_silt[,2]) #storing the mean of the 100-200 value
  
  random_silt_0.5_means <- c(random_silt_0.5_means, random_silt_0.5_mean) #adding the 0-5 mean to the list of means
  random_silt_100.200_means <- c(random_silt_100.200_means, random_silt_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_silt_0.5_mean <- mean(all_known_pop_soils$silt.0.5)
all_known_silt_100.200_mean <- mean(all_known_pop_soils$silt.100.200)

#for silt 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_silt_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50)+
  geom_vline(xintercept=all_known_silt_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean silt 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_silt_0.5_means <- na.omit(random_silt_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_silt_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_silt_0.5_means[i] < all_known_silt_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
silt_0.5_random_p.value <- (total / length(random_silt_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for silt 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_silt_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_silt_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean silt 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_silt_100.200_means <- na.omit(random_silt_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_silt_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_silt_100.200_means[i] < all_known_silt_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
silt_100.200_random_p.value <- (total / length(random_silt_100.200_means)) #the proportion of random ANNs that are less than our ANN
silt_100.200_random_p.value


#sand

#extracting means from randomly selected 20 points 

random_sand_0.5_means <- c() #creating empty list to collect means
random_sand_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_sand <- raster::extract(soil_stack_sand, random_20) #extracting the soil metrics for the random points
  
  random_sand_0.5_mean <- mean(random_20_pop_soil_sand[,1]) #storing the mean of the 0-5 value
  random_sand_100.200_mean <- mean(random_20_pop_soil_sand[,2]) #storing the mean of the 100-200 value
  
  random_sand_0.5_means <- c(random_sand_0.5_means, random_sand_0.5_mean) #adding the 0-5 mean to the list of means
  random_sand_100.200_means <- c(random_sand_100.200_means, random_sand_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_sand_0.5_mean <- mean(all_known_pop_soils$sand.0.5)
all_known_sand_100.200_mean <- mean(all_known_pop_soils$sand.100.200)

#for sand 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_sand_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_sand_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean sand 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_sand_0.5_means <- na.omit(random_sand_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_sand_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_sand_0.5_means[i] > all_known_sand_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
sand_0.5_random_p.value <- (total / length(random_sand_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for sand 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_sand_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_sand_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean sand 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_sand_100.200_means <- na.omit(random_sand_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_sand_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_sand_100.200_means[i] > all_known_sand_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
sand_100.200_random_p.value <- (total / length(random_sand_100.200_means)) #the proportion of random ANNs that are less than our ANN


#ph

#extracting means from randomly selected 20 points 

random_ph_0.5_means <- c() #creating empty list to collect means
random_ph_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_ph <- raster::extract(soil_stack_ph, random_20) #extracting the soil metrics for the random points

  random_ph_0.5_mean <- mean(random_20_pop_soil_ph[,1]) #storing the mean of the 0-5 value
  random_ph_100.200_mean <- mean(random_20_pop_soil_ph[,2]) #storing the mean of the 100-200 value
  
  random_ph_0.5_means <- c(random_ph_0.5_means, random_ph_0.5_mean) #adding the 0-5 mean to the list of means
  random_ph_100.200_means <- c(random_ph_100.200_means, random_ph_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_ph_0.5_mean <- mean(all_known_pop_soils$ph_0.5)
all_known_ph_100.200_mean <- mean(all_known_pop_soils$ph_100.200)

#for ph 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_ph_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_ph_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean ph 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_ph_0.5_means <- na.omit(random_ph_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_ph_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_ph_0.5_means[i] > all_known_ph_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
ph_0.5_random_p.value <- (total / length(random_ph_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for ph 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_ph_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_ph_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean ph 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_ph_100.200_means <- na.omit(random_ph_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_ph_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_ph_100.200_means[i] > all_known_ph_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
ph_100.200_random_p.value <- (total / length(random_ph_100.200_means)) #the proportion of random ANNs that are less than our ANN


#soc

#extracting means from randomly selected 20 points 

random_soc_0.5_means <- c() #creating empty list to collect means
random_soc_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>%
    st_as_sf()
  random_20_pop_soil_soc <- raster::extract(soil_stack_soc, random_20) #extracting the soil metrics for the random points
  
  random_soc_0.5_mean <- mean(random_20_pop_soil_soc[,1]) #storing the mean of the 0-5 value
  random_soc_100.200_mean <- mean(random_20_pop_soil_soc[,2]) #storing the mean of the 100-200 value
  
  random_soc_0.5_means <- c(random_soc_0.5_means, random_soc_0.5_mean) #adding the 0-5 mean to the list of means
  random_soc_100.200_means <- c(random_soc_100.200_means, random_soc_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_soc_0.5_mean <- mean(all_known_pop_soils$SOC.0.5)
all_known_soc_100.200_mean <- mean(all_known_pop_soils$SOC.100.200)

#for soc 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_soc_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_soc_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean soc 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_soc_0.5_means <- na.omit(random_soc_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_soc_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_soc_0.5_means[i] < all_known_soc_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
soc_0.5_random_p.value <- (total / length(random_soc_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for soc 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_soc_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_soc_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean soc 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_soc_100.200_means <- na.omit(random_soc_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_soc_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_soc_100.200_means[i] < all_known_soc_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
soc_100.200_random_p.value <- (total / length(random_soc_100.200_means)) #the proportion of random ANNs that are less than our ANN


#vol_wat_10kpa

#extracting means from randomly selected 20 points 

random_vol_wat_10kpa_0.5_means <- c() #creating empty list to collect means
random_vol_wat_10kpa_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_vol_wat_10kpa <- raster::extract(soil_stack_vol_wat_10kpa, random_20) #extracting the soil metrics for the random points

  random_vol_wat_10kpa_0.5_mean <- mean(random_20_pop_soil_vol_wat_10kpa[,1]) #storing the mean of the 0-5 value
  random_vol_wat_10kpa_100.200_mean <- mean(random_20_pop_soil_vol_wat_10kpa[,2]) #storing the mean of the 100-200 value
  
  random_vol_wat_10kpa_0.5_means <- c(random_vol_wat_10kpa_0.5_means, random_vol_wat_10kpa_0.5_mean) #adding the 0-5 mean to the list of means
  random_vol_wat_10kpa_100.200_means <- c(random_vol_wat_10kpa_100.200_means, random_vol_wat_10kpa_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_vol_wat_10kpa_0.5_mean <- mean(all_known_pop_soils$vol_water_.10_0.5)
all_known_vol_wat_10kpa_100.200_mean <- mean(all_known_pop_soils$vol_water_.10_100.200)

#for vol_wat_10kpa 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_10kpa_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_10kpa_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_10kpa 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_10kpa_0.5_means <- na.omit(random_vol_wat_10kpa_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_10kpa_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_10kpa_0.5_means[i] > all_known_vol_wat_10kpa_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_10kpa_0.5_random_p.value <- (total / length(random_vol_wat_10kpa_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for vol_wat_10kpa 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_10kpa_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_10kpa_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_10kpa 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_10kpa_100.200_means <- na.omit(random_vol_wat_10kpa_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_10kpa_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_10kpa_100.200_means[i] > all_known_vol_wat_10kpa_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_10kpa_100.200_random_p.value <- (total / length(random_vol_wat_10kpa_100.200_means)) #the proportion of random ANNs that are less than our ANN


#vol_wat_33kpa

#extracting means from randomly selected 20 points 

random_vol_wat_33kpa_0.5_means <- c() #creating empty list to collect means
random_vol_wat_33kpa_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_vol_wat_33kpa <- raster::extract(soil_stack_vol_wat_33kpa, random_20) #extracting the soil metrics for the random points
  
  random_vol_wat_33kpa_0.5_mean <- mean(random_20_pop_soil_vol_wat_33kpa[,1]) #storing the mean of the 0-5 value
  random_vol_wat_33kpa_100.200_mean <- mean(random_20_pop_soil_vol_wat_33kpa[,2]) #storing the mean of the 100-200 value
  
  random_vol_wat_33kpa_0.5_means <- c(random_vol_wat_33kpa_0.5_means, random_vol_wat_33kpa_0.5_mean) #adding the 0-5 mean to the list of means
  random_vol_wat_33kpa_100.200_means <- c(random_vol_wat_33kpa_100.200_means, random_vol_wat_33kpa_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")


#storing the real means
all_known_vol_wat_33kpa_0.5_mean <- mean(all_known_pop_soils$vol_water_0.5)
all_known_vol_wat_33kpa_100.200_mean <- mean(all_known_pop_soils$vol_water_100.200)

#for vol_wat_33kpa 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_33kpa_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_33kpa_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_33kpa 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_33kpa_0.5_means <- na.omit(random_vol_wat_33kpa_0.5_means) #remove NAs

#for presentation
as_tibble(random_vol_wat_33kpa_0.5_means) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(random_vol_wat_33kpa_0.5_means, all_known_vol_wat_33kpa_0.5_mean)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=all_known_vol_wat_33kpa_0.5_mean, col = "red", size = 1.2) + #plotting our tree's mean ANN
  xlab("Mean Clay/Loam Field Capacity 0-5 cm of Random Populations vs. Known Populations (n=20)") +
  ylab("Frequency") +
  theme_classic()+
  theme(axis.text=element_text(size=15),  axis.title.x =element_text(size= 13),
        axis.title.y =element_text(size= 13))


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_33kpa_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_33kpa_0.5_means[i] > all_known_vol_wat_33kpa_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_33kpa_0.5_random_p.value <- (total / length(random_vol_wat_33kpa_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for vol_wat_33kpa 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_33kpa_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_33kpa_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_33kpa 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_33kpa_100.200_means <- na.omit(random_vol_wat_33kpa_100.200_means) #remove NAs

#for presentation
as_tibble(random_vol_wat_33kpa_100.200_means) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(random_vol_wat_33kpa_100.200_means, all_known_vol_wat_33kpa_100.200_mean)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=all_known_vol_wat_33kpa_100.200_mean, col = "red", size = 1.2) + #plotting our tree's mean ANN
  xlab("Mean Clay/Loam Field Capacity 100-200 cm of Random Populations vs.Known Populations (n=20)") +
  ylab("Frequency") +
  theme_classic()+
  theme(axis.text=element_text(size=15),  axis.title.x =element_text(size= 13),
        axis.title.y =element_text(size= 13))

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_33kpa_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_33kpa_100.200_means[i] > all_known_vol_wat_33kpa_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_33kpa_100.200_random_p.value <- (total / length(random_vol_wat_33kpa_100.200_means)) #the proportion of random ANNs that are less than our ANN


#vol_wat_1500kpa

#extracting means from randomly selected 20 points 

random_vol_wat_1500kpa_0.5_means <- c() #creating empty list to collect means
random_vol_wat_1500kpa_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select rando 16 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_vol_wat_1500kpa <- raster::extract(soil_stack_vol_wat_1500kpa, random_20) #extracting the soil metrics for the random points

  
  random_vol_wat_1500kpa_0.5_mean <- mean(random_20_pop_soil_vol_wat_1500kpa[,1]) #storing the mean of the 0-5 value
  random_vol_wat_1500kpa_100.200_mean <- mean(random_20_pop_soil_vol_wat_1500kpa[,2]) #storing the mean of the 100-200 value
  
  random_vol_wat_1500kpa_0.5_means <- c(random_vol_wat_1500kpa_0.5_means, random_vol_wat_1500kpa_0.5_mean) #adding the 0-5 mean to the list of means
  random_vol_wat_1500kpa_100.200_means <- c(random_vol_wat_1500kpa_100.200_means, random_vol_wat_1500kpa_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_vol_wat_1500kpa_0.5_mean <- mean(all_known_pop_soils$vol_water_.1500kPa_0.5)
all_known_vol_wat_1500kpa_100.200_mean <- mean(all_known_pop_soils$vol_water_.1500_100.200)

#for vol_wat_1500kpa 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_1500kpa_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_1500kpa_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_1500kpa 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_1500kpa_0.5_means <- na.omit(random_vol_wat_1500kpa_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_1500kpa_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_1500kpa_0.5_means[i] > all_known_vol_wat_1500kpa_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_1500kpa_0.5_random_p.value <- (total / length(random_vol_wat_1500kpa_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for vol_wat_1500kpa 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_1500kpa_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_1500kpa_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_1500kpa 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_vol_wat_1500kpa_100.200_means <- na.omit(random_vol_wat_1500kpa_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_1500kpa_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_1500kpa_100.200_means[i] > all_known_vol_wat_1500kpa_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
vol_wat_1500kpa_100.200_random_p.value <- (total / length(random_vol_wat_1500kpa_100.200_means)) #the proportion of random ANNs that are less than our ANN

#nitrogen

#extracting means from randomly selected 20 points 

random_nitrogen_0.5_means <- c() #creating empty list to collect means
random_nitrogen_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_nitrogen <- raster::extract(soil_stack_nitrogen, random_20) #extracting the soil metrics for the random points

  
  random_nitrogen_0.5_mean <- mean(random_20_pop_soil_nitrogen[,1]) #storing the mean of the 0-5 value
  random_nitrogen_100.200_mean <- mean(random_20_pop_soil_nitrogen[,2]) #storing the mean of the 100-200 value
  
  random_nitrogen_0.5_means <- c(random_nitrogen_0.5_means, random_nitrogen_0.5_mean) #adding the 0-5 mean to the list of means
  random_nitrogen_100.200_means <- c(random_nitrogen_100.200_means, random_nitrogen_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_nitrogen_0.5_mean <- mean(all_known_pop_soils$nitrogen.0.5)
all_known_nitrogen_100.200_mean <- mean(all_known_pop_soils$nitrogen.100.200)

#for nitrogen 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_nitrogen_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_nitrogen_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean nitrogen 0-5 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_nitrogen_0.5_means <- na.omit(random_nitrogen_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_nitrogen_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_nitrogen_0.5_means[i] < all_known_nitrogen_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
nitrogen_0.5_random_p.value <- (total / length(random_nitrogen_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for nitrogen 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_nitrogen_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_nitrogen_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean nitrogen 100-200 of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_nitrogen_100.200_means <- na.omit(random_nitrogen_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_nitrogen_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_nitrogen_100.200_means[i] < all_known_nitrogen_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
nitrogen_100.200_random_p.value <- (total / length(random_nitrogen_100.200_means)) #the proportion of random ANNs that are less than our ANN

#for sandy available water 

#extracting means from randomly selected 20 points 

random_sandy_avail_water_0.5_means <- c() #creating empty list to collect means
random_sandy_avail_water_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_sandy_water <- raster::extract(soil_stack_sandy_water, random_20) #extracting the soil metrics for the random points
  
  
  random_sandy_avail_water_0.5_mean <- mean(random_20_pop_soil_sandy_water[,1]) #storing the mean of the 0-5 value
  random_sandy_avail_water_100.200_mean <- mean(random_20_pop_soil_sandy_water[,2]) #storing the mean of the 100-200 value
  
  random_sandy_avail_water_0.5_means <- c(random_sandy_avail_water_0.5_means, random_sandy_avail_water_0.5_mean) #adding the 0-5 mean to the list of means
  random_sandy_avail_water_100.200_means <- c(random_sandy_avail_water_100.200_means, random_sandy_avail_water_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_sandy_avail_water_0.5_mean <- mean(all_known_pop_soils$sandy_avail_water_0.5)
all_known_sandy_avail_water_100.200_mean <- mean(all_known_pop_soils$sandy_avail_water_100.200)

#for sandy avail water 0.5 cm

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_sandy_avail_water_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_sandy_avail_water_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean sandy available water 0-5 cm of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

#for presentation
as_tibble(random_sandy_avail_water_0.5_means) %>% #turning the ann.r vector as a tibble
  ggplot()+
  geom_histogram(aes(x = value), fill = "skyblue", color = "black", bins = 50) + 
  xlim(range(random_sandy_avail_water_0.5_means, all_known_sandy_avail_water_0.5_mean)) + #setting the range of the graph to include both the simulated ANN and our tree's mean ANN
  geom_vline(xintercept=all_known_sandy_avail_water_0.5_mean, col = "red", size = 1.2) + #plotting our tree's mean ANN
  xlab("Mean Sand Available Water 0-5 cm of Random Populations vs.Known Populations (n=20)")+
  ylab("Frequency") +
  theme_classic()+
  theme(axis.text=element_text(size=15),  axis.title.x =element_text(size= 13),
        axis.title.y =element_text(size= 13))

random_sandy_avail_water_0.5_means <- na.omit(random_sandy_avail_water_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_sandy_avail_water_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_sandy_avail_water_0.5_means[i] > all_known_sandy_avail_water_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
sandy_avail_water_0.5_random_p.value <- (total / length(random_sandy_avail_water_0.5_means)) #the proportion of random ANNs that are less than our ANN


#for sandy available water 100-200 cm

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_sandy_avail_water_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_sandy_avail_water_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean sandy available water 100-200 cm of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_sandy_avail_water_100.200_means <- na.omit(random_sandy_avail_water_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_sandy_avail_water_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_sandy_avail_water_100.200_means[i] > all_known_sandy_avail_water_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
sandy_avail_water_100.200_random_p.value <- (total / length(random_sandy_avail_water_100.200_means)) #the proportion of random ANNs that are less than our ANN

#clay loam available water 

#extracting means from randomly selected 20 points 

random_clay_loam_avail_water_0.5_means <- c() #creating empty list to collect means
random_clay_loam_avail_water_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
  random_20 <- random_20 %>% #turning the points into an sf object
    st_as_sf()
  random_20_pop_soil_clay_loam_water <- raster::extract(soil_stack_clay_loam_water, random_20) #extracting the soil metrics for the random points
  
  
  random_clay_loam_avail_water_0.5_mean <- mean(random_20_pop_soil_clay_loam_water[,1]) #storing the mean of the 0-5 value
  random_clay_loam_avail_water_100.200_mean <- mean(random_20_pop_soil_clay_loam_water[,2]) #storing the mean of the 100-200 value
  
  random_clay_loam_avail_water_0.5_means <- c(random_clay_loam_avail_water_0.5_means, random_clay_loam_avail_water_0.5_mean) #adding the 0-5 mean to the list of means
  random_clay_loam_avail_water_100.200_means <- c(random_clay_loam_avail_water_100.200_means, random_clay_loam_avail_water_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points on the Baja polygon
ggplot()+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#plotting the randomly selected points just on the cropped polygon
ggplot()+
  geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_20, color ="blue")

#storing the real means
all_known_clay_loam_avail_water_0.5_mean <- mean(all_known_pop_soils$clay_loam_avail_water_0.5)
all_known_clay_loam_avail_water_100.200_mean <- mean(all_known_pop_soils$clay_loam_avail_water_100.200)


#for clay loam available water 0-5 cm

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_loam_avail_water_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_loam_avail_water_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean clay and loam available water 0-5 cm of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_clay_loam_avail_water_0.5_means <- na.omit(random_clay_loam_avail_water_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_clay_loam_avail_water_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_loam_avail_water_0.5_means[i] < all_known_clay_loam_avail_water_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
clay_loam_avail_water_0.5_random_p.value <- (total / length(random_clay_loam_avail_water_0.5_means)) #the proportion of random ANNs that are less than our ANN


#for clay loam available water 100-200 cm

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_loam_avail_water_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_loam_avail_water_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean clay and loam available water 100-200 cm of Random Populations vs.Known Populations (n=20)")+
  theme_classic()

random_clay_loam_avail_water_100.200_means <- na.omit(random_clay_loam_avail_water_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty value
for (i in 1:length(random_clay_loam_avail_water_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_loam_avail_water_100.200_means[i] > all_known_clay_loam_avail_water_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
clay_loam_avail_water_100.200_random_p.value <- (total / length(random_clay_loam_avail_water_100.200_means)) #the proportion of random ANNs that are less than our ANN



#Heat Map 


p_value_random <- c(clay_0.5_random_p.value, clay_100.200_random_p.value, 
                    silt_0.5_random_p.value,
                    silt_100.200_random_p.value,
                    sand_0.5_random_p.value,
                    sand_100.200_random_p.value,
                    ph_0.5_random_p.value,
                    ph_100.200_random_p.value,
                    soc_0.5_random_p.value,
                    soc_100.200_random_p.value,
                    vol_wat_10kpa_0.5_random_p.value,
                    vol_wat_10kpa_100.200_random_p.value,
                    vol_wat_33kpa_0.5_random_p.value,
                    vol_wat_33kpa_100.200_random_p.value,
                    vol_wat_1500kpa_0.5_random_p.value,
                    vol_wat_1500kpa_100.200_random_p.value,
                    nitrogen_0.5_random_p.value,
                    nitrogen_100.200_random_p.value,
                    sandy_avail_water_0.5_random_p.value,
                    sandy_avail_water_100.200_random_p.value,
                    clay_loam_avail_water_0.5_random_p.value,
                    clay_loam_avail_water_100.200_random_p.value
)

# Bonferroni correcting for multiple testing
p_bonf_corrected <- p.adjust(p_value_random, method = "bonferroni")
p_bonf_corrected

#creating empty dataframe for inputting the function into
random_pop.df <- data.frame("Shape.Size" = rep(c("Clay 0-5 cm", "Clay 100-200 cm", "Silt 0-5 cm", "Silt 100-200 cm", "Sand 0-5 cm", "Sand 100-200 cm",
                                                      "pH 0-5 cm", "pH 100-200 cm", "Soil Organic Carbon 0-5 cm", "Soil Organic Carbon 100-200 cm", 
                                                 "Volume of water content -10 kPa 0-5 cm", "Volume of water content -10 kPa 100-200 cm",
                                                 "Volume of water content -33 kPa 0-5 cm", "Volume of water content -33 kPa 100-200 cm",
                                                       "Volume of water content -1500 kPa 0-5 cm", "Volume of water content -1500 kPa 100-200 cm", 
                                                      "Nitrogen 0-5 cm", "Nitrogen 100-200 cm", "Sand Available Water 0-5 cm", "Sand Available Water 100-200 cm",
                                                 "Clay/Loam Available Water 0-5 cm", "Clay/Loam Available Water 100-200 cm")),
                                 "P_Value" = p_bonf_corrected,
                                 "Significance" = c(rep(NA, 22)))   #ifelse(p_values < 0.05, "Y", "N")

#creating the significance column for the p-values
random_pop.df <- random_pop.df %>%
  mutate(Significance = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                  p_bonf_corrected >= 0.05 ~ "N"))

#labeled p-values
ggplot(aes(x = fct_reorder(Shape.Size, P_Value), y = Significance, fill = P_Value), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Association Between Soil Metrics and Population Locations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_Value < 0.5, round(P_Value, 4), NA)), col = "white") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))

