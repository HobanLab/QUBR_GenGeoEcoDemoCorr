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
library(rstatix) #to run the Games-Howell Test


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


soil_stack_LM.df <- as.data.frame(getValues(soil_stack_LM))

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
LM_soil_text_raster_250_data_pts <- extract(soil_stack_LM_soil_text, LM_fixed_field_data_processed) #extracting soil textures for each point value
LM_soil_other_raster_250_data_pts <- extract(soil_stack_LM_other, LM_fixed_field_data_processed) #extracting the other soil variables for each point value
LM_soil_extra_raster_250_data_pts <- extract(soil_stack_LM_extra, LM_fixed_field_data_processed) #extracting the extra soil variables for each point value
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed, LM_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LM point dataframe
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LM point dataframe
LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LM point dataframe


#LC
LC_soil_text_raster_250_data_pts <- extract(soil_stack_LC_soil_text, LC_fixed_field_data_processed) #extracting soil textures for each point value
LC_soil_other_raster_250_data_pts <- extract(soil_stack_LC_other, LC_fixed_field_data_processed) #extracting the other soil variables for each point value
LC_soil_extra_raster_250_data_pts <- extract(soil_stack_LC_extra, LC_fixed_field_data_processed) #extracting the extra soil variables for each point value
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed, LC_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe


#SD
SD_soil_text_raster_250_data_pts <- extract(soil_stack_SD_soil_text, SD_fixed_field_data_processed) #extracting soil textures for each point value
SD_soil_other_raster_250_data_pts <- extract(soil_stack_SD_other, SD_fixed_field_data_processed) #extracting the other soil variables for each point value
SD_soil_extra_raster_250_data_pts <- extract(soil_stack_SD_extra, SD_fixed_field_data_processed) #extracting the extra soil variables for each point value
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed, SD_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe


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
  geom_sf(data = LM_fixed_field_data_processed_trees, color = "red")



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


#combining the LM, LC, and SD tree randomly chosen tree point data into one dataframe

fixed_field_data_processed_trees_soils <- rbind(LM_fixed_field_data_processed_trees_soils, LC_fixed_field_data_processed_trees_soils) #combining the LM and LC soil and randomly chosen tree data
fixed_field_data_processed_trees_soils <- rbind(fixed_field_data_processed_trees_soils, SD_fixed_field_data_processed_trees_soils) #combining the SD tree point data to the LM and LC soil and randomly chosen tree point data


#creating a locality as factor column to be able to use Tamhane's T2 Test later
fixed_field_data_processed_trees_soils$Locality_Factor <- as.factor(fixed_field_data_processed_trees_soils$Locality)


#ANOVA comparing mean soil values between population 

##clay 0-5 cm

anova_clay_0_5 <- aov(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.0.5))+
  theme_minimal()

# checking to see if residuals are normal
hist(anova_clay_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content vs. Population")

qqnorm(anova_clay_0_5$residuals) #qqnorm plot

shapiro.test(anova_clay_0_5$residuals) #Shapiro-Wilk test, not significant, meaning residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$clay.content.0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_clay_0_5 <- tapply(fixed_field_data_processed_trees_soils$clay.content.0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_clay_0_5, na.rm = T) / min(thumb_test_clay_0_5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the levene's and rule of thumb test, the data does not meet the condition of equal variance, meaning we will use a Welch test

#Welch's ANOVA, does not assume equal variances 
oneway.test(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils, var.equal = F)

#post hoc Welch's ANOVA test: Tamhane's T2 Test

tamhaneT2Test(clay.content.0.5 ~ Locality_Factor, data = fixed_field_data_processed_trees_soils)


#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$clay.content.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$clay.content.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



##clay 100-200 

anova_clay_100_200 <- aov(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.100.200))

# checking to see if residuals are normal
hist(anova_clay_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_clay_100_200$residuals) #qqnorm plot

shapiro.test(anova_clay_100_200$residuals) #Shapiro-Wilk test, not significant, meaning residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$clay.content.100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_clay_100_200 <- tapply(fixed_field_data_processed_trees_soils$clay.content.100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_clay_100_200, na.rm = T) / min(thumb_test_clay_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shaprio test and fligner-killeen test, the data does not meet the condition of normal residuals and equal variance, meaning we will use a Kruskal-Wallis test

#kruskall wallis test
kruskal.test(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$clay.content.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$clay.content.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#silt 0-5

anova_silt_0_5 <- aov(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.0.5))

# checking to see if residuals are normal
hist(anova_silt_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_silt_0_5$residuals) #qqnorm plot

shapiro.test(anova_silt_0_5$residuals) #Shapiro-Wilk test, significant, meaning residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils) #not significant so semi equal variance

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$silt.0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_silt_0_5 <- tapply(fixed_field_data_processed_trees_soils$silt.0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_silt_0_5, na.rm = T) / min(thumb_test_silt_0_5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shaprio test, the data does not meet the condition of normal residuals, meaning we will use a Kruskal-Wallis test

#kruskall wallis test
kruskal.test(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$silt.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$silt.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



##silt 100-200

anova_silt_100_200 <- aov(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.100.200))

# checking to see if residuals are normal
hist(anova_silt_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_silt_100_200$residuals) #qqnorm plot

shapiro.test(anova_silt_100_200$residuals) #Shapiro-Wilk test, not significant, meaning residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$silt.100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_clay_100_200 <- tapply(fixed_field_data_processed_trees_soils$silt.100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_clay_100_200, na.rm = T) / min(thumb_test_clay_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shaprio test and fligner-killeen test, the data appears to meet all the conditions and we can use a regular anova and XX Test

#ANOVA test 
anova(anova_silt_100_200)

#post-hoc pairwise t tests

pairwise.t.test(fixed_field_data_processed_trees_soils$silt.100.200, fixed_field_data_processed_trees_soils$Locality, p.adj.method = "bonf")


#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$silt.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$silt.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



##sand  0-5 

anova_sand_0_5 <- aov(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.0.5))

# checking to see if residuals are normal
hist(anova_sand_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_sand_0_5$residuals) #qqnorm plot

shapiro.test(anova_sand_0_5$residuals) #Shapiro-Wilk test, not significant, meaning residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$sand.0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_sand_0_5 <- tapply(fixed_field_data_processed_trees_soils$sand.0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_sand_0_5, na.rm = T) / min(thumb_test_sand_0_5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shaprio test and fligner-killeen test, the data appears to meet all the conditions and we can use a regular anova and XX Test

#ANOVA test 
anova(anova_sand_0_5)

#post-hoc pairwise t tests

pairwise.t.test(fixed_field_data_processed_trees_soils$sand.0.5, fixed_field_data_processed_trees_soils$Locality, p.adj.method = "bonf")

#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$sand.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$sand.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method




## sand 100-200

anova_sand_100_200 <- aov(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.100.200))

# checking to see if residuals are normal
hist(anova_sand_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_sand_100_200$residuals) #qqnorm plot

shapiro.test(anova_sand_100_200$residuals) #Shapiro-Wilk test, not significant, meaning residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$sand.100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_sand_100_200 <- tapply(fixed_field_data_processed_trees_soils$sand.100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_sand_100_200, na.rm = T) / min(thumb_test_sand_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the fligner-killeen test, the data appears to meet he equal variance condition, so we will use the welch's anova and tamhane's t2 post hoc test

#Welch's ANOVA, does not assume equal variances 
oneway.test(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils, var.equal = F)

#post hoc Welch's ANOVA test: Tamhane's T2 Test

tamhaneT2Test(sand.100.200 ~ Locality_Factor, data = fixed_field_data_processed_trees_soils)

#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$sand.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$sand.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method




## ph 0-5

anova_ph_0_5 <- aov(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_0.5))

# checking to see if residuals are normal
hist(anova_ph_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_ph_0_5$residuals) #qqnorm plot

shapiro.test(anova_ph_0_5$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$ph_0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_ph_100_200 <- tapply(fixed_field_data_processed_trees_soils$ph_0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_ph_100_200, na.rm = T) / min(thumb_test_ph_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions, so we will use the kruskal wallis and  wilcox post hoc test

#kruskall wallis test
kruskal.test(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$ph_0.5 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$ph_0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



##ph 100-200

anova_ph_100_200 <- aov(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_100.200))

# checking to see if residuals are normal
hist(anova_ph_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_ph_100_200$residuals) #qqnorm plot

shapiro.test(anova_ph_100_200$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$ph_100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_ph_100_200 <- tapply(fixed_field_data_processed_trees_soils$ph_0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_ph_100_200, na.rm = T) / min(thumb_test_ph_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions, so we will use the kruskal wallis and  wilcox post hoc test

#kruskall wallis test
kruskal.test(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$ph_100.200 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$ph_100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



##soil organic carbon 0-5

anova_soc_0_5 <- aov(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.0.5))

# checking to see if residuals are normal
hist(anova_soc_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

qqnorm(anova_soc_0_5$residuals) #qqnorm plot

shapiro.test(anova_soc_0_5$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$SOC.0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_soc_0_5 <- tapply(fixed_field_data_processed_trees_soils$ph_0.5, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_soc_0_5, na.rm = T) / min(thumb_test_soc_0_5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions, so we will use the kruskal wallis and  wilcox post hoc test

#kruskall wallis test
kruskal.test(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$SOC.0.5 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$SOC.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#soil organic carbon 100-200

anova_soc_100_200 <- aov(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.100.200))

# checking to see if residuals are normal
hist(anova_soc_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_soc_100_200$residuals) #qqnorm plot

shapiro.test(anova_soc_100_200$residuals) #Shapiro-Wilk test, is not significant, meaning the residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$SOC.100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_soc_100_200 <- tapply(fixed_field_data_processed_trees_soils$SOC.100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_soc_100_200, na.rm = T) / min(thumb_test_soc_100_200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data meets the conditions and we can use a regular ANOVA and pairwise t test

#ANOVA test 
anova(anova_soc_100_200)

#post-hoc pairwise t tests

pairwise.t.test(fixed_field_data_processed_trees_soils$SOC.100.200, fixed_field_data_processed_trees_soils$Locality, p.adj.method = "bonf")

#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$SOC.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$SOC.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



#volume of water content at -10 kpa 0-5

anova_vol_water_10_0.5 <- aov(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_0.5))

# checking to see if residuals are normal
hist(anova_vol_water_10_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_vol_water_10_0.5$residuals) #qqnorm plot

shapiro.test(anova_vol_water_10_0.5$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$vol_water_.10_0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_vol_water_10_0.5 <- tapply(fixed_field_data_processed_trees_soils$SOC.100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_vol_water_10_0.5, na.rm = T) / min(thumb_test_vol_water_10_0.5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions of normal residuals and so we have to use the kruskal wallis test

#kruskall wallis test
kruskal.test(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.10_0.5 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.10_0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



#volume of water content at -10 kpa 100-200

anova_vol_water_10_100.200 <- aov(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_100.200))

# checking to see if residuals are normal
hist(anova_vol_water_10_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_vol_water_10_100.200$residuals) #qqnorm plot

shapiro.test(anova_vol_water_10_100.200$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are not normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$vol_water_.10_100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_vol_water_10_100.200 <- tapply(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_vol_water_10_100.200, na.rm = T) / min(thumb_test_vol_water_10_100.200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions of normal residuals and so we have to use the kruskal wallis test

#kruskall wallis test
kruskal.test(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.10_100.200 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method



#volume of water content at -1500 kpa 0-5

anova_vol_water_1500_0.5 <- aov(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500kPa_0.5))

# checking to see if residuals are normal
hist(anova_vol_water_1500_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_vol_water_1500_0.5$residuals) #qqnorm plot

shapiro.test(anova_vol_water_1500_0.5$residuals) #Shapiro-Wilk test, is not significant, meaning the residuals are normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$vol_water_.1500kPa_0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_vol_water_1500_0.5 <- tapply(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_vol_water_1500_0.5, na.rm = T) / min(thumb_test_vol_water_1500_0.5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data meets the condition of normal residuals but not equal variance so we will use a welch's anova and tamhanes t2 posthoc test

#based on the levene's and rule of thumb test, the data does not meet the condition of equal variance, meaning we will use a Welch test

#Welch's ANOVA, does not assume equal variances 
oneway.test(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils, var.equal = F)

#post hoc Welch's ANOVA test: Tamhane's T2 Test

tamhaneT2Test(vol_water_.1500kPa_0.5 ~ Locality_Factor, data = fixed_field_data_processed_trees_soils)


#Despite that the data did not meet the condition of equal variance we will also add results of the kruskal-wallis

#kruskall wallis test
kruskal.test(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.1500kPa_0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.1500kPa_0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method




#volume of water content at -1500 kpa 100-200

anova_vol_water_1500_100.200 <- aov(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500_100.200))

# checking to see if residuals are normal
hist(anova_vol_water_1500_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_vol_water_1500_100.200$residuals) #qqnorm plot

shapiro.test(anova_vol_water_1500_100.200$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are NOT normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$vol_water_.1500_100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_vol_water_1500_100.200 <- tapply(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_vol_water_1500_100.200, na.rm = T) / min(thumb_test_vol_water_1500_100.200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions of normal so we will use kruskal wallice

#based on the levene's and rule of thumb test, the data does not meet the condition of equal variance, meaning we will use a Welch test

#kruskall wallis test
kruskal.test(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.1500_100.200 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$vol_water_.1500_100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


#nitrogen 05-

anova_nitrogen_0.5 <- aov(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.0.5))

# checking to see if residuals are normal
hist(anova_nitrogen_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_nitrogen_0.5$residuals) #qqnorm plot

shapiro.test(anova_nitrogen_0.5$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are NOT normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$nitrogen.0.5 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_nitrogen_0.5 <- tapply(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_nitrogen_0.5, na.rm = T) / min(thumb_test_nitrogen_0.5, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions of normal so we will use kruskal wallice

#based on the levene's and rule of thumb test, the data does not meet the condition of equal variance, meaning we will use a Welch test

#kruskall wallis test
kruskal.test(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$nitrogen.0.5 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$nitrogen.0.5, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


# nitrogen 100-200


anova_nitrogen_100.200 <- aov(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.100.200))

# checking to see if residuals are normal
hist(anova_nitrogen_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Oranic Carbon at 100-200 cm vs. Population")

qqnorm(anova_nitrogen_100.200$residuals) #qqnorm plot

shapiro.test(anova_nitrogen_100.200$residuals) #Shapiro-Wilk test, is significant, meaning the residuals are NOT normal

# checking equal variances with levene's test and rule of thumb

#Fligner-Killeen, more useful when data is not normal or there are outliers 
fligner.test(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#bartlett's test for equal variances when data is normal, which in this case it is
bartlett.test(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#levene's test, not super robust to strong differences to normality
leveneTest(fixed_field_data_processed_trees_soils$nitrogen.100.200 ~ fixed_field_data_processed_trees_soils$Locality)

#rule of thumb test
thumb_test_nitrogen_100.200 <- tapply(fixed_field_data_processed_trees_soils$vol_water_.10_100.200, fixed_field_data_processed_trees_soils$Locality, sd)
max(thumb_test_nitrogen_100.200, na.rm = T) / min(thumb_test_nitrogen_100.200, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#based on the shapiro test and fligner-killeen test, the data does not meet the conditions of normal so we will use kruskal wallice

#based on the levene's and rule of thumb test, the data does not meet the condition of equal variance, meaning we will use a Welch test

#kruskall wallis test
kruskal.test(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)

#post-hoc Wilcoxon rank sum tests
pairwise.wilcox.test(fixed_field_data_processed_trees_soils$nitrogen.100.200 , fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(fixed_field_data_processed_trees_soils$nitrogen.100.200, fixed_field_data_processed_trees_soils$Locality,
                     p.adjust.method = "fdr") #p value adjusted using false discovery rate method


### Comparing the soil vs. size values ###


#FOR THIS WE ARE GOING TO SHUFFLE THE SHAPE/SIZE VALUES BETWEEN ALL OF THE POINTS IN 
#EACH PERMUTATION AND CALCULATE THE SLOPE OF THE SHAPE AND SOIL VALUE, WE THEN COMPARE THE RANDOMIZED SLOPE 
#TO OUR REAL SLOPE AND EXTRACT THE P-VALUE


### LM

## sca

# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_clay_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LM_sca_clay_0_5_lm_sum <- summary(LM_sca_clay_0_5_lm) #extracting the linear regression information
  LM_sca_clay_0.5_slopes <- c(LM_sca_clay_0.5_slopes, LM_sca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}



#extracting the slope of our points
LM_sca_clay_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LM_sca_clay_0_5_lm_real_sum <- summary(LM_sca_clay_0_5_lm_real) #extract the summary 
LM_sca_clay_0_5_lm_real_slope <- LM_sca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope



#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_soils, (aes(x=clay.content.0.5, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Clay Content 0-5")+
  ylab("Canopy Short")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_sca_clay_0_5_lm_real, aes(x= LM_sca_clay_0_5_lm_real$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Short vs. Clay Content 0-5")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_sca_clay_0_5_lm_real, aes(sample = LM_sca_clay_0_5_lm_real$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_sca_clay_0_5_lm_real$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_sca_clay_0_5_lm_real, aes(x = LM_sca_clay_0_5_lm_real$fitted.values, y = LM_sca_clay_0_5_lm_real$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for Canopy Short vs. Clay Content 0-5")

#Slope Test visible in summary of the lm
summary(LM_sca_clay_0_5_lm_real)



#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_clay_0.5_slopes[i] < LM_sca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN
  
# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_clay_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LM_sca_clay_100_200_lm_sum <- summary(LM_sca_clay_100_200_lm) #extracting the linear regression information
  LM_sca_clay_100_200_slopes <- c(LM_sca_clay_100_200_slopes, LM_sca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_clay_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LM_sca_clay_100_200_lm_real_sum <- summary(LM_sca_clay_100_200_lm_real) #extract the summary 
LM_sca_clay_100_200_lm_real_slope <- LM_sca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_soils, (aes(x=clay.content.100.200, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Clay Content 100-200")+
  ylab("Canopy Short")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_sca_clay_100_200_lm_real, aes(x= LM_sca_clay_100_200_lm_real$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Short vs. Clay Content 100-200")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_sca_clay_100_200_lm_real, aes(sample = LM_sca_clay_100_200_lm_real$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_sca_clay_100_200_lm_real$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_sca_clay_100_200_lm_real, aes(x = LM_sca_clay_100_200_lm_real$fitted.values, y = LM_sca_clay_0_5_lm_real$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for Canopy Short vs. Clay Content 100-200")


#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_clay_100_200_slopes[i] < LM_sca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_silt_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LM_sca_silt_0_5_lm_sum <- summary(LM_sca_silt_0_5_lm) #extracting the linear regression information
  LM_sca_silt_0.5_slopes <- c(LM_sca_silt_0.5_slopes, LM_sca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_silt_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LM_sca_silt_0_5_lm_real_sum <- summary(LM_sca_silt_0_5_lm_real) #extract the summary 
LM_sca_silt_0_5_lm_real_slope <- LM_sca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_soils, (aes(x=silt.0.5, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Silt Content 0-5")+
  ylab("Canopy Short")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_sca_silt_0_5_lm_real, aes(x= LM_sca_silt_0_5_lm_real$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Short vs. Silt Content 0-5")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_sca_silt_0_5_lm_real, aes(sample = LM_sca_silt_0_5_lm_real$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_sca_silt_0_5_lm_real$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_sca_silt_0_5_lm_real, aes(x = LM_sca_silt_0_5_lm_real$fitted.values, y = LM_sca_clay_0_5_lm_real$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for Canopy Short vs. Silt Content 0-5")


#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_silt_0.5_slopes[i] < LM_sca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_silt_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LM_sca_silt_100_200_lm_sum <- summary(LM_sca_silt_100_200_lm) #extracting the linear regression information
  LM_sca_silt_100_200_slopes <- c(LM_sca_silt_100_200_slopes, LM_sca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_silt_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LM_sca_silt_100_200_lm_real_sum <- summary(LM_sca_silt_100_200_lm_real) #extract the summary 
LM_sca_silt_100_200_lm_real_slope <- LM_sca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_soils, (aes(x=silt.100.200, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Silt Content 100-200")+
  ylab("Canopy Short")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_sca_silt_100_200_lm_real, aes(x= LM_sca_silt_100_200_lm_real$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Short vs. Silt Content 100-200")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_sca_silt_100_200_lm_real, aes(sample = LM_sca_silt_100_200_lm_real$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_sca_silt_100_200_lm_real$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_sca_silt_100_200_lm_real, aes(x = LM_sca_silt_100_200_lm_real$fitted.values, y = LM_sca_clay_0_5_lm_real$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for Canopy Short vs. Silt Content 100-200")


#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_silt_100_200_slopes[i] < LM_sca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_sand_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LM_sca_sand_0_5_lm_sum <- summary(LM_sca_sand_0_5_lm) #extracting the linear regression information
  LM_sca_sand_0.5_slopes <- c(LM_sca_sand_0.5_slopes, LM_sca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_sand_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LM_sca_sand_0_5_lm_real_sum <- summary(LM_sca_sand_0_5_lm_real) #extract the summary 
LM_sca_sand_0_5_lm_real_slope <- LM_sca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_soils, (aes(x=sand.0.5, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Sand Content 0-5")+
  ylab("Canopy Short")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_sca_sand_0_5_lm_real, aes(x= LM_sca_sand_0_5_lm_real$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Short vs. Sand Content 0-5")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_sca_sand_0_5_lm_real, aes(sample = LM_sca_sand_0_5_lm_real$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_sca_sand_0_5_lm_real$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_sca_sand_0_5_lm_real, aes(x = LM_sca_sand_0_5_lm_real$fitted.values, y = LM_sca_clay_0_5_lm_real$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for Canopy Short vs. Sand Content 0-5")


#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_sand_0.5_slopes[i] > LM_sca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_sand_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LM_sca_sand_100_200_lm_sum <- summary(LM_sca_sand_100_200_lm) #extracting the linear regression information
  LM_sca_sand_100_200_slopes <- c(LM_sca_sand_100_200_slopes, LM_sca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_sand_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LM_sca_sand_100_200_lm_real_sum <- summary(LM_sca_sand_100_200_lm_real) #extract the summary 
LM_sca_sand_100_200_lm_real_slope <- LM_sca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_sand_100_200_slopes[i] > LM_sca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_ph_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LM_sca_ph_0_5_lm_sum <- summary(LM_sca_ph_0_5_lm) #extracting the linear regression information
  LM_sca_ph_0.5_slopes <- c(LM_sca_ph_0.5_slopes, LM_sca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_ph_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LM_sca_ph_0_5_lm_real_sum <- summary(LM_sca_ph_0_5_lm_real) #extract the summary 
LM_sca_ph_0_5_lm_real_slope <- LM_sca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_ph_0.5_slopes[i] > LM_sca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_ph_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LM_sca_ph_100_200_lm_sum <- summary(LM_sca_ph_100_200_lm) #extracting the linear regression information
  LM_sca_ph_100_200_slopes <- c(LM_sca_ph_100_200_slopes, LM_sca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_ph_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LM_sca_ph_100_200_lm_real_sum <- summary(LM_sca_ph_100_200_lm_real) #extract the summary 
LM_sca_ph_100_200_lm_real_slope <- LM_sca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_ph_100_200_slopes[i] > LM_sca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_soc_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LM_sca_soc_0_5_lm_sum <- summary(LM_sca_soc_0_5_lm) #extracting the linear regression information
  LM_sca_soc_0.5_slopes <- c(LM_sca_soc_0.5_slopes, LM_sca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_soc_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LM_sca_soc_0_5_lm_real_sum <- summary(LM_sca_soc_0_5_lm_real) #extract the summary 
LM_sca_soc_0_5_lm_real_slope <- LM_sca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_soc_0.5_slopes[i] < LM_sca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_soc_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LM_sca_soc_100_200_lm_sum <- summary(LM_sca_soc_100_200_lm) #extracting the linear regression information
  LM_sca_soc_100_200_slopes <- c(LM_sca_soc_100_200_slopes, LM_sca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_soc_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LM_sca_soc_100_200_lm_real_sum <- summary(LM_sca_soc_100_200_lm_real) #extract the summary 
LM_sca_soc_100_200_lm_real_slope <- LM_sca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_soc_100_200_slopes[i] < LM_sca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_vol_10_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LM_sca_vol_10_0_5_lm_sum <- summary(LM_sca_vol_10_0_5_lm) #extracting the linear regression information
  LM_sca_vol_10_0.5_slopes <- c(LM_sca_vol_10_0.5_slopes, LM_sca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_vol_10_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LM_sca_vol_10_0_5_lm_real_sum <- summary(LM_sca_vol_10_0_5_lm_real) #extract the summary 
LM_sca_vol_10_0_5_lm_real_slope <- LM_sca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_vol_10_0.5_slopes[i] > LM_sca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_vol_10_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LM_sca_vol_10_100_200_lm_sum <- summary(LM_sca_vol_10_100_200_lm) #extracting the linear regression information
  LM_sca_vol_10_100_200_slopes <- c(LM_sca_vol_10_100_200_slopes, LM_sca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_vol_10_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LM_sca_vol_10_100_200_lm_real_sum <- summary(LM_sca_vol_10_100_200_lm_real) #extract the summary 
LM_sca_vol_10_100_200_lm_real_slope <- LM_sca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_vol_10_100_200_slopes[i] < LM_sca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_vol_1500_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LM_sca_vol_1500_0_5_lm_sum <- summary(LM_sca_vol_1500_0_5_lm) #extracting the linear regression information
  LM_sca_vol_1500_0.5_slopes <- c(LM_sca_vol_1500_0.5_slopes, LM_sca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_vol_1500_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LM_sca_vol_1500_0_5_lm_real_sum <- summary(LM_sca_vol_1500_0_5_lm_real) #extract the summary 
LM_sca_vol_1500_0_5_lm_real_slope <- LM_sca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_vol_1500_0.5_slopes[i] > LM_sca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_vol_1500_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LM_sca_vol_1500_100_200_lm_sum <- summary(LM_sca_vol_1500_100_200_lm) #extracting the linear regression information
  LM_sca_vol_1500_100_200_slopes <- c(LM_sca_vol_1500_100_200_slopes, LM_sca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_vol_1500_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LM_sca_vol_1500_100_200_lm_real_sum <- summary(LM_sca_vol_1500_100_200_lm_real) #extract the summary 
LM_sca_vol_1500_100_200_lm_real_slope <- LM_sca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_vol_1500_100_200_slopes[i] < LM_sca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_nitrogen_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LM_sca_nitrogen_0_5_lm_sum <- summary(LM_sca_nitrogen_0_5_lm) #extracting the linear regression information
  LM_sca_nitrogen_0.5_slopes <- c(LM_sca_nitrogen_0.5_slopes, LM_sca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_nitrogen_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LM_sca_nitrogen_0_5_lm_real_sum <- summary(LM_sca_nitrogen_0_5_lm_real) #extract the summary 
LM_sca_nitrogen_0_5_lm_real_slope <- LM_sca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_nitrogen_0.5_slopes[i] > LM_sca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_sca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LM_sca_nitrogen_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LM_sca_nitrogen_100_200_lm_sum <- summary(LM_sca_nitrogen_100_200_lm) #extracting the linear regression information
  LM_sca_nitrogen_100_200_slopes <- c(LM_sca_nitrogen_100_200_slopes, LM_sca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LM_sca_nitrogen_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_short~LM_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LM_sca_nitrogen_100_200_lm_real_sum <- summary(LM_sca_nitrogen_100_200_lm_real) #extract the summary 
LM_sca_nitrogen_100_200_lm_real_slope <- LM_sca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_sca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_sca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_sca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_sca_nitrogen_100_200_slopes[i] > LM_sca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_sca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# LCA



# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_clay_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LM_lca_clay_0_5_lm_sum <- summary(LM_lca_clay_0_5_lm) #extracting the linear regression information
  LM_lca_clay_0.5_slopes <- c(LM_lca_clay_0.5_slopes, LM_lca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_clay_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LM_lca_clay_0_5_lm_real_sum <- summary(LM_lca_clay_0_5_lm_real) #extract the summary 
LM_lca_clay_0_5_lm_real_slope <- LM_lca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_clay_0.5_slopes[i] < LM_lca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_clay_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LM_lca_clay_100_200_lm_sum <- summary(LM_lca_clay_100_200_lm) #extracting the linear regression information
  LM_lca_clay_100_200_slopes <- c(LM_lca_clay_100_200_slopes, LM_lca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_clay_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LM_lca_clay_100_200_lm_real_sum <- summary(LM_lca_clay_100_200_lm_real) #extract the summary 
LM_lca_clay_100_200_lm_real_slope <- LM_lca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_clay_100_200_slopes[i] < LM_lca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_silt_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LM_lca_silt_0_5_lm_sum <- summary(LM_lca_silt_0_5_lm) #extracting the linear regression information
  LM_lca_silt_0.5_slopes <- c(LM_lca_silt_0.5_slopes, LM_lca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_silt_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LM_lca_silt_0_5_lm_real_sum <- summary(LM_lca_silt_0_5_lm_real) #extract the summary 
LM_lca_silt_0_5_lm_real_slope <- LM_lca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_silt_0.5_slopes[i] < LM_lca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_silt_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LM_lca_silt_100_200_lm_sum <- summary(LM_lca_silt_100_200_lm) #extracting the linear regression information
  LM_lca_silt_100_200_slopes <- c(LM_lca_silt_100_200_slopes, LM_lca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_silt_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LM_lca_silt_100_200_lm_real_sum <- summary(LM_lca_silt_100_200_lm_real) #extract the summary 
LM_lca_silt_100_200_lm_real_slope <- LM_lca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_silt_100_200_slopes[i] < LM_lca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_sand_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LM_lca_sand_0_5_lm_sum <- summary(LM_lca_sand_0_5_lm) #extracting the linear regression information
  LM_lca_sand_0.5_slopes <- c(LM_lca_sand_0.5_slopes, LM_lca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_sand_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LM_lca_sand_0_5_lm_real_sum <- summary(LM_lca_sand_0_5_lm_real) #extract the summary 
LM_lca_sand_0_5_lm_real_slope <- LM_lca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_sand_0.5_slopes[i] > LM_lca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_sand_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LM_lca_sand_100_200_lm_sum <- summary(LM_lca_sand_100_200_lm) #extracting the linear regression information
  LM_lca_sand_100_200_slopes <- c(LM_lca_sand_100_200_slopes, LM_lca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_sand_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LM_lca_sand_100_200_lm_real_sum <- summary(LM_lca_sand_100_200_lm_real) #extract the summary 
LM_lca_sand_100_200_lm_real_slope <- LM_lca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_sand_100_200_slopes[i] > LM_lca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_ph_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LM_lca_ph_0_5_lm_sum <- summary(LM_lca_ph_0_5_lm) #extracting the linear regression information
  LM_lca_ph_0.5_slopes <- c(LM_lca_ph_0.5_slopes, LM_lca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_ph_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LM_lca_ph_0_5_lm_real_sum <- summary(LM_lca_ph_0_5_lm_real) #extract the summary 
LM_lca_ph_0_5_lm_real_slope <- LM_lca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_ph_0.5_slopes[i] > LM_lca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_ph_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LM_lca_ph_100_200_lm_sum <- summary(LM_lca_ph_100_200_lm) #extracting the linear regression information
  LM_lca_ph_100_200_slopes <- c(LM_lca_ph_100_200_slopes, LM_lca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_ph_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LM_lca_ph_100_200_lm_real_sum <- summary(LM_lca_ph_100_200_lm_real) #extract the summary 
LM_lca_ph_100_200_lm_real_slope <- LM_lca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_ph_100_200_slopes[i] > LM_lca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_soc_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LM_lca_soc_0_5_lm_sum <- summary(LM_lca_soc_0_5_lm) #extracting the linear regression information
  LM_lca_soc_0.5_slopes <- c(LM_lca_soc_0.5_slopes, LM_lca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_soc_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LM_lca_soc_0_5_lm_real_sum <- summary(LM_lca_soc_0_5_lm_real) #extract the summary 
LM_lca_soc_0_5_lm_real_slope <- LM_lca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_soc_0.5_slopes[i] < LM_lca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_soc_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LM_lca_soc_100_200_lm_sum <- summary(LM_lca_soc_100_200_lm) #extracting the linear regression information
  LM_lca_soc_100_200_slopes <- c(LM_lca_soc_100_200_slopes, LM_lca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_soc_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LM_lca_soc_100_200_lm_real_sum <- summary(LM_lca_soc_100_200_lm_real) #extract the summary 
LM_lca_soc_100_200_lm_real_slope <- LM_lca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_soc_100_200_slopes[i] < LM_lca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_vol_10_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LM_lca_vol_10_0_5_lm_sum <- summary(LM_lca_vol_10_0_5_lm) #extracting the linear regression information
  LM_lca_vol_10_0.5_slopes <- c(LM_lca_vol_10_0.5_slopes, LM_lca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_vol_10_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LM_lca_vol_10_0_5_lm_real_sum <- summary(LM_lca_vol_10_0_5_lm_real) #extract the summary 
LM_lca_vol_10_0_5_lm_real_slope <- LM_lca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_vol_10_0.5_slopes[i] > LM_lca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_vol_10_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LM_lca_vol_10_100_200_lm_sum <- summary(LM_lca_vol_10_100_200_lm) #extracting the linear regression information
  LM_lca_vol_10_100_200_slopes <- c(LM_lca_vol_10_100_200_slopes, LM_lca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_vol_10_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LM_lca_vol_10_100_200_lm_real_sum <- summary(LM_lca_vol_10_100_200_lm_real) #extract the summary 
LM_lca_vol_10_100_200_lm_real_slope <- LM_lca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_vol_10_100_200_slopes[i] < LM_lca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_vol_1500_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LM_lca_vol_1500_0_5_lm_sum <- summary(LM_lca_vol_1500_0_5_lm) #extracting the linear regression information
  LM_lca_vol_1500_0.5_slopes <- c(LM_lca_vol_1500_0.5_slopes, LM_lca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_vol_1500_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LM_lca_vol_1500_0_5_lm_real_sum <- summary(LM_lca_vol_1500_0_5_lm_real) #extract the summary 
LM_lca_vol_1500_0_5_lm_real_slope <- LM_lca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_vol_1500_0.5_slopes[i] > LM_lca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_vol_1500_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LM_lca_vol_1500_100_200_lm_sum <- summary(LM_lca_vol_1500_100_200_lm) #extracting the linear regression information
  LM_lca_vol_1500_100_200_slopes <- c(LM_lca_vol_1500_100_200_slopes, LM_lca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_vol_1500_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LM_lca_vol_1500_100_200_lm_real_sum <- summary(LM_lca_vol_1500_100_200_lm_real) #extract the summary 
LM_lca_vol_1500_100_200_lm_real_slope <- LM_lca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_vol_1500_100_200_slopes[i] < LM_lca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_nitrogen_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LM_lca_nitrogen_0_5_lm_sum <- summary(LM_lca_nitrogen_0_5_lm) #extracting the linear regression information
  LM_lca_nitrogen_0.5_slopes <- c(LM_lca_nitrogen_0.5_slopes, LM_lca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_nitrogen_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LM_lca_nitrogen_0_5_lm_real_sum <- summary(LM_lca_nitrogen_0_5_lm_real) #extract the summary 
LM_lca_nitrogen_0_5_lm_real_slope <- LM_lca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_nitrogen_0.5_slopes[i] > LM_lca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_lca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LM_lca_nitrogen_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LM_lca_nitrogen_100_200_lm_sum <- summary(LM_lca_nitrogen_100_200_lm) #extracting the linear regression information
  LM_lca_nitrogen_100_200_slopes <- c(LM_lca_nitrogen_100_200_slopes, LM_lca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LM_lca_nitrogen_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_long~LM_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LM_lca_nitrogen_100_200_lm_real_sum <- summary(LM_lca_nitrogen_100_200_lm_real) #extract the summary 
LM_lca_nitrogen_100_200_lm_real_slope <- LM_lca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_lca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_lca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_lca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_lca_nitrogen_100_200_slopes[i] > LM_lca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_lca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CA


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_clay_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LM_ca_clay_0_5_lm_sum <- summary(LM_ca_clay_0_5_lm) #extracting the linear regression information
  LM_ca_clay_0.5_slopes <- c(LM_ca_clay_0.5_slopes, LM_ca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_clay_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LM_ca_clay_0_5_lm_real_sum <- summary(LM_ca_clay_0_5_lm_real) #extract the summary 
LM_ca_clay_0_5_lm_real_slope <- LM_ca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_clay_0.5_slopes[i] < LM_ca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_clay_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LM_ca_clay_100_200_lm_sum <- summary(LM_ca_clay_100_200_lm) #extracting the linear regression information
  LM_ca_clay_100_200_slopes <- c(LM_ca_clay_100_200_slopes, LM_ca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_clay_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LM_ca_clay_100_200_lm_real_sum <- summary(LM_ca_clay_100_200_lm_real) #extract the summary 
LM_ca_clay_100_200_lm_real_slope <- LM_ca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_clay_100_200_slopes[i] < LM_ca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_silt_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LM_ca_silt_0_5_lm_sum <- summary(LM_ca_silt_0_5_lm) #extracting the linear regression information
  LM_ca_silt_0.5_slopes <- c(LM_ca_silt_0.5_slopes, LM_ca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_silt_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LM_ca_silt_0_5_lm_real_sum <- summary(LM_ca_silt_0_5_lm_real) #extract the summary 
LM_ca_silt_0_5_lm_real_slope <- LM_ca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_silt_0.5_slopes[i] < LM_ca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_silt_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LM_ca_silt_100_200_lm_sum <- summary(LM_ca_silt_100_200_lm) #extracting the linear regression information
  LM_ca_silt_100_200_slopes <- c(LM_ca_silt_100_200_slopes, LM_ca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_silt_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LM_ca_silt_100_200_lm_real_sum <- summary(LM_ca_silt_100_200_lm_real) #extract the summary 
LM_ca_silt_100_200_lm_real_slope <- LM_ca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_silt_100_200_slopes[i] < LM_ca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_sand_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LM_ca_sand_0_5_lm_sum <- summary(LM_ca_sand_0_5_lm) #extracting the linear regression information
  LM_ca_sand_0.5_slopes <- c(LM_ca_sand_0.5_slopes, LM_ca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_sand_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LM_ca_sand_0_5_lm_real_sum <- summary(LM_ca_sand_0_5_lm_real) #extract the summary 
LM_ca_sand_0_5_lm_real_slope <- LM_ca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_sand_0.5_slopes[i] > LM_ca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_sand_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LM_ca_sand_100_200_lm_sum <- summary(LM_ca_sand_100_200_lm) #extracting the linear regression information
  LM_ca_sand_100_200_slopes <- c(LM_ca_sand_100_200_slopes, LM_ca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_sand_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LM_ca_sand_100_200_lm_real_sum <- summary(LM_ca_sand_100_200_lm_real) #extract the summary 
LM_ca_sand_100_200_lm_real_slope <- LM_ca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_sand_100_200_slopes[i] > LM_ca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_ph_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LM_ca_ph_0_5_lm_sum <- summary(LM_ca_ph_0_5_lm) #extracting the linear regression information
  LM_ca_ph_0.5_slopes <- c(LM_ca_ph_0.5_slopes, LM_ca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_ph_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LM_ca_ph_0_5_lm_real_sum <- summary(LM_ca_ph_0_5_lm_real) #extract the summary 
LM_ca_ph_0_5_lm_real_slope <- LM_ca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_ph_0.5_slopes[i] > LM_ca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_ph_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LM_ca_ph_100_200_lm_sum <- summary(LM_ca_ph_100_200_lm) #extracting the linear regression information
  LM_ca_ph_100_200_slopes <- c(LM_ca_ph_100_200_slopes, LM_ca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_ph_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LM_ca_ph_100_200_lm_real_sum <- summary(LM_ca_ph_100_200_lm_real) #extract the summary 
LM_ca_ph_100_200_lm_real_slope <- LM_ca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_ph_100_200_slopes[i] > LM_ca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_soc_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LM_ca_soc_0_5_lm_sum <- summary(LM_ca_soc_0_5_lm) #extracting the linear regression information
  LM_ca_soc_0.5_slopes <- c(LM_ca_soc_0.5_slopes, LM_ca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_soc_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LM_ca_soc_0_5_lm_real_sum <- summary(LM_ca_soc_0_5_lm_real) #extract the summary 
LM_ca_soc_0_5_lm_real_slope <- LM_ca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_soc_0.5_slopes[i] < LM_ca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_soc_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LM_ca_soc_100_200_lm_sum <- summary(LM_ca_soc_100_200_lm) #extracting the linear regression information
  LM_ca_soc_100_200_slopes <- c(LM_ca_soc_100_200_slopes, LM_ca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_soc_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LM_ca_soc_100_200_lm_real_sum <- summary(LM_ca_soc_100_200_lm_real) #extract the summary 
LM_ca_soc_100_200_lm_real_slope <- LM_ca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_soc_100_200_slopes[i] < LM_ca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_vol_10_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LM_ca_vol_10_0_5_lm_sum <- summary(LM_ca_vol_10_0_5_lm) #extracting the linear regression information
  LM_ca_vol_10_0.5_slopes <- c(LM_ca_vol_10_0.5_slopes, LM_ca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_vol_10_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LM_ca_vol_10_0_5_lm_real_sum <- summary(LM_ca_vol_10_0_5_lm_real) #extract the summary 
LM_ca_vol_10_0_5_lm_real_slope <- LM_ca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_vol_10_0.5_slopes[i] > LM_ca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_vol_10_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LM_ca_vol_10_100_200_lm_sum <- summary(LM_ca_vol_10_100_200_lm) #extracting the linear regression information
  LM_ca_vol_10_100_200_slopes <- c(LM_ca_vol_10_100_200_slopes, LM_ca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_vol_10_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LM_ca_vol_10_100_200_lm_real_sum <- summary(LM_ca_vol_10_100_200_lm_real) #extract the summary 
LM_ca_vol_10_100_200_lm_real_slope <- LM_ca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_vol_10_100_200_slopes[i] < LM_ca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_vol_1500_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LM_ca_vol_1500_0_5_lm_sum <- summary(LM_ca_vol_1500_0_5_lm) #extracting the linear regression information
  LM_ca_vol_1500_0.5_slopes <- c(LM_ca_vol_1500_0.5_slopes, LM_ca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_vol_1500_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LM_ca_vol_1500_0_5_lm_real_sum <- summary(LM_ca_vol_1500_0_5_lm_real) #extract the summary 
LM_ca_vol_1500_0_5_lm_real_slope <- LM_ca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_vol_1500_0.5_slopes[i] > LM_ca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_vol_1500_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LM_ca_vol_1500_100_200_lm_sum <- summary(LM_ca_vol_1500_100_200_lm) #extracting the linear regression information
  LM_ca_vol_1500_100_200_slopes <- c(LM_ca_vol_1500_100_200_slopes, LM_ca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_vol_1500_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LM_ca_vol_1500_100_200_lm_real_sum <- summary(LM_ca_vol_1500_100_200_lm_real) #extract the summary 
LM_ca_vol_1500_100_200_lm_real_slope <- LM_ca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_vol_1500_100_200_slopes[i] < LM_ca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_nitrogen_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LM_ca_nitrogen_0_5_lm_sum <- summary(LM_ca_nitrogen_0_5_lm) #extracting the linear regression information
  LM_ca_nitrogen_0.5_slopes <- c(LM_ca_nitrogen_0.5_slopes, LM_ca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_nitrogen_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LM_ca_nitrogen_0_5_lm_real_sum <- summary(LM_ca_nitrogen_0_5_lm_real) #extract the summary 
LM_ca_nitrogen_0_5_lm_real_slope <- LM_ca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_nitrogen_0.5_slopes[i] > LM_ca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_ca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LM_ca_nitrogen_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LM_ca_nitrogen_100_200_lm_sum <- summary(LM_ca_nitrogen_100_200_lm) #extracting the linear regression information
  LM_ca_nitrogen_100_200_slopes <- c(LM_ca_nitrogen_100_200_slopes, LM_ca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LM_ca_nitrogen_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Canopy_area~LM_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LM_ca_nitrogen_100_200_lm_real_sum <- summary(LM_ca_nitrogen_100_200_lm_real) #extract the summary 
LM_ca_nitrogen_100_200_lm_real_slope <- LM_ca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_ca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_ca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_ca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_ca_nitrogen_100_200_slopes[i] > LM_ca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_ca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CS


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_clay_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LM_cs_clay_0_5_lm_sum <- summary(LM_cs_clay_0_5_lm) #extracting the linear regression information
  LM_cs_clay_0.5_slopes <- c(LM_cs_clay_0.5_slopes, LM_cs_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_clay_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LM_cs_clay_0_5_lm_real_sum <- summary(LM_cs_clay_0_5_lm_real) #extract the summary 
LM_cs_clay_0_5_lm_real_slope <- LM_cs_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_clay_0.5_slopes[i] < LM_cs_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_clay_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LM_cs_clay_100_200_lm_sum <- summary(LM_cs_clay_100_200_lm) #extracting the linear regression information
  LM_cs_clay_100_200_slopes <- c(LM_cs_clay_100_200_slopes, LM_cs_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_clay_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LM_cs_clay_100_200_lm_real_sum <- summary(LM_cs_clay_100_200_lm_real) #extract the summary 
LM_cs_clay_100_200_lm_real_slope <- LM_cs_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_clay_100_200_slopes[i] < LM_cs_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_silt_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LM_cs_silt_0_5_lm_sum <- summary(LM_cs_silt_0_5_lm) #extracting the linear regression information
  LM_cs_silt_0.5_slopes <- c(LM_cs_silt_0.5_slopes, LM_cs_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_silt_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LM_cs_silt_0_5_lm_real_sum <- summary(LM_cs_silt_0_5_lm_real) #extract the summary 
LM_cs_silt_0_5_lm_real_slope <- LM_cs_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_silt_0.5_slopes[i] < LM_cs_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_silt_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LM_cs_silt_100_200_lm_sum <- summary(LM_cs_silt_100_200_lm) #extracting the linear regression information
  LM_cs_silt_100_200_slopes <- c(LM_cs_silt_100_200_slopes, LM_cs_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_silt_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LM_cs_silt_100_200_lm_real_sum <- summary(LM_cs_silt_100_200_lm_real) #extract the summary 
LM_cs_silt_100_200_lm_real_slope <- LM_cs_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_silt_100_200_slopes[i] < LM_cs_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_sand_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LM_cs_sand_0_5_lm_sum <- summary(LM_cs_sand_0_5_lm) #extracting the linear regression information
  LM_cs_sand_0.5_slopes <- c(LM_cs_sand_0.5_slopes, LM_cs_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_sand_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LM_cs_sand_0_5_lm_real_sum <- summary(LM_cs_sand_0_5_lm_real) #extract the summary 
LM_cs_sand_0_5_lm_real_slope <- LM_cs_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_sand_0.5_slopes[i] > LM_cs_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_sand_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LM_cs_sand_100_200_lm_sum <- summary(LM_cs_sand_100_200_lm) #extracting the linear regression information
  LM_cs_sand_100_200_slopes <- c(LM_cs_sand_100_200_slopes, LM_cs_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_sand_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LM_cs_sand_100_200_lm_real_sum <- summary(LM_cs_sand_100_200_lm_real) #extract the summary 
LM_cs_sand_100_200_lm_real_slope <- LM_cs_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_sand_100_200_slopes[i] > LM_cs_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_ph_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LM_cs_ph_0_5_lm_sum <- summary(LM_cs_ph_0_5_lm) #extracting the linear regression information
  LM_cs_ph_0.5_slopes <- c(LM_cs_ph_0.5_slopes, LM_cs_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_ph_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LM_cs_ph_0_5_lm_real_sum <- summary(LM_cs_ph_0_5_lm_real) #extract the summary 
LM_cs_ph_0_5_lm_real_slope <- LM_cs_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_ph_0.5_slopes[i] > LM_cs_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_ph_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LM_cs_ph_100_200_lm_sum <- summary(LM_cs_ph_100_200_lm) #extracting the linear regression information
  LM_cs_ph_100_200_slopes <- c(LM_cs_ph_100_200_slopes, LM_cs_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_ph_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LM_cs_ph_100_200_lm_real_sum <- summary(LM_cs_ph_100_200_lm_real) #extract the summary 
LM_cs_ph_100_200_lm_real_slope <- LM_cs_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_ph_100_200_slopes[i] > LM_cs_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_soc_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LM_cs_soc_0_5_lm_sum <- summary(LM_cs_soc_0_5_lm) #extracting the linear regression information
  LM_cs_soc_0.5_slopes <- c(LM_cs_soc_0.5_slopes, LM_cs_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_soc_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LM_cs_soc_0_5_lm_real_sum <- summary(LM_cs_soc_0_5_lm_real) #extract the summary 
LM_cs_soc_0_5_lm_real_slope <- LM_cs_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_soc_0.5_slopes[i] < LM_cs_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_soc_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LM_cs_soc_100_200_lm_sum <- summary(LM_cs_soc_100_200_lm) #extracting the linear regression information
  LM_cs_soc_100_200_slopes <- c(LM_cs_soc_100_200_slopes, LM_cs_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_soc_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LM_cs_soc_100_200_lm_real_sum <- summary(LM_cs_soc_100_200_lm_real) #extract the summary 
LM_cs_soc_100_200_lm_real_slope <- LM_cs_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_soc_100_200_slopes[i] < LM_cs_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_vol_10_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LM_cs_vol_10_0_5_lm_sum <- summary(LM_cs_vol_10_0_5_lm) #extracting the linear regression information
  LM_cs_vol_10_0.5_slopes <- c(LM_cs_vol_10_0.5_slopes, LM_cs_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_vol_10_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LM_cs_vol_10_0_5_lm_real_sum <- summary(LM_cs_vol_10_0_5_lm_real) #extract the summary 
LM_cs_vol_10_0_5_lm_real_slope <- LM_cs_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_vol_10_0.5_slopes[i] > LM_cs_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_vol_10_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LM_cs_vol_10_100_200_lm_sum <- summary(LM_cs_vol_10_100_200_lm) #extracting the linear regression information
  LM_cs_vol_10_100_200_slopes <- c(LM_cs_vol_10_100_200_slopes, LM_cs_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_vol_10_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LM_cs_vol_10_100_200_lm_real_sum <- summary(LM_cs_vol_10_100_200_lm_real) #extract the summary 
LM_cs_vol_10_100_200_lm_real_slope <- LM_cs_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_vol_10_100_200_slopes[i] < LM_cs_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_vol_1500_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LM_cs_vol_1500_0_5_lm_sum <- summary(LM_cs_vol_1500_0_5_lm) #extracting the linear regression information
  LM_cs_vol_1500_0.5_slopes <- c(LM_cs_vol_1500_0.5_slopes, LM_cs_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_vol_1500_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LM_cs_vol_1500_0_5_lm_real_sum <- summary(LM_cs_vol_1500_0_5_lm_real) #extract the summary 
LM_cs_vol_1500_0_5_lm_real_slope <- LM_cs_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_vol_1500_0.5_slopes[i] > LM_cs_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_vol_1500_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LM_cs_vol_1500_100_200_lm_sum <- summary(LM_cs_vol_1500_100_200_lm) #extracting the linear regression information
  LM_cs_vol_1500_100_200_slopes <- c(LM_cs_vol_1500_100_200_slopes, LM_cs_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_vol_1500_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LM_cs_vol_1500_100_200_lm_real_sum <- summary(LM_cs_vol_1500_100_200_lm_real) #extract the summary 
LM_cs_vol_1500_100_200_lm_real_slope <- LM_cs_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_vol_1500_100_200_slopes[i] < LM_cs_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_nitrogen_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LM_cs_nitrogen_0_5_lm_sum <- summary(LM_cs_nitrogen_0_5_lm) #extracting the linear regression information
  LM_cs_nitrogen_0.5_slopes <- c(LM_cs_nitrogen_0.5_slopes, LM_cs_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_nitrogen_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LM_cs_nitrogen_0_5_lm_real_sum <- summary(LM_cs_nitrogen_0_5_lm_real) #extract the summary 
LM_cs_nitrogen_0_5_lm_real_slope <- LM_cs_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_nitrogen_0.5_slopes[i] > LM_cs_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_cs_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LM_cs_nitrogen_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LM_cs_nitrogen_100_200_lm_sum <- summary(LM_cs_nitrogen_100_200_lm) #extracting the linear regression information
  LM_cs_nitrogen_100_200_slopes <- c(LM_cs_nitrogen_100_200_slopes, LM_cs_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LM_cs_nitrogen_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$Crown_spread~LM_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LM_cs_nitrogen_100_200_lm_real_sum <- summary(LM_cs_nitrogen_100_200_lm_real) #extract the summary 
LM_cs_nitrogen_100_200_lm_real_slope <- LM_cs_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_cs_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_cs_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_cs_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_cs_nitrogen_100_200_slopes[i] > LM_cs_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_cs_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



#DBH_ag


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_clay_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LM_dbh_clay_0_5_lm_sum <- summary(LM_dbh_clay_0_5_lm) #extracting the linear regression information
  LM_dbh_clay_0.5_slopes <- c(LM_dbh_clay_0.5_slopes, LM_dbh_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_clay_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LM_dbh_clay_0_5_lm_real_sum <- summary(LM_dbh_clay_0_5_lm_real) #extract the summary 
LM_dbh_clay_0_5_lm_real_slope <- LM_dbh_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_clay_0.5_slopes[i] < LM_dbh_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_clay_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LM_dbh_clay_100_200_lm_sum <- summary(LM_dbh_clay_100_200_lm) #extracting the linear regression information
  LM_dbh_clay_100_200_slopes <- c(LM_dbh_clay_100_200_slopes, LM_dbh_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_clay_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LM_dbh_clay_100_200_lm_real_sum <- summary(LM_dbh_clay_100_200_lm_real) #extract the summary 
LM_dbh_clay_100_200_lm_real_slope <- LM_dbh_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_clay_100_200_slopes[i] < LM_dbh_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_silt_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LM_dbh_silt_0_5_lm_sum <- summary(LM_dbh_silt_0_5_lm) #extracting the linear regression information
  LM_dbh_silt_0.5_slopes <- c(LM_dbh_silt_0.5_slopes, LM_dbh_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_silt_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LM_dbh_silt_0_5_lm_real_sum <- summary(LM_dbh_silt_0_5_lm_real) #extract the summary 
LM_dbh_silt_0_5_lm_real_slope <- LM_dbh_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_silt_0.5_slopes[i] < LM_dbh_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_silt_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LM_dbh_silt_100_200_lm_sum <- summary(LM_dbh_silt_100_200_lm) #extracting the linear regression information
  LM_dbh_silt_100_200_slopes <- c(LM_dbh_silt_100_200_slopes, LM_dbh_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_silt_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LM_dbh_silt_100_200_lm_real_sum <- summary(LM_dbh_silt_100_200_lm_real) #extract the summary 
LM_dbh_silt_100_200_lm_real_slope <- LM_dbh_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_silt_100_200_slopes[i] < LM_dbh_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_sand_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LM_dbh_sand_0_5_lm_sum <- summary(LM_dbh_sand_0_5_lm) #extracting the linear regression information
  LM_dbh_sand_0.5_slopes <- c(LM_dbh_sand_0.5_slopes, LM_dbh_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_sand_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LM_dbh_sand_0_5_lm_real_sum <- summary(LM_dbh_sand_0_5_lm_real) #extract the summary 
LM_dbh_sand_0_5_lm_real_slope <- LM_dbh_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_sand_0.5_slopes[i] > LM_dbh_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_sand_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LM_dbh_sand_100_200_lm_sum <- summary(LM_dbh_sand_100_200_lm) #extracting the linear regression information
  LM_dbh_sand_100_200_slopes <- c(LM_dbh_sand_100_200_slopes, LM_dbh_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_sand_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LM_dbh_sand_100_200_lm_real_sum <- summary(LM_dbh_sand_100_200_lm_real) #extract the summary 
LM_dbh_sand_100_200_lm_real_slope <- LM_dbh_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_sand_100_200_slopes[i] > LM_dbh_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_ph_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LM_dbh_ph_0_5_lm_sum <- summary(LM_dbh_ph_0_5_lm) #extracting the linear regression information
  LM_dbh_ph_0.5_slopes <- c(LM_dbh_ph_0.5_slopes, LM_dbh_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_ph_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LM_dbh_ph_0_5_lm_real_sum <- summary(LM_dbh_ph_0_5_lm_real) #extract the summary 
LM_dbh_ph_0_5_lm_real_slope <- LM_dbh_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_ph_0.5_slopes[i] > LM_dbh_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_ph_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LM_dbh_ph_100_200_lm_sum <- summary(LM_dbh_ph_100_200_lm) #extracting the linear regression information
  LM_dbh_ph_100_200_slopes <- c(LM_dbh_ph_100_200_slopes, LM_dbh_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_ph_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LM_dbh_ph_100_200_lm_real_sum <- summary(LM_dbh_ph_100_200_lm_real) #extract the summary 
LM_dbh_ph_100_200_lm_real_slope <- LM_dbh_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_ph_100_200_slopes[i] > LM_dbh_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_soc_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LM_dbh_soc_0_5_lm_sum <- summary(LM_dbh_soc_0_5_lm) #extracting the linear regression information
  LM_dbh_soc_0.5_slopes <- c(LM_dbh_soc_0.5_slopes, LM_dbh_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_soc_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LM_dbh_soc_0_5_lm_real_sum <- summary(LM_dbh_soc_0_5_lm_real) #extract the summary 
LM_dbh_soc_0_5_lm_real_slope <- LM_dbh_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_soc_0.5_slopes[i] < LM_dbh_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_soc_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LM_dbh_soc_100_200_lm_sum <- summary(LM_dbh_soc_100_200_lm) #extracting the linear regression information
  LM_dbh_soc_100_200_slopes <- c(LM_dbh_soc_100_200_slopes, LM_dbh_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_soc_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LM_dbh_soc_100_200_lm_real_sum <- summary(LM_dbh_soc_100_200_lm_real) #extract the summary 
LM_dbh_soc_100_200_lm_real_slope <- LM_dbh_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_soc_100_200_slopes[i] < LM_dbh_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_vol_10_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LM_dbh_vol_10_0_5_lm_sum <- summary(LM_dbh_vol_10_0_5_lm) #extracting the linear regression information
  LM_dbh_vol_10_0.5_slopes <- c(LM_dbh_vol_10_0.5_slopes, LM_dbh_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_vol_10_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LM_dbh_vol_10_0_5_lm_real_sum <- summary(LM_dbh_vol_10_0_5_lm_real) #extract the summary 
LM_dbh_vol_10_0_5_lm_real_slope <- LM_dbh_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_vol_10_0.5_slopes[i] > LM_dbh_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_vol_10_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LM_dbh_vol_10_100_200_lm_sum <- summary(LM_dbh_vol_10_100_200_lm) #extracting the linear regression information
  LM_dbh_vol_10_100_200_slopes <- c(LM_dbh_vol_10_100_200_slopes, LM_dbh_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_vol_10_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LM_dbh_vol_10_100_200_lm_real_sum <- summary(LM_dbh_vol_10_100_200_lm_real) #extract the summary 
LM_dbh_vol_10_100_200_lm_real_slope <- LM_dbh_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_vol_10_100_200_slopes[i] < LM_dbh_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_vol_1500_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LM_dbh_vol_1500_0_5_lm_sum <- summary(LM_dbh_vol_1500_0_5_lm) #extracting the linear regression information
  LM_dbh_vol_1500_0.5_slopes <- c(LM_dbh_vol_1500_0.5_slopes, LM_dbh_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_vol_1500_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LM_dbh_vol_1500_0_5_lm_real_sum <- summary(LM_dbh_vol_1500_0_5_lm_real) #extract the summary 
LM_dbh_vol_1500_0_5_lm_real_slope <- LM_dbh_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_vol_1500_0.5_slopes[i] > LM_dbh_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_vol_1500_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LM_dbh_vol_1500_100_200_lm_sum <- summary(LM_dbh_vol_1500_100_200_lm) #extracting the linear regression information
  LM_dbh_vol_1500_100_200_slopes <- c(LM_dbh_vol_1500_100_200_slopes, LM_dbh_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_vol_1500_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LM_dbh_vol_1500_100_200_lm_real_sum <- summary(LM_dbh_vol_1500_100_200_lm_real) #extract the summary 
LM_dbh_vol_1500_100_200_lm_real_slope <- LM_dbh_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_vol_1500_100_200_slopes[i] < LM_dbh_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_nitrogen_0_5_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LM_dbh_nitrogen_0_5_lm_sum <- summary(LM_dbh_nitrogen_0_5_lm) #extracting the linear regression information
  LM_dbh_nitrogen_0.5_slopes <- c(LM_dbh_nitrogen_0.5_slopes, LM_dbh_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_nitrogen_0_5_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LM_dbh_nitrogen_0_5_lm_real_sum <- summary(LM_dbh_nitrogen_0_5_lm_real) #extract the summary 
LM_dbh_nitrogen_0_5_lm_real_slope <- LM_dbh_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_nitrogen_0.5_slopes[i] > LM_dbh_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LM_dbh_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LM_fixed_field_data_processed_soils_shuffled <- transform(LM_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LM_dbh_nitrogen_100_200_lm <- lm(LM_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LM_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LM_dbh_nitrogen_100_200_lm_sum <- summary(LM_dbh_nitrogen_100_200_lm) #extracting the linear regression information
  LM_dbh_nitrogen_100_200_slopes <- c(LM_dbh_nitrogen_100_200_slopes, LM_dbh_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LM_dbh_nitrogen_100_200_lm_real <- lm(LM_fixed_field_data_processed_soils$DBH_ag~LM_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LM_dbh_nitrogen_100_200_lm_real_sum <- summary(LM_dbh_nitrogen_100_200_lm_real) #extract the summary 
LM_dbh_nitrogen_100_200_lm_real_slope <- LM_dbh_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LM_dbh_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LM_dbh_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LM_dbh_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LM_dbh_nitrogen_100_200_slopes[i] > LM_dbh_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LM_dbh_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN





### LC ###

## sca

# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_clay_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LC_sca_clay_0_5_LC_sum <- summary(LC_sca_clay_0_5_lm) #extracting the linear regression information
  LC_sca_clay_0.5_slopes <- c(LC_sca_clay_0.5_slopes, LC_sca_clay_0_5_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_clay_0_5_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LC_sca_clay_0_5_LC_real_sum <- summary(LC_sca_clay_0_5_LC_real) #extract the summary 
LC_sca_clay_0_5_LC_real_slope <- LC_sca_clay_0_5_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_clay_0_5_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_clay_0.5_slopes[i] > LC_sca_clay_0_5_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN
  
# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_clay_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LC_sca_clay_100_200_LC_sum <- summary(LC_sca_clay_100_200_lm) #extracting the linear regression information
  LC_sca_clay_100_200_slopes <- c(LC_sca_clay_100_200_slopes, LC_sca_clay_100_200_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_clay_100_200_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LC_sca_clay_100_200_LC_real_sum <- summary(LC_sca_clay_100_200_LC_real) #extract the summary 
LC_sca_clay_100_200_LC_real_slope <- LC_sca_clay_100_200_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_clay_100_200_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_clay_100_200_slopes[i] < LC_sca_clay_100_200_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_silt_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LC_sca_silt_0_5_LC_sum <- summary(LC_sca_silt_0_5_lm) #extracting the linear regression information
  LC_sca_silt_0.5_slopes <- c(LC_sca_silt_0.5_slopes, LC_sca_silt_0_5_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_silt_0_5_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LC_sca_silt_0_5_LC_real_sum <- summary(LC_sca_silt_0_5_LC_real) #extract the summary 
LC_sca_silt_0_5_LC_real_slope <- LC_sca_silt_0_5_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_silt_0_5_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_silt_0.5_slopes[i] > LC_sca_silt_0_5_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_silt_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LC_sca_silt_100_200_LC_sum <- summary(LC_sca_silt_100_200_lm) #extracting the linear regression information
  LC_sca_silt_100_200_slopes <- c(LC_sca_silt_100_200_slopes, LC_sca_silt_100_200_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_silt_100_200_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LC_sca_silt_100_200_LC_real_sum <- summary(LC_sca_silt_100_200_LC_real) #extract the summary 
LC_sca_silt_100_200_LC_real_slope <- LC_sca_silt_100_200_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_silt_100_200_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_silt_100_200_slopes[i] > LC_sca_silt_100_200_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_sand_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LC_sca_sand_0_5_LC_sum <- summary(LC_sca_sand_0_5_lm) #extracting the linear regression information
  LC_sca_sand_0.5_slopes <- c(LC_sca_sand_0.5_slopes, LC_sca_sand_0_5_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_sand_0_5_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LC_sca_sand_0_5_LC_real_sum <- summary(LC_sca_sand_0_5_LC_real) #extract the summary 
LC_sca_sand_0_5_LC_real_slope <- LC_sca_sand_0_5_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_sand_0_5_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_sand_0.5_slopes[i] < LC_sca_sand_0_5_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_sand_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LC_sca_sand_100_200_LC_sum <- summary(LC_sca_sand_100_200_lm) #extracting the linear regression information
  LC_sca_sand_100_200_slopes <- c(LC_sca_sand_100_200_slopes, LC_sca_sand_100_200_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_sand_100_200_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LC_sca_sand_100_200_LC_real_sum <- summary(LC_sca_sand_100_200_LC_real) #extract the summary 
LC_sca_sand_100_200_LC_real_slope <- LC_sca_sand_100_200_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_sand_100_200_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_sand_100_200_slopes[i] < LC_sca_sand_100_200_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_ph_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LC_sca_ph_0_5_LC_sum <- summary(LC_sca_ph_0_5_lm) #extracting the linear regression information
  LC_sca_ph_0.5_slopes <- c(LC_sca_ph_0.5_slopes, LC_sca_ph_0_5_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_ph_0_5_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LC_sca_ph_0_5_LC_real_sum <- summary(LC_sca_ph_0_5_LC_real) #extract the summary 
LC_sca_ph_0_5_LC_real_slope <- LC_sca_ph_0_5_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_ph_0_5_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_ph_0.5_slopes[i] > LC_sca_ph_0_5_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_ph_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LC_sca_ph_100_200_LC_sum <- summary(LC_sca_ph_100_200_lm) #extracting the linear regression information
  LC_sca_ph_100_200_slopes <- c(LC_sca_ph_100_200_slopes, LC_sca_ph_100_200_LC_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_ph_100_200_LC_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LC_sca_ph_100_200_LC_real_sum <- summary(LC_sca_ph_100_200_LC_real) #extract the summary 
LC_sca_ph_100_200_LC_real_slope <- LC_sca_ph_100_200_LC_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_ph_100_200_LC_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_ph_100_200_slopes[i] < LC_sca_ph_100_200_LC_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_soc_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LC_sca_soc_0_5_lm_sum <- summary(LC_sca_soc_0_5_lm) #extracting the linear regression information
  LC_sca_soc_0.5_slopes <- c(LC_sca_soc_0.5_slopes, LC_sca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_soc_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LC_sca_soc_0_5_lm_real_sum <- summary(LC_sca_soc_0_5_lm_real) #extract the summary 
LC_sca_soc_0_5_lm_real_slope <- LC_sca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_soc_0.5_slopes[i] > LC_sca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_soc_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LC_sca_soc_100_200_lm_sum <- summary(LC_sca_soc_100_200_lm) #extracting the linear regression information
  LC_sca_soc_100_200_slopes <- c(LC_sca_soc_100_200_slopes, LC_sca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_soc_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LC_sca_soc_100_200_lm_real_sum <- summary(LC_sca_soc_100_200_lm_real) #extract the summary 
LC_sca_soc_100_200_lm_real_slope <- LC_sca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_soc_100_200_slopes[i] > LC_sca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_vol_10_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LC_sca_vol_10_0_5_lm_sum <- summary(LC_sca_vol_10_0_5_lm) #extracting the linear regression information
  LC_sca_vol_10_0.5_slopes <- c(LC_sca_vol_10_0.5_slopes, LC_sca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_vol_10_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LC_sca_vol_10_0_5_lm_real_sum <- summary(LC_sca_vol_10_0_5_lm_real) #extract the summary 
LC_sca_vol_10_0_5_lm_real_slope <- LC_sca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_vol_10_0.5_slopes[i] < LC_sca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_vol_10_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LC_sca_vol_10_100_200_lm_sum <- summary(LC_sca_vol_10_100_200_lm) #extracting the linear regression information
  LC_sca_vol_10_100_200_slopes <- c(LC_sca_vol_10_100_200_slopes, LC_sca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_vol_10_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LC_sca_vol_10_100_200_lm_real_sum <- summary(LC_sca_vol_10_100_200_lm_real) #extract the summary 
LC_sca_vol_10_100_200_lm_real_slope <- LC_sca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_vol_10_100_200_slopes[i] < LC_sca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_vol_1500_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LC_sca_vol_1500_0_5_lm_sum <- summary(LC_sca_vol_1500_0_5_lm) #extracting the linear regression information
  LC_sca_vol_1500_0.5_slopes <- c(LC_sca_vol_1500_0.5_slopes, LC_sca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_vol_1500_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LC_sca_vol_1500_0_5_lm_real_sum <- summary(LC_sca_vol_1500_0_5_lm_real) #extract the summary 
LC_sca_vol_1500_0_5_lm_real_slope <- LC_sca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_vol_1500_0.5_slopes[i] > LC_sca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_vol_1500_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LC_sca_vol_1500_100_200_lm_sum <- summary(LC_sca_vol_1500_100_200_lm) #extracting the linear regression information
  LC_sca_vol_1500_100_200_slopes <- c(LC_sca_vol_1500_100_200_slopes, LC_sca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_vol_1500_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LC_sca_vol_1500_100_200_lm_real_sum <- summary(LC_sca_vol_1500_100_200_lm_real) #extract the summary 
LC_sca_vol_1500_100_200_lm_real_slope <- LC_sca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_vol_1500_100_200_slopes[i] < LC_sca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_nitrogen_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LC_sca_nitrogen_0_5_lm_sum <- summary(LC_sca_nitrogen_0_5_lm) #extracting the linear regression information
  LC_sca_nitrogen_0.5_slopes <- c(LC_sca_nitrogen_0.5_slopes, LC_sca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_nitrogen_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LC_sca_nitrogen_0_5_lm_real_sum <- summary(LC_sca_nitrogen_0_5_lm_real) #extract the summary 
LC_sca_nitrogen_0_5_lm_real_slope <- LC_sca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_nitrogen_0.5_slopes[i] > LC_sca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_sca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  LC_sca_nitrogen_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LC_sca_nitrogen_100_200_lm_sum <- summary(LC_sca_nitrogen_100_200_lm) #extracting the linear regression information
  LC_sca_nitrogen_100_200_slopes <- c(LC_sca_nitrogen_100_200_slopes, LC_sca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
LC_sca_nitrogen_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_short~LC_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LC_sca_nitrogen_100_200_lm_real_sum <- summary(LC_sca_nitrogen_100_200_lm_real) #extract the summary 
LC_sca_nitrogen_100_200_lm_real_slope <- LC_sca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_sca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_sca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_sca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_sca_nitrogen_100_200_slopes[i] < LC_sca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_sca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# LCA



# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_clay_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LC_lca_clay_0_5_lm_sum <- summary(LC_lca_clay_0_5_lm) #extracting the linear regression information
  LC_lca_clay_0.5_slopes <- c(LC_lca_clay_0.5_slopes, LC_lca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_clay_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LC_lca_clay_0_5_lm_real_sum <- summary(LC_lca_clay_0_5_lm_real) #extract the summary 
LC_lca_clay_0_5_lm_real_slope <- LC_lca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_clay_0.5_slopes[i] > LC_lca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_clay_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LC_lca_clay_100_200_lm_sum <- summary(LC_lca_clay_100_200_lm) #extracting the linear regression information
  LC_lca_clay_100_200_slopes <- c(LC_lca_clay_100_200_slopes, LC_lca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_clay_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LC_lca_clay_100_200_lm_real_sum <- summary(LC_lca_clay_100_200_lm_real) #extract the summary 
LC_lca_clay_100_200_lm_real_slope <- LC_lca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_clay_100_200_slopes[i] < LC_lca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_silt_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LC_lca_silt_0_5_lm_sum <- summary(LC_lca_silt_0_5_lm) #extracting the linear regression information
  LC_lca_silt_0.5_slopes <- c(LC_lca_silt_0.5_slopes, LC_lca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_silt_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LC_lca_silt_0_5_lm_real_sum <- summary(LC_lca_silt_0_5_lm_real) #extract the summary 
LC_lca_silt_0_5_lm_real_slope <- LC_lca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_silt_0.5_slopes[i] > LC_lca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_silt_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LC_lca_silt_100_200_lm_sum <- summary(LC_lca_silt_100_200_lm) #extracting the linear regression information
  LC_lca_silt_100_200_slopes <- c(LC_lca_silt_100_200_slopes, LC_lca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_silt_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LC_lca_silt_100_200_lm_real_sum <- summary(LC_lca_silt_100_200_lm_real) #extract the summary 
LC_lca_silt_100_200_lm_real_slope <- LC_lca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_silt_100_200_slopes[i] > LC_lca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_sand_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LC_lca_sand_0_5_lm_sum <- summary(LC_lca_sand_0_5_lm) #extracting the linear regression information
  LC_lca_sand_0.5_slopes <- c(LC_lca_sand_0.5_slopes, LC_lca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_sand_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LC_lca_sand_0_5_lm_real_sum <- summary(LC_lca_sand_0_5_lm_real) #extract the summary 
LC_lca_sand_0_5_lm_real_slope <- LC_lca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_sand_0.5_slopes[i] < LC_lca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_sand_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LC_lca_sand_100_200_lm_sum <- summary(LC_lca_sand_100_200_lm) #extracting the linear regression information
  LC_lca_sand_100_200_slopes <- c(LC_lca_sand_100_200_slopes, LC_lca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_sand_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LC_lca_sand_100_200_lm_real_sum <- summary(LC_lca_sand_100_200_lm_real) #extract the summary 
LC_lca_sand_100_200_lm_real_slope <- LC_lca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_sand_100_200_slopes[i] < LC_lca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_ph_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LC_lca_ph_0_5_lm_sum <- summary(LC_lca_ph_0_5_lm) #extracting the linear regression information
  LC_lca_ph_0.5_slopes <- c(LC_lca_ph_0.5_slopes, LC_lca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_ph_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LC_lca_ph_0_5_lm_real_sum <- summary(LC_lca_ph_0_5_lm_real) #extract the summary 
LC_lca_ph_0_5_lm_real_slope <- LC_lca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_ph_0.5_slopes[i] > LC_lca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_ph_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LC_lca_ph_100_200_lm_sum <- summary(LC_lca_ph_100_200_lm) #extracting the linear regression information
  LC_lca_ph_100_200_slopes <- c(LC_lca_ph_100_200_slopes, LC_lca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_ph_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LC_lca_ph_100_200_lm_real_sum <- summary(LC_lca_ph_100_200_lm_real) #extract the summary 
LC_lca_ph_100_200_lm_real_slope <- LC_lca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_ph_100_200_slopes[i] < LC_lca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_soc_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LC_lca_soc_0_5_lm_sum <- summary(LC_lca_soc_0_5_lm) #extracting the linear regression information
  LC_lca_soc_0.5_slopes <- c(LC_lca_soc_0.5_slopes, LC_lca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_soc_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LC_lca_soc_0_5_lm_real_sum <- summary(LC_lca_soc_0_5_lm_real) #extract the summary 
LC_lca_soc_0_5_lm_real_slope <- LC_lca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_soc_0.5_slopes[i] > LC_lca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_soc_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LC_lca_soc_100_200_lm_sum <- summary(LC_lca_soc_100_200_lm) #extracting the linear regression information
  LC_lca_soc_100_200_slopes <- c(LC_lca_soc_100_200_slopes, LC_lca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_soc_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LC_lca_soc_100_200_lm_real_sum <- summary(LC_lca_soc_100_200_lm_real) #extract the summary 
LC_lca_soc_100_200_lm_real_slope <- LC_lca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_soc_100_200_slopes[i] > LC_lca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_vol_10_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LC_lca_vol_10_0_5_lm_sum <- summary(LC_lca_vol_10_0_5_lm) #extracting the linear regression information
  LC_lca_vol_10_0.5_slopes <- c(LC_lca_vol_10_0.5_slopes, LC_lca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_vol_10_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LC_lca_vol_10_0_5_lm_real_sum <- summary(LC_lca_vol_10_0_5_lm_real) #extract the summary 
LC_lca_vol_10_0_5_lm_real_slope <- LC_lca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_vol_10_0.5_slopes[i] < LC_lca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_vol_10_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LC_lca_vol_10_100_200_lm_sum <- summary(LC_lca_vol_10_100_200_lm) #extracting the linear regression information
  LC_lca_vol_10_100_200_slopes <- c(LC_lca_vol_10_100_200_slopes, LC_lca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_vol_10_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LC_lca_vol_10_100_200_lm_real_sum <- summary(LC_lca_vol_10_100_200_lm_real) #extract the summary 
LC_lca_vol_10_100_200_lm_real_slope <- LC_lca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_vol_10_100_200_slopes[i] < LC_lca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_vol_1500_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LC_lca_vol_1500_0_5_lm_sum <- summary(LC_lca_vol_1500_0_5_lm) #extracting the linear regression information
  LC_lca_vol_1500_0.5_slopes <- c(LC_lca_vol_1500_0.5_slopes, LC_lca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_vol_1500_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LC_lca_vol_1500_0_5_lm_real_sum <- summary(LC_lca_vol_1500_0_5_lm_real) #extract the summary 
LC_lca_vol_1500_0_5_lm_real_slope <- LC_lca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_vol_1500_0.5_slopes[i] > LC_lca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_vol_1500_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LC_lca_vol_1500_100_200_lm_sum <- summary(LC_lca_vol_1500_100_200_lm) #extracting the linear regression information
  LC_lca_vol_1500_100_200_slopes <- c(LC_lca_vol_1500_100_200_slopes, LC_lca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_vol_1500_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LC_lca_vol_1500_100_200_lm_real_sum <- summary(LC_lca_vol_1500_100_200_lm_real) #extract the summary 
LC_lca_vol_1500_100_200_lm_real_slope <- LC_lca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_vol_1500_100_200_slopes[i] < LC_lca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_nitrogen_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LC_lca_nitrogen_0_5_lm_sum <- summary(LC_lca_nitrogen_0_5_lm) #extracting the linear regression information
  LC_lca_nitrogen_0.5_slopes <- c(LC_lca_nitrogen_0.5_slopes, LC_lca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_nitrogen_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LC_lca_nitrogen_0_5_lm_real_sum <- summary(LC_lca_nitrogen_0_5_lm_real) #extract the summary 
LC_lca_nitrogen_0_5_lm_real_slope <- LC_lca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_nitrogen_0.5_slopes[i] > LC_lca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_lca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  LC_lca_nitrogen_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LC_lca_nitrogen_100_200_lm_sum <- summary(LC_lca_nitrogen_100_200_lm) #extracting the linear regression information
  LC_lca_nitrogen_100_200_slopes <- c(LC_lca_nitrogen_100_200_slopes, LC_lca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
LC_lca_nitrogen_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_long~LC_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LC_lca_nitrogen_100_200_lm_real_sum <- summary(LC_lca_nitrogen_100_200_lm_real) #extract the summary 
LC_lca_nitrogen_100_200_lm_real_slope <- LC_lca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_lca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_lca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_lca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_lca_nitrogen_100_200_slopes[i] > LC_lca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_lca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CA


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_clay_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LC_ca_clay_0_5_lm_sum <- summary(LC_ca_clay_0_5_lm) #extracting the linear regression information
  LC_ca_clay_0.5_slopes <- c(LC_ca_clay_0.5_slopes, LC_ca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_clay_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LC_ca_clay_0_5_lm_real_sum <- summary(LC_ca_clay_0_5_lm_real) #extract the summary 
LC_ca_clay_0_5_lm_real_slope <- LC_ca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_clay_0.5_slopes[i] > LC_ca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_clay_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LC_ca_clay_100_200_lm_sum <- summary(LC_ca_clay_100_200_lm) #extracting the linear regression information
  LC_ca_clay_100_200_slopes <- c(LC_ca_clay_100_200_slopes, LC_ca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_clay_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LC_ca_clay_100_200_lm_real_sum <- summary(LC_ca_clay_100_200_lm_real) #extract the summary 
LC_ca_clay_100_200_lm_real_slope <- LC_ca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_clay_100_200_slopes[i] < LC_ca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_silt_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LC_ca_silt_0_5_lm_sum <- summary(LC_ca_silt_0_5_lm) #extracting the linear regression information
  LC_ca_silt_0.5_slopes <- c(LC_ca_silt_0.5_slopes, LC_ca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_silt_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LC_ca_silt_0_5_lm_real_sum <- summary(LC_ca_silt_0_5_lm_real) #extract the summary 
LC_ca_silt_0_5_lm_real_slope <- LC_ca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_silt_0.5_slopes[i] > LC_ca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_silt_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LC_ca_silt_100_200_lm_sum <- summary(LC_ca_silt_100_200_lm) #extracting the linear regression information
  LC_ca_silt_100_200_slopes <- c(LC_ca_silt_100_200_slopes, LC_ca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_silt_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LC_ca_silt_100_200_lm_real_sum <- summary(LC_ca_silt_100_200_lm_real) #extract the summary 
LC_ca_silt_100_200_lm_real_slope <- LC_ca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_silt_100_200_slopes[i] > LC_ca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_sand_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LC_ca_sand_0_5_lm_sum <- summary(LC_ca_sand_0_5_lm) #extracting the linear regression information
  LC_ca_sand_0.5_slopes <- c(LC_ca_sand_0.5_slopes, LC_ca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_sand_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LC_ca_sand_0_5_lm_real_sum <- summary(LC_ca_sand_0_5_lm_real) #extract the summary 
LC_ca_sand_0_5_lm_real_slope <- LC_ca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_sand_0.5_slopes[i] < LC_ca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_sand_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LC_ca_sand_100_200_lm_sum <- summary(LC_ca_sand_100_200_lm) #extracting the linear regression information
  LC_ca_sand_100_200_slopes <- c(LC_ca_sand_100_200_slopes, LC_ca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_sand_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LC_ca_sand_100_200_lm_real_sum <- summary(LC_ca_sand_100_200_lm_real) #extract the summary 
LC_ca_sand_100_200_lm_real_slope <- LC_ca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_sand_100_200_slopes[i] < LC_ca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_ph_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LC_ca_ph_0_5_lm_sum <- summary(LC_ca_ph_0_5_lm) #extracting the linear regression information
  LC_ca_ph_0.5_slopes <- c(LC_ca_ph_0.5_slopes, LC_ca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_ph_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LC_ca_ph_0_5_lm_real_sum <- summary(LC_ca_ph_0_5_lm_real) #extract the summary 
LC_ca_ph_0_5_lm_real_slope <- LC_ca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_ph_0.5_slopes[i] > LC_ca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_ph_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LC_ca_ph_100_200_lm_sum <- summary(LC_ca_ph_100_200_lm) #extracting the linear regression information
  LC_ca_ph_100_200_slopes <- c(LC_ca_ph_100_200_slopes, LC_ca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_ph_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LC_ca_ph_100_200_lm_real_sum <- summary(LC_ca_ph_100_200_lm_real) #extract the summary 
LC_ca_ph_100_200_lm_real_slope <- LC_ca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_ph_100_200_slopes[i] < LC_ca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_soc_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LC_ca_soc_0_5_lm_sum <- summary(LC_ca_soc_0_5_lm) #extracting the linear regression information
  LC_ca_soc_0.5_slopes <- c(LC_ca_soc_0.5_slopes, LC_ca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_soc_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LC_ca_soc_0_5_lm_real_sum <- summary(LC_ca_soc_0_5_lm_real) #extract the summary 
LC_ca_soc_0_5_lm_real_slope <- LC_ca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_soc_0.5_slopes[i] > LC_ca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_soc_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LC_ca_soc_100_200_lm_sum <- summary(LC_ca_soc_100_200_lm) #extracting the linear regression information
  LC_ca_soc_100_200_slopes <- c(LC_ca_soc_100_200_slopes, LC_ca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_soc_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LC_ca_soc_100_200_lm_real_sum <- summary(LC_ca_soc_100_200_lm_real) #extract the summary 
LC_ca_soc_100_200_lm_real_slope <- LC_ca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_soc_100_200_slopes[i] > LC_ca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_vol_10_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LC_ca_vol_10_0_5_lm_sum <- summary(LC_ca_vol_10_0_5_lm) #extracting the linear regression information
  LC_ca_vol_10_0.5_slopes <- c(LC_ca_vol_10_0.5_slopes, LC_ca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_vol_10_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LC_ca_vol_10_0_5_lm_real_sum <- summary(LC_ca_vol_10_0_5_lm_real) #extract the summary 
LC_ca_vol_10_0_5_lm_real_slope <- LC_ca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_vol_10_0.5_slopes[i] < LC_ca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_vol_10_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LC_ca_vol_10_100_200_lm_sum <- summary(LC_ca_vol_10_100_200_lm) #extracting the linear regression information
  LC_ca_vol_10_100_200_slopes <- c(LC_ca_vol_10_100_200_slopes, LC_ca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_vol_10_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LC_ca_vol_10_100_200_lm_real_sum <- summary(LC_ca_vol_10_100_200_lm_real) #extract the summary 
LC_ca_vol_10_100_200_lm_real_slope <- LC_ca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_vol_10_100_200_slopes[i] < LC_ca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_vol_1500_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LC_ca_vol_1500_0_5_lm_sum <- summary(LC_ca_vol_1500_0_5_lm) #extracting the linear regression information
  LC_ca_vol_1500_0.5_slopes <- c(LC_ca_vol_1500_0.5_slopes, LC_ca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_vol_1500_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LC_ca_vol_1500_0_5_lm_real_sum <- summary(LC_ca_vol_1500_0_5_lm_real) #extract the summary 
LC_ca_vol_1500_0_5_lm_real_slope <- LC_ca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_vol_1500_0.5_slopes[i] > LC_ca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_vol_1500_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LC_ca_vol_1500_100_200_lm_sum <- summary(LC_ca_vol_1500_100_200_lm) #extracting the linear regression information
  LC_ca_vol_1500_100_200_slopes <- c(LC_ca_vol_1500_100_200_slopes, LC_ca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_vol_1500_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LC_ca_vol_1500_100_200_lm_real_sum <- summary(LC_ca_vol_1500_100_200_lm_real) #extract the summary 
LC_ca_vol_1500_100_200_lm_real_slope <- LC_ca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_vol_1500_100_200_slopes[i] < LC_ca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_nitrogen_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LC_ca_nitrogen_0_5_lm_sum <- summary(LC_ca_nitrogen_0_5_lm) #extracting the linear regression information
  LC_ca_nitrogen_0.5_slopes <- c(LC_ca_nitrogen_0.5_slopes, LC_ca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_nitrogen_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LC_ca_nitrogen_0_5_lm_real_sum <- summary(LC_ca_nitrogen_0_5_lm_real) #extract the summary 
LC_ca_nitrogen_0_5_lm_real_slope <- LC_ca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_nitrogen_0.5_slopes[i] > LC_ca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_ca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  LC_ca_nitrogen_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LC_ca_nitrogen_100_200_lm_sum <- summary(LC_ca_nitrogen_100_200_lm) #extracting the linear regression information
  LC_ca_nitrogen_100_200_slopes <- c(LC_ca_nitrogen_100_200_slopes, LC_ca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
LC_ca_nitrogen_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Canopy_area~LC_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LC_ca_nitrogen_100_200_lm_real_sum <- summary(LC_ca_nitrogen_100_200_lm_real) #extract the summary 
LC_ca_nitrogen_100_200_lm_real_slope <- LC_ca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_ca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_ca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_ca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_ca_nitrogen_100_200_slopes[i] < LC_ca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_ca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CS



# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_clay_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LC_cs_clay_0_5_lm_sum <- summary(LC_cs_clay_0_5_lm) #extracting the linear regression information
  LC_cs_clay_0.5_slopes <- c(LC_cs_clay_0.5_slopes, LC_cs_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_clay_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LC_cs_clay_0_5_lm_real_sum <- summary(LC_cs_clay_0_5_lm_real) #extract the summary 
LC_cs_clay_0_5_lm_real_slope <- LC_cs_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_clay_0.5_slopes[i] > LC_cs_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_clay_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LC_cs_clay_100_200_lm_sum <- summary(LC_cs_clay_100_200_lm) #extracting the linear regression information
  LC_cs_clay_100_200_slopes <- c(LC_cs_clay_100_200_slopes, LC_cs_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_clay_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LC_cs_clay_100_200_lm_real_sum <- summary(LC_cs_clay_100_200_lm_real) #extract the summary 
LC_cs_clay_100_200_lm_real_slope <- LC_cs_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_clay_100_200_slopes[i] < LC_cs_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_silt_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LC_cs_silt_0_5_lm_sum <- summary(LC_cs_silt_0_5_lm) #extracting the linear regression information
  LC_cs_silt_0.5_slopes <- c(LC_cs_silt_0.5_slopes, LC_cs_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_silt_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LC_cs_silt_0_5_lm_real_sum <- summary(LC_cs_silt_0_5_lm_real) #extract the summary 
LC_cs_silt_0_5_lm_real_slope <- LC_cs_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_silt_0.5_slopes[i] > LC_cs_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_silt_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LC_cs_silt_100_200_lm_sum <- summary(LC_cs_silt_100_200_lm) #extracting the linear regression information
  LC_cs_silt_100_200_slopes <- c(LC_cs_silt_100_200_slopes, LC_cs_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_silt_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LC_cs_silt_100_200_lm_real_sum <- summary(LC_cs_silt_100_200_lm_real) #extract the summary 
LC_cs_silt_100_200_lm_real_slope <- LC_cs_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_silt_100_200_slopes[i] > LC_cs_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_sand_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LC_cs_sand_0_5_lm_sum <- summary(LC_cs_sand_0_5_lm) #extracting the linear regression information
  LC_cs_sand_0.5_slopes <- c(LC_cs_sand_0.5_slopes, LC_cs_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_sand_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LC_cs_sand_0_5_lm_real_sum <- summary(LC_cs_sand_0_5_lm_real) #extract the summary 
LC_cs_sand_0_5_lm_real_slope <- LC_cs_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_sand_0.5_slopes[i] < LC_cs_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_sand_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LC_cs_sand_100_200_lm_sum <- summary(LC_cs_sand_100_200_lm) #extracting the linear regression information
  LC_cs_sand_100_200_slopes <- c(LC_cs_sand_100_200_slopes, LC_cs_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_sand_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LC_cs_sand_100_200_lm_real_sum <- summary(LC_cs_sand_100_200_lm_real) #extract the summary 
LC_cs_sand_100_200_lm_real_slope <- LC_cs_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_sand_100_200_slopes[i] < LC_cs_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_ph_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LC_cs_ph_0_5_lm_sum <- summary(LC_cs_ph_0_5_lm) #extracting the linear regression information
  LC_cs_ph_0.5_slopes <- c(LC_cs_ph_0.5_slopes, LC_cs_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_ph_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LC_cs_ph_0_5_lm_real_sum <- summary(LC_cs_ph_0_5_lm_real) #extract the summary 
LC_cs_ph_0_5_lm_real_slope <- LC_cs_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_ph_0.5_slopes[i] > LC_cs_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_ph_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LC_cs_ph_100_200_lm_sum <- summary(LC_cs_ph_100_200_lm) #extracting the linear regression information
  LC_cs_ph_100_200_slopes <- c(LC_cs_ph_100_200_slopes, LC_cs_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_ph_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LC_cs_ph_100_200_lm_real_sum <- summary(LC_cs_ph_100_200_lm_real) #extract the summary 
LC_cs_ph_100_200_lm_real_slope <- LC_cs_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_ph_100_200_slopes[i] < LC_cs_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_soc_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LC_cs_soc_0_5_lm_sum <- summary(LC_cs_soc_0_5_lm) #extracting the linear regression information
  LC_cs_soc_0.5_slopes <- c(LC_cs_soc_0.5_slopes, LC_cs_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_soc_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LC_cs_soc_0_5_lm_real_sum <- summary(LC_cs_soc_0_5_lm_real) #extract the summary 
LC_cs_soc_0_5_lm_real_slope <- LC_cs_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_soc_0.5_slopes[i] > LC_cs_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_soc_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LC_cs_soc_100_200_lm_sum <- summary(LC_cs_soc_100_200_lm) #extracting the linear regression information
  LC_cs_soc_100_200_slopes <- c(LC_cs_soc_100_200_slopes, LC_cs_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_soc_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LC_cs_soc_100_200_lm_real_sum <- summary(LC_cs_soc_100_200_lm_real) #extract the summary 
LC_cs_soc_100_200_lm_real_slope <- LC_cs_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_soc_100_200_slopes[i] > LC_cs_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_vol_10_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LC_cs_vol_10_0_5_lm_sum <- summary(LC_cs_vol_10_0_5_lm) #extracting the linear regression information
  LC_cs_vol_10_0.5_slopes <- c(LC_cs_vol_10_0.5_slopes, LC_cs_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_vol_10_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LC_cs_vol_10_0_5_lm_real_sum <- summary(LC_cs_vol_10_0_5_lm_real) #extract the summary 
LC_cs_vol_10_0_5_lm_real_slope <- LC_cs_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_vol_10_0.5_slopes[i] < LC_cs_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_vol_10_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LC_cs_vol_10_100_200_lm_sum <- summary(LC_cs_vol_10_100_200_lm) #extracting the linear regression information
  LC_cs_vol_10_100_200_slopes <- c(LC_cs_vol_10_100_200_slopes, LC_cs_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_vol_10_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LC_cs_vol_10_100_200_lm_real_sum <- summary(LC_cs_vol_10_100_200_lm_real) #extract the summary 
LC_cs_vol_10_100_200_lm_real_slope <- LC_cs_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_vol_10_100_200_slopes[i] < LC_cs_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_vol_1500_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LC_cs_vol_1500_0_5_lm_sum <- summary(LC_cs_vol_1500_0_5_lm) #extracting the linear regression information
  LC_cs_vol_1500_0.5_slopes <- c(LC_cs_vol_1500_0.5_slopes, LC_cs_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_vol_1500_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LC_cs_vol_1500_0_5_lm_real_sum <- summary(LC_cs_vol_1500_0_5_lm_real) #extract the summary 
LC_cs_vol_1500_0_5_lm_real_slope <- LC_cs_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_vol_1500_0.5_slopes[i] > LC_cs_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_vol_1500_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LC_cs_vol_1500_100_200_lm_sum <- summary(LC_cs_vol_1500_100_200_lm) #extracting the linear regression information
  LC_cs_vol_1500_100_200_slopes <- c(LC_cs_vol_1500_100_200_slopes, LC_cs_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_vol_1500_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LC_cs_vol_1500_100_200_lm_real_sum <- summary(LC_cs_vol_1500_100_200_lm_real) #extract the summary 
LC_cs_vol_1500_100_200_lm_real_slope <- LC_cs_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_vol_1500_100_200_slopes[i] < LC_cs_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_nitrogen_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LC_cs_nitrogen_0_5_lm_sum <- summary(LC_cs_nitrogen_0_5_lm) #extracting the linear regression information
  LC_cs_nitrogen_0.5_slopes <- c(LC_cs_nitrogen_0.5_slopes, LC_cs_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_nitrogen_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LC_cs_nitrogen_0_5_lm_real_sum <- summary(LC_cs_nitrogen_0_5_lm_real) #extract the summary 
LC_cs_nitrogen_0_5_lm_real_slope <- LC_cs_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_nitrogen_0.5_slopes[i] > LC_cs_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_cs_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  LC_cs_nitrogen_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LC_cs_nitrogen_100_200_lm_sum <- summary(LC_cs_nitrogen_100_200_lm) #extracting the linear regression information
  LC_cs_nitrogen_100_200_slopes <- c(LC_cs_nitrogen_100_200_slopes, LC_cs_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
LC_cs_nitrogen_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$Crown_spread~LC_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LC_cs_nitrogen_100_200_lm_real_sum <- summary(LC_cs_nitrogen_100_200_lm_real) #extract the summary 
LC_cs_nitrogen_100_200_lm_real_slope <- LC_cs_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_cs_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_cs_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_cs_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_cs_nitrogen_100_200_slopes[i] > LC_cs_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_cs_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



#DBH_ag


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_clay_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  LC_dbh_clay_0_5_lm_sum <- summary(LC_dbh_clay_0_5_lm) #extracting the linear regression information
  LC_dbh_clay_0.5_slopes <- c(LC_dbh_clay_0.5_slopes, LC_dbh_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_clay_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
LC_dbh_clay_0_5_lm_real_sum <- summary(LC_dbh_clay_0_5_lm_real) #extract the summary 
LC_dbh_clay_0_5_lm_real_slope <- LC_dbh_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_clay_0.5_slopes[i] > LC_dbh_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_clay_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  LC_dbh_clay_100_200_lm_sum <- summary(LC_dbh_clay_100_200_lm) #extracting the linear regression information
  LC_dbh_clay_100_200_slopes <- c(LC_dbh_clay_100_200_slopes, LC_dbh_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_clay_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
LC_dbh_clay_100_200_lm_real_sum <- summary(LC_dbh_clay_100_200_lm_real) #extract the summary 
LC_dbh_clay_100_200_lm_real_slope <- LC_dbh_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_clay_100_200_slopes[i] > LC_dbh_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_silt_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.0.5)
  LC_dbh_silt_0_5_lm_sum <- summary(LC_dbh_silt_0_5_lm) #extracting the linear regression information
  LC_dbh_silt_0.5_slopes <- c(LC_dbh_silt_0.5_slopes, LC_dbh_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_silt_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
LC_dbh_silt_0_5_lm_real_sum <- summary(LC_dbh_silt_0_5_lm_real) #extract the summary 
LC_dbh_silt_0_5_lm_real_slope <- LC_dbh_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_silt_0.5_slopes[i] > LC_dbh_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_silt_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$silt.100.200)
  LC_dbh_silt_100_200_lm_sum <- summary(LC_dbh_silt_100_200_lm) #extracting the linear regression information
  LC_dbh_silt_100_200_slopes <- c(LC_dbh_silt_100_200_slopes, LC_dbh_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_silt_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
LC_dbh_silt_100_200_lm_real_sum <- summary(LC_dbh_silt_100_200_lm_real) #extract the summary 
LC_dbh_silt_100_200_lm_real_slope <- LC_dbh_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_silt_100_200_slopes[i] > LC_dbh_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_sand_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.0.5)
  LC_dbh_sand_0_5_lm_sum <- summary(LC_dbh_sand_0_5_lm) #extracting the linear regression information
  LC_dbh_sand_0.5_slopes <- c(LC_dbh_sand_0.5_slopes, LC_dbh_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_sand_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
LC_dbh_sand_0_5_lm_real_sum <- summary(LC_dbh_sand_0_5_lm_real) #extract the summary 
LC_dbh_sand_0_5_lm_real_slope <- LC_dbh_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_sand_0.5_slopes[i] < LC_dbh_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_sand_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$sand.100.200)
  LC_dbh_sand_100_200_lm_sum <- summary(LC_dbh_sand_100_200_lm) #extracting the linear regression information
  LC_dbh_sand_100_200_slopes <- c(LC_dbh_sand_100_200_slopes, LC_dbh_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_sand_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
LC_dbh_sand_100_200_lm_real_sum <- summary(LC_dbh_sand_100_200_lm_real) #extract the summary 
LC_dbh_sand_100_200_lm_real_slope <- LC_dbh_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_sand_100_200_slopes[i] < LC_dbh_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_ph_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_0.5)
  LC_dbh_ph_0_5_lm_sum <- summary(LC_dbh_ph_0_5_lm) #extracting the linear regression information
  LC_dbh_ph_0.5_slopes <- c(LC_dbh_ph_0.5_slopes, LC_dbh_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_ph_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
LC_dbh_ph_0_5_lm_real_sum <- summary(LC_dbh_ph_0_5_lm_real) #extract the summary 
LC_dbh_ph_0_5_lm_real_slope <- LC_dbh_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_ph_0.5_slopes[i] > LC_dbh_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_ph_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$ph_100.200)
  LC_dbh_ph_100_200_lm_sum <- summary(LC_dbh_ph_100_200_lm) #extracting the linear regression information
  LC_dbh_ph_100_200_slopes <- c(LC_dbh_ph_100_200_slopes, LC_dbh_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_ph_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
LC_dbh_ph_100_200_lm_real_sum <- summary(LC_dbh_ph_100_200_lm_real) #extract the summary 
LC_dbh_ph_100_200_lm_real_slope <- LC_dbh_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_ph_100_200_slopes[i] < LC_dbh_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_soc_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  LC_dbh_soc_0_5_lm_sum <- summary(LC_dbh_soc_0_5_lm) #extracting the linear regression information
  LC_dbh_soc_0.5_slopes <- c(LC_dbh_soc_0.5_slopes, LC_dbh_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_soc_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
LC_dbh_soc_0_5_lm_real_sum <- summary(LC_dbh_soc_0_5_lm_real) #extract the summary 
LC_dbh_soc_0_5_lm_real_slope <- LC_dbh_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_soc_0.5_slopes[i] > LC_dbh_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_soc_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  LC_dbh_soc_100_200_lm_sum <- summary(LC_dbh_soc_100_200_lm) #extracting the linear regression information
  LC_dbh_soc_100_200_slopes <- c(LC_dbh_soc_100_200_slopes, LC_dbh_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_soc_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
LC_dbh_soc_100_200_lm_real_sum <- summary(LC_dbh_soc_100_200_lm_real) #extract the summary 
LC_dbh_soc_100_200_lm_real_slope <- LC_dbh_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_soc_100_200_slopes[i] > LC_dbh_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_vol_10_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  LC_dbh_vol_10_0_5_lm_sum <- summary(LC_dbh_vol_10_0_5_lm) #extracting the linear regression information
  LC_dbh_vol_10_0.5_slopes <- c(LC_dbh_vol_10_0.5_slopes, LC_dbh_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_vol_10_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
LC_dbh_vol_10_0_5_lm_real_sum <- summary(LC_dbh_vol_10_0_5_lm_real) #extract the summary 
LC_dbh_vol_10_0_5_lm_real_slope <- LC_dbh_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_vol_10_0.5_slopes[i] < LC_dbh_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_vol_10_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  LC_dbh_vol_10_100_200_lm_sum <- summary(LC_dbh_vol_10_100_200_lm) #extracting the linear regression information
  LC_dbh_vol_10_100_200_slopes <- c(LC_dbh_vol_10_100_200_slopes, LC_dbh_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_vol_10_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
LC_dbh_vol_10_100_200_lm_real_sum <- summary(LC_dbh_vol_10_100_200_lm_real) #extract the summary 
LC_dbh_vol_10_100_200_lm_real_slope <- LC_dbh_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_vol_10_100_200_slopes[i] < LC_dbh_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_vol_1500_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  LC_dbh_vol_1500_0_5_lm_sum <- summary(LC_dbh_vol_1500_0_5_lm) #extracting the linear regression information
  LC_dbh_vol_1500_0.5_slopes <- c(LC_dbh_vol_1500_0.5_slopes, LC_dbh_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_vol_1500_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
LC_dbh_vol_1500_0_5_lm_real_sum <- summary(LC_dbh_vol_1500_0_5_lm_real) #extract the summary 
LC_dbh_vol_1500_0_5_lm_real_slope <- LC_dbh_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_vol_1500_0.5_slopes[i] > LC_dbh_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_vol_1500_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  LC_dbh_vol_1500_100_200_lm_sum <- summary(LC_dbh_vol_1500_100_200_lm) #extracting the linear regression information
  LC_dbh_vol_1500_100_200_slopes <- c(LC_dbh_vol_1500_100_200_slopes, LC_dbh_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_vol_1500_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
LC_dbh_vol_1500_100_200_lm_real_sum <- summary(LC_dbh_vol_1500_100_200_lm_real) #extract the summary 
LC_dbh_vol_1500_100_200_lm_real_slope <- LC_dbh_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_vol_1500_100_200_slopes[i] > LC_dbh_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_nitrogen_0_5_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  LC_dbh_nitrogen_0_5_lm_sum <- summary(LC_dbh_nitrogen_0_5_lm) #extracting the linear regression information
  LC_dbh_nitrogen_0.5_slopes <- c(LC_dbh_nitrogen_0.5_slopes, LC_dbh_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_nitrogen_0_5_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
LC_dbh_nitrogen_0_5_lm_real_sum <- summary(LC_dbh_nitrogen_0_5_lm_real) #extract the summary 
LC_dbh_nitrogen_0_5_lm_real_slope <- LC_dbh_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_nitrogen_0.5_slopes[i] < LC_dbh_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
LC_dbh_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  LC_fixed_field_data_processed_soils_shuffled <- transform(LC_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  LC_dbh_nitrogen_100_200_lm <- lm(LC_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~LC_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  LC_dbh_nitrogen_100_200_lm_sum <- summary(LC_dbh_nitrogen_100_200_lm) #extracting the linear regression information
  LC_dbh_nitrogen_100_200_slopes <- c(LC_dbh_nitrogen_100_200_slopes, LC_dbh_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
LC_dbh_nitrogen_100_200_lm_real <- lm(LC_fixed_field_data_processed_soils$DBH_ag~LC_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
LC_dbh_nitrogen_100_200_lm_real_sum <- summary(LC_dbh_nitrogen_100_200_lm_real) #extract the summary 
LC_dbh_nitrogen_100_200_lm_real_slope <- LC_dbh_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=LC_dbh_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=LC_dbh_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(LC_dbh_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (LC_dbh_nitrogen_100_200_slopes[i] < LC_dbh_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(LC_dbh_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN




### SD ###

## sca

# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_clay_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  SD_sca_clay_0_5_SD_sum <- summary(SD_sca_clay_0_5_lm) #extracting the linear regression information
  SD_sca_clay_0.5_slopes <- c(SD_sca_clay_0.5_slopes, SD_sca_clay_0_5_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_clay_0_5_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
SD_sca_clay_0_5_SD_real_sum <- summary(SD_sca_clay_0_5_SD_real) #extract the summary 
SD_sca_clay_0_5_SD_real_slope <- SD_sca_clay_0_5_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_clay_0_5_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_clay_0.5_slopes[i] > SD_sca_clay_0_5_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_clay_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  SD_sca_clay_100_200_SD_sum <- summary(SD_sca_clay_100_200_lm) #extracting the linear regression information
  SD_sca_clay_100_200_slopes <- c(SD_sca_clay_100_200_slopes, SD_sca_clay_100_200_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_clay_100_200_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
SD_sca_clay_100_200_SD_real_sum <- summary(SD_sca_clay_100_200_SD_real) #extract the summary 
SD_sca_clay_100_200_SD_real_slope <- SD_sca_clay_100_200_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_clay_100_200_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_clay_100_200_slopes[i] > SD_sca_clay_100_200_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_silt_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.0.5)
  SD_sca_silt_0_5_SD_sum <- summary(SD_sca_silt_0_5_lm) #extracting the linear regression information
  SD_sca_silt_0.5_slopes <- c(SD_sca_silt_0.5_slopes, SD_sca_silt_0_5_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_silt_0_5_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
SD_sca_silt_0_5_SD_real_sum <- summary(SD_sca_silt_0_5_SD_real) #extract the summary 
SD_sca_silt_0_5_SD_real_slope <- SD_sca_silt_0_5_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_silt_0_5_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_silt_0.5_slopes[i] > SD_sca_silt_0_5_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_silt_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.100.200)
  SD_sca_silt_100_200_SD_sum <- summary(SD_sca_silt_100_200_lm) #extracting the linear regression information
  SD_sca_silt_100_200_slopes <- c(SD_sca_silt_100_200_slopes, SD_sca_silt_100_200_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_silt_100_200_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
SD_sca_silt_100_200_SD_real_sum <- summary(SD_sca_silt_100_200_SD_real) #extract the summary 
SD_sca_silt_100_200_SD_real_slope <- SD_sca_silt_100_200_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_silt_100_200_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_silt_100_200_slopes[i] > SD_sca_silt_100_200_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_sand_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.0.5)
  SD_sca_sand_0_5_SD_sum <- summary(SD_sca_sand_0_5_lm) #extracting the linear regression information
  SD_sca_sand_0.5_slopes <- c(SD_sca_sand_0.5_slopes, SD_sca_sand_0_5_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_sand_0_5_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
SD_sca_sand_0_5_SD_real_sum <- summary(SD_sca_sand_0_5_SD_real) #extract the summary 
SD_sca_sand_0_5_SD_real_slope <- SD_sca_sand_0_5_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_sand_0_5_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_sand_0.5_slopes[i] < SD_sca_sand_0_5_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_sand_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.100.200)
  SD_sca_sand_100_200_SD_sum <- summary(SD_sca_sand_100_200_lm) #extracting the linear regression information
  SD_sca_sand_100_200_slopes <- c(SD_sca_sand_100_200_slopes, SD_sca_sand_100_200_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_sand_100_200_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
SD_sca_sand_100_200_SD_real_sum <- summary(SD_sca_sand_100_200_SD_real) #extract the summary 
SD_sca_sand_100_200_SD_real_slope <- SD_sca_sand_100_200_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_sand_100_200_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_sand_100_200_slopes[i] < SD_sca_sand_100_200_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_ph_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_0.5)
  SD_sca_ph_0_5_SD_sum <- summary(SD_sca_ph_0_5_lm) #extracting the linear regression information
  SD_sca_ph_0.5_slopes <- c(SD_sca_ph_0.5_slopes, SD_sca_ph_0_5_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_ph_0_5_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
SD_sca_ph_0_5_SD_real_sum <- summary(SD_sca_ph_0_5_SD_real) #extract the summary 
SD_sca_ph_0_5_SD_real_slope <- SD_sca_ph_0_5_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_ph_0_5_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_ph_0.5_slopes[i] < SD_sca_ph_0_5_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_ph_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_100.200)
  SD_sca_ph_100_200_SD_sum <- summary(SD_sca_ph_100_200_lm) #extracting the linear regression information
  SD_sca_ph_100_200_slopes <- c(SD_sca_ph_100_200_slopes, SD_sca_ph_100_200_SD_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_ph_100_200_SD_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
SD_sca_ph_100_200_SD_real_sum <- summary(SD_sca_ph_100_200_SD_real) #extract the summary 
SD_sca_ph_100_200_SD_real_slope <- SD_sca_ph_100_200_SD_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_ph_100_200_SD_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_ph_100_200_slopes[i] > SD_sca_ph_100_200_SD_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_soc_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  SD_sca_soc_0_5_lm_sum <- summary(SD_sca_soc_0_5_lm) #extracting the linear regression information
  SD_sca_soc_0.5_slopes <- c(SD_sca_soc_0.5_slopes, SD_sca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_soc_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
SD_sca_soc_0_5_lm_real_sum <- summary(SD_sca_soc_0_5_lm_real) #extract the summary 
SD_sca_soc_0_5_lm_real_slope <- SD_sca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_soc_0.5_slopes[i] < SD_sca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_soc_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  SD_sca_soc_100_200_lm_sum <- summary(SD_sca_soc_100_200_lm) #extracting the linear regression information
  SD_sca_soc_100_200_slopes <- c(SD_sca_soc_100_200_slopes, SD_sca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_soc_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
SD_sca_soc_100_200_lm_real_sum <- summary(SD_sca_soc_100_200_lm_real) #extract the summary 
SD_sca_soc_100_200_lm_real_slope <- SD_sca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_soc_100_200_slopes[i] < SD_sca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_vol_10_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  SD_sca_vol_10_0_5_lm_sum <- summary(SD_sca_vol_10_0_5_lm) #extracting the linear regression information
  SD_sca_vol_10_0.5_slopes <- c(SD_sca_vol_10_0.5_slopes, SD_sca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_vol_10_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
SD_sca_vol_10_0_5_lm_real_sum <- summary(SD_sca_vol_10_0_5_lm_real) #extract the summary 
SD_sca_vol_10_0_5_lm_real_slope <- SD_sca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_vol_10_0.5_slopes[i] > SD_sca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_vol_10_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  SD_sca_vol_10_100_200_lm_sum <- summary(SD_sca_vol_10_100_200_lm) #extracting the linear regression information
  SD_sca_vol_10_100_200_slopes <- c(SD_sca_vol_10_100_200_slopes, SD_sca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_vol_10_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
SD_sca_vol_10_100_200_lm_real_sum <- summary(SD_sca_vol_10_100_200_lm_real) #extract the summary 
SD_sca_vol_10_100_200_lm_real_slope <- SD_sca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_vol_10_100_200_slopes[i] < SD_sca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_vol_1500_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  SD_sca_vol_1500_0_5_lm_sum <- summary(SD_sca_vol_1500_0_5_lm) #extracting the linear regression information
  SD_sca_vol_1500_0.5_slopes <- c(SD_sca_vol_1500_0.5_slopes, SD_sca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_vol_1500_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
SD_sca_vol_1500_0_5_lm_real_sum <- summary(SD_sca_vol_1500_0_5_lm_real) #extract the summary 
SD_sca_vol_1500_0_5_lm_real_slope <- SD_sca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_vol_1500_0.5_slopes[i] < SD_sca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_vol_1500_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  SD_sca_vol_1500_100_200_lm_sum <- summary(SD_sca_vol_1500_100_200_lm) #extracting the linear regression information
  SD_sca_vol_1500_100_200_slopes <- c(SD_sca_vol_1500_100_200_slopes, SD_sca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_vol_1500_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
SD_sca_vol_1500_100_200_lm_real_sum <- summary(SD_sca_vol_1500_100_200_lm_real) #extract the summary 
SD_sca_vol_1500_100_200_lm_real_slope <- SD_sca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_vol_1500_100_200_slopes[i] > SD_sca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_nitrogen_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  SD_sca_nitrogen_0_5_lm_sum <- summary(SD_sca_nitrogen_0_5_lm) #extracting the linear regression information
  SD_sca_nitrogen_0.5_slopes <- c(SD_sca_nitrogen_0.5_slopes, SD_sca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_nitrogen_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
SD_sca_nitrogen_0_5_lm_real_sum <- summary(SD_sca_nitrogen_0_5_lm_real) #extract the summary 
SD_sca_nitrogen_0_5_lm_real_slope <- SD_sca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_nitrogen_0.5_slopes[i] > SD_sca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_sca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_short.shuffled = sample(Canopy_short)) #create a data frame with a shuffled 
  SD_sca_nitrogen_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_short.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  SD_sca_nitrogen_100_200_lm_sum <- summary(SD_sca_nitrogen_100_200_lm) #extracting the linear regression information
  SD_sca_nitrogen_100_200_slopes <- c(SD_sca_nitrogen_100_200_slopes, SD_sca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized sca values to the list of stored slopes
}

#extracting the slope of our points
SD_sca_nitrogen_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_short~SD_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression

SD_sca_nitrogen_100_200_lm_real_sum <- summary(SD_sca_nitrogen_100_200_lm_real) #extract the summary 
SD_sca_nitrogen_100_200_lm_real_slope <- SD_sca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_sca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_sca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled sca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_sca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_sca_nitrogen_100_200_slopes[i] < SD_sca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_sca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# LCA



# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_clay_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  SD_lca_clay_0_5_lm_sum <- summary(SD_lca_clay_0_5_lm) #extracting the linear regression information
  SD_lca_clay_0.5_slopes <- c(SD_lca_clay_0.5_slopes, SD_lca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_clay_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
SD_lca_clay_0_5_lm_real_sum <- summary(SD_lca_clay_0_5_lm_real) #extract the summary 
SD_lca_clay_0_5_lm_real_slope <- SD_lca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_clay_0.5_slopes[i] > SD_lca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_clay_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  SD_lca_clay_100_200_lm_sum <- summary(SD_lca_clay_100_200_lm) #extracting the linear regression information
  SD_lca_clay_100_200_slopes <- c(SD_lca_clay_100_200_slopes, SD_lca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_clay_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
SD_lca_clay_100_200_lm_real_sum <- summary(SD_lca_clay_100_200_lm_real) #extract the summary 
SD_lca_clay_100_200_lm_real_slope <- SD_lca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_clay_100_200_slopes[i] > SD_lca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_silt_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.0.5)
  SD_lca_silt_0_5_lm_sum <- summary(SD_lca_silt_0_5_lm) #extracting the linear regression information
  SD_lca_silt_0.5_slopes <- c(SD_lca_silt_0.5_slopes, SD_lca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_silt_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
SD_lca_silt_0_5_lm_real_sum <- summary(SD_lca_silt_0_5_lm_real) #extract the summary 
SD_lca_silt_0_5_lm_real_slope <- SD_lca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_silt_0.5_slopes[i] > SD_lca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_silt_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.100.200)
  SD_lca_silt_100_200_lm_sum <- summary(SD_lca_silt_100_200_lm) #extracting the linear regression information
  SD_lca_silt_100_200_slopes <- c(SD_lca_silt_100_200_slopes, SD_lca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_silt_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
SD_lca_silt_100_200_lm_real_sum <- summary(SD_lca_silt_100_200_lm_real) #extract the summary 
SD_lca_silt_100_200_lm_real_slope <- SD_lca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_silt_100_200_slopes[i] > SD_lca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_sand_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.0.5)
  SD_lca_sand_0_5_lm_sum <- summary(SD_lca_sand_0_5_lm) #extracting the linear regression information
  SD_lca_sand_0.5_slopes <- c(SD_lca_sand_0.5_slopes, SD_lca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_sand_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
SD_lca_sand_0_5_lm_real_sum <- summary(SD_lca_sand_0_5_lm_real) #extract the summary 
SD_lca_sand_0_5_lm_real_slope <- SD_lca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_sand_0.5_slopes[i] < SD_lca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_sand_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.100.200)
  SD_lca_sand_100_200_lm_sum <- summary(SD_lca_sand_100_200_lm) #extracting the linear regression information
  SD_lca_sand_100_200_slopes <- c(SD_lca_sand_100_200_slopes, SD_lca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_sand_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
SD_lca_sand_100_200_lm_real_sum <- summary(SD_lca_sand_100_200_lm_real) #extract the summary 
SD_lca_sand_100_200_lm_real_slope <- SD_lca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_sand_100_200_slopes[i] < SD_lca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_ph_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_0.5)
  SD_lca_ph_0_5_lm_sum <- summary(SD_lca_ph_0_5_lm) #extracting the linear regression information
  SD_lca_ph_0.5_slopes <- c(SD_lca_ph_0.5_slopes, SD_lca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_ph_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
SD_lca_ph_0_5_lm_real_sum <- summary(SD_lca_ph_0_5_lm_real) #extract the summary 
SD_lca_ph_0_5_lm_real_slope <- SD_lca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_ph_0.5_slopes[i] < SD_lca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_ph_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_100.200)
  SD_lca_ph_100_200_lm_sum <- summary(SD_lca_ph_100_200_lm) #extracting the linear regression information
  SD_lca_ph_100_200_slopes <- c(SD_lca_ph_100_200_slopes, SD_lca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_ph_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
SD_lca_ph_100_200_lm_real_sum <- summary(SD_lca_ph_100_200_lm_real) #extract the summary 
SD_lca_ph_100_200_lm_real_slope <- SD_lca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_ph_100_200_slopes[i] < SD_lca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_soc_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  SD_lca_soc_0_5_lm_sum <- summary(SD_lca_soc_0_5_lm) #extracting the linear regression information
  SD_lca_soc_0.5_slopes <- c(SD_lca_soc_0.5_slopes, SD_lca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_soc_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
SD_lca_soc_0_5_lm_real_sum <- summary(SD_lca_soc_0_5_lm_real) #extract the summary 
SD_lca_soc_0_5_lm_real_slope <- SD_lca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_soc_0.5_slopes[i] < SD_lca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_soc_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  SD_lca_soc_100_200_lm_sum <- summary(SD_lca_soc_100_200_lm) #extracting the linear regression information
  SD_lca_soc_100_200_slopes <- c(SD_lca_soc_100_200_slopes, SD_lca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_soc_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
SD_lca_soc_100_200_lm_real_sum <- summary(SD_lca_soc_100_200_lm_real) #extract the summary 
SD_lca_soc_100_200_lm_real_slope <- SD_lca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_soc_100_200_slopes[i] < SD_lca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_vol_10_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  SD_lca_vol_10_0_5_lm_sum <- summary(SD_lca_vol_10_0_5_lm) #extracting the linear regression information
  SD_lca_vol_10_0.5_slopes <- c(SD_lca_vol_10_0.5_slopes, SD_lca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_vol_10_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
SD_lca_vol_10_0_5_lm_real_sum <- summary(SD_lca_vol_10_0_5_lm_real) #extract the summary 
SD_lca_vol_10_0_5_lm_real_slope <- SD_lca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_vol_10_0.5_slopes[i] > SD_lca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_vol_10_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  SD_lca_vol_10_100_200_lm_sum <- summary(SD_lca_vol_10_100_200_lm) #extracting the linear regression information
  SD_lca_vol_10_100_200_slopes <- c(SD_lca_vol_10_100_200_slopes, SD_lca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_vol_10_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
SD_lca_vol_10_100_200_lm_real_sum <- summary(SD_lca_vol_10_100_200_lm_real) #extract the summary 
SD_lca_vol_10_100_200_lm_real_slope <- SD_lca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_vol_10_100_200_slopes[i] < SD_lca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_vol_1500_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  SD_lca_vol_1500_0_5_lm_sum <- summary(SD_lca_vol_1500_0_5_lm) #extracting the linear regression information
  SD_lca_vol_1500_0.5_slopes <- c(SD_lca_vol_1500_0.5_slopes, SD_lca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_vol_1500_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
SD_lca_vol_1500_0_5_lm_real_sum <- summary(SD_lca_vol_1500_0_5_lm_real) #extract the summary 
SD_lca_vol_1500_0_5_lm_real_slope <- SD_lca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_vol_1500_0.5_slopes[i] < SD_lca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_vol_1500_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  SD_lca_vol_1500_100_200_lm_sum <- summary(SD_lca_vol_1500_100_200_lm) #extracting the linear regression information
  SD_lca_vol_1500_100_200_slopes <- c(SD_lca_vol_1500_100_200_slopes, SD_lca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_vol_1500_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
SD_lca_vol_1500_100_200_lm_real_sum <- summary(SD_lca_vol_1500_100_200_lm_real) #extract the summary 
SD_lca_vol_1500_100_200_lm_real_slope <- SD_lca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_vol_1500_100_200_slopes[i] > SD_lca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_nitrogen_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  SD_lca_nitrogen_0_5_lm_sum <- summary(SD_lca_nitrogen_0_5_lm) #extracting the linear regression information
  SD_lca_nitrogen_0.5_slopes <- c(SD_lca_nitrogen_0.5_slopes, SD_lca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_nitrogen_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
SD_lca_nitrogen_0_5_lm_real_sum <- summary(SD_lca_nitrogen_0_5_lm_real) #extract the summary 
SD_lca_nitrogen_0_5_lm_real_slope <- SD_lca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_nitrogen_0.5_slopes[i] > SD_lca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_lca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_long.shuffled = sample(Canopy_long)) #create a data frame with a shuffled 
  SD_lca_nitrogen_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_long.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  SD_lca_nitrogen_100_200_lm_sum <- summary(SD_lca_nitrogen_100_200_lm) #extracting the linear regression information
  SD_lca_nitrogen_100_200_slopes <- c(SD_lca_nitrogen_100_200_slopes, SD_lca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized lca values to the list of stored slopes
}

#extracting the slope of our points
SD_lca_nitrogen_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_long~SD_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
SD_lca_nitrogen_100_200_lm_real_sum <- summary(SD_lca_nitrogen_100_200_lm_real) #extract the summary 
SD_lca_nitrogen_100_200_lm_real_slope <- SD_lca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_lca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_lca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled lca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_lca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_lca_nitrogen_100_200_slopes[i] > SD_lca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_lca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CA


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_clay_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  SD_ca_clay_0_5_lm_sum <- summary(SD_ca_clay_0_5_lm) #extracting the linear regression information
  SD_ca_clay_0.5_slopes <- c(SD_ca_clay_0.5_slopes, SD_ca_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_clay_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
SD_ca_clay_0_5_lm_real_sum <- summary(SD_ca_clay_0_5_lm_real) #extract the summary 
SD_ca_clay_0_5_lm_real_slope <- SD_ca_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_clay_0.5_slopes[i] > SD_ca_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_clay_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  SD_ca_clay_100_200_lm_sum <- summary(SD_ca_clay_100_200_lm) #extracting the linear regression information
  SD_ca_clay_100_200_slopes <- c(SD_ca_clay_100_200_slopes, SD_ca_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_clay_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
SD_ca_clay_100_200_lm_real_sum <- summary(SD_ca_clay_100_200_lm_real) #extract the summary 
SD_ca_clay_100_200_lm_real_slope <- SD_ca_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_clay_100_200_slopes[i] > SD_ca_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_silt_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.0.5)
  SD_ca_silt_0_5_lm_sum <- summary(SD_ca_silt_0_5_lm) #extracting the linear regression information
  SD_ca_silt_0.5_slopes <- c(SD_ca_silt_0.5_slopes, SD_ca_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_silt_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
SD_ca_silt_0_5_lm_real_sum <- summary(SD_ca_silt_0_5_lm_real) #extract the summary 
SD_ca_silt_0_5_lm_real_slope <- SD_ca_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_silt_0.5_slopes[i] > SD_ca_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_silt_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.100.200)
  SD_ca_silt_100_200_lm_sum <- summary(SD_ca_silt_100_200_lm) #extracting the linear regression information
  SD_ca_silt_100_200_slopes <- c(SD_ca_silt_100_200_slopes, SD_ca_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_silt_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
SD_ca_silt_100_200_lm_real_sum <- summary(SD_ca_silt_100_200_lm_real) #extract the summary 
SD_ca_silt_100_200_lm_real_slope <- SD_ca_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_silt_100_200_slopes[i] > SD_ca_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_sand_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.0.5)
  SD_ca_sand_0_5_lm_sum <- summary(SD_ca_sand_0_5_lm) #extracting the linear regression information
  SD_ca_sand_0.5_slopes <- c(SD_ca_sand_0.5_slopes, SD_ca_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_sand_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
SD_ca_sand_0_5_lm_real_sum <- summary(SD_ca_sand_0_5_lm_real) #extract the summary 
SD_ca_sand_0_5_lm_real_slope <- SD_ca_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_sand_0.5_slopes[i] < SD_ca_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_sand_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.100.200)
  SD_ca_sand_100_200_lm_sum <- summary(SD_ca_sand_100_200_lm) #extracting the linear regression information
  SD_ca_sand_100_200_slopes <- c(SD_ca_sand_100_200_slopes, SD_ca_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_sand_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
SD_ca_sand_100_200_lm_real_sum <- summary(SD_ca_sand_100_200_lm_real) #extract the summary 
SD_ca_sand_100_200_lm_real_slope <- SD_ca_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_sand_100_200_slopes[i] < SD_ca_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_ph_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_0.5)
  SD_ca_ph_0_5_lm_sum <- summary(SD_ca_ph_0_5_lm) #extracting the linear regression information
  SD_ca_ph_0.5_slopes <- c(SD_ca_ph_0.5_slopes, SD_ca_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_ph_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
SD_ca_ph_0_5_lm_real_sum <- summary(SD_ca_ph_0_5_lm_real) #extract the summary 
SD_ca_ph_0_5_lm_real_slope <- SD_ca_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_ph_0.5_slopes[i] < SD_ca_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_ph_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_100.200)
  SD_ca_ph_100_200_lm_sum <- summary(SD_ca_ph_100_200_lm) #extracting the linear regression information
  SD_ca_ph_100_200_slopes <- c(SD_ca_ph_100_200_slopes, SD_ca_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_ph_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
SD_ca_ph_100_200_lm_real_sum <- summary(SD_ca_ph_100_200_lm_real) #extract the summary 
SD_ca_ph_100_200_lm_real_slope <- SD_ca_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_ph_100_200_slopes[i] < SD_ca_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_soc_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  SD_ca_soc_0_5_lm_sum <- summary(SD_ca_soc_0_5_lm) #extracting the linear regression information
  SD_ca_soc_0.5_slopes <- c(SD_ca_soc_0.5_slopes, SD_ca_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_soc_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
SD_ca_soc_0_5_lm_real_sum <- summary(SD_ca_soc_0_5_lm_real) #extract the summary 
SD_ca_soc_0_5_lm_real_slope <- SD_ca_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_soc_0.5_slopes[i] < SD_ca_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_soc_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  SD_ca_soc_100_200_lm_sum <- summary(SD_ca_soc_100_200_lm) #extracting the linear regression information
  SD_ca_soc_100_200_slopes <- c(SD_ca_soc_100_200_slopes, SD_ca_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_soc_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
SD_ca_soc_100_200_lm_real_sum <- summary(SD_ca_soc_100_200_lm_real) #extract the summary 
SD_ca_soc_100_200_lm_real_slope <- SD_ca_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_soc_100_200_slopes[i] < SD_ca_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_vol_10_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  SD_ca_vol_10_0_5_lm_sum <- summary(SD_ca_vol_10_0_5_lm) #extracting the linear regression information
  SD_ca_vol_10_0.5_slopes <- c(SD_ca_vol_10_0.5_slopes, SD_ca_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_vol_10_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
SD_ca_vol_10_0_5_lm_real_sum <- summary(SD_ca_vol_10_0_5_lm_real) #extract the summary 
SD_ca_vol_10_0_5_lm_real_slope <- SD_ca_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_vol_10_0.5_slopes[i] > SD_ca_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_vol_10_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  SD_ca_vol_10_100_200_lm_sum <- summary(SD_ca_vol_10_100_200_lm) #extracting the linear regression information
  SD_ca_vol_10_100_200_slopes <- c(SD_ca_vol_10_100_200_slopes, SD_ca_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_vol_10_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
SD_ca_vol_10_100_200_lm_real_sum <- summary(SD_ca_vol_10_100_200_lm_real) #extract the summary 
SD_ca_vol_10_100_200_lm_real_slope <- SD_ca_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_vol_10_100_200_slopes[i] < SD_ca_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_vol_1500_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  SD_ca_vol_1500_0_5_lm_sum <- summary(SD_ca_vol_1500_0_5_lm) #extracting the linear regression information
  SD_ca_vol_1500_0.5_slopes <- c(SD_ca_vol_1500_0.5_slopes, SD_ca_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_vol_1500_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
SD_ca_vol_1500_0_5_lm_real_sum <- summary(SD_ca_vol_1500_0_5_lm_real) #extract the summary 
SD_ca_vol_1500_0_5_lm_real_slope <- SD_ca_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_vol_1500_0.5_slopes[i] < SD_ca_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_vol_1500_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  SD_ca_vol_1500_100_200_lm_sum <- summary(SD_ca_vol_1500_100_200_lm) #extracting the linear regression information
  SD_ca_vol_1500_100_200_slopes <- c(SD_ca_vol_1500_100_200_slopes, SD_ca_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_vol_1500_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
SD_ca_vol_1500_100_200_lm_real_sum <- summary(SD_ca_vol_1500_100_200_lm_real) #extract the summary 
SD_ca_vol_1500_100_200_lm_real_slope <- SD_ca_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_vol_1500_100_200_slopes[i] > SD_ca_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_nitrogen_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  SD_ca_nitrogen_0_5_lm_sum <- summary(SD_ca_nitrogen_0_5_lm) #extracting the linear regression information
  SD_ca_nitrogen_0.5_slopes <- c(SD_ca_nitrogen_0.5_slopes, SD_ca_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_nitrogen_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
SD_ca_nitrogen_0_5_lm_real_sum <- summary(SD_ca_nitrogen_0_5_lm_real) #extract the summary 
SD_ca_nitrogen_0_5_lm_real_slope <- SD_ca_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_nitrogen_0.5_slopes[i] > SD_ca_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_ca_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Canopy_area.shuffled = sample(Canopy_area)) #create a data frame with a shuffled 
  SD_ca_nitrogen_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Canopy_area.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  SD_ca_nitrogen_100_200_lm_sum <- summary(SD_ca_nitrogen_100_200_lm) #extracting the linear regression information
  SD_ca_nitrogen_100_200_slopes <- c(SD_ca_nitrogen_100_200_slopes, SD_ca_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized ca values to the list of stored slopes
}

#extracting the slope of our points
SD_ca_nitrogen_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Canopy_area~SD_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
SD_ca_nitrogen_100_200_lm_real_sum <- summary(SD_ca_nitrogen_100_200_lm_real) #extract the summary 
SD_ca_nitrogen_100_200_lm_real_slope <- SD_ca_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_ca_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_ca_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled ca vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_ca_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_ca_nitrogen_100_200_slopes[i] > SD_ca_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_ca_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



# CS



# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_clay_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  SD_cs_clay_0_5_lm_sum <- summary(SD_cs_clay_0_5_lm) #extracting the linear regression information
  SD_cs_clay_0.5_slopes <- c(SD_cs_clay_0.5_slopes, SD_cs_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_clay_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
SD_cs_clay_0_5_lm_real_sum <- summary(SD_cs_clay_0_5_lm_real) #extract the summary 
SD_cs_clay_0_5_lm_real_slope <- SD_cs_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_clay_0.5_slopes[i] > SD_cs_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_clay_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  SD_cs_clay_100_200_lm_sum <- summary(SD_cs_clay_100_200_lm) #extracting the linear regression information
  SD_cs_clay_100_200_slopes <- c(SD_cs_clay_100_200_slopes, SD_cs_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_clay_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
SD_cs_clay_100_200_lm_real_sum <- summary(SD_cs_clay_100_200_lm_real) #extract the summary 
SD_cs_clay_100_200_lm_real_slope <- SD_cs_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_clay_100_200_slopes[i] > SD_cs_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_silt_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.0.5)
  SD_cs_silt_0_5_lm_sum <- summary(SD_cs_silt_0_5_lm) #extracting the linear regression information
  SD_cs_silt_0.5_slopes <- c(SD_cs_silt_0.5_slopes, SD_cs_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_silt_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
SD_cs_silt_0_5_lm_real_sum <- summary(SD_cs_silt_0_5_lm_real) #extract the summary 
SD_cs_silt_0_5_lm_real_slope <- SD_cs_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_silt_0.5_slopes[i] > SD_cs_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_silt_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.100.200)
  SD_cs_silt_100_200_lm_sum <- summary(SD_cs_silt_100_200_lm) #extracting the linear regression information
  SD_cs_silt_100_200_slopes <- c(SD_cs_silt_100_200_slopes, SD_cs_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_silt_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
SD_cs_silt_100_200_lm_real_sum <- summary(SD_cs_silt_100_200_lm_real) #extract the summary 
SD_cs_silt_100_200_lm_real_slope <- SD_cs_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_silt_100_200_slopes[i] > SD_cs_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_sand_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.0.5)
  SD_cs_sand_0_5_lm_sum <- summary(SD_cs_sand_0_5_lm) #extracting the linear regression information
  SD_cs_sand_0.5_slopes <- c(SD_cs_sand_0.5_slopes, SD_cs_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_sand_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
SD_cs_sand_0_5_lm_real_sum <- summary(SD_cs_sand_0_5_lm_real) #extract the summary 
SD_cs_sand_0_5_lm_real_slope <- SD_cs_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_sand_0.5_slopes[i] < SD_cs_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_sand_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.100.200)
  SD_cs_sand_100_200_lm_sum <- summary(SD_cs_sand_100_200_lm) #extracting the linear regression information
  SD_cs_sand_100_200_slopes <- c(SD_cs_sand_100_200_slopes, SD_cs_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_sand_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
SD_cs_sand_100_200_lm_real_sum <- summary(SD_cs_sand_100_200_lm_real) #extract the summary 
SD_cs_sand_100_200_lm_real_slope <- SD_cs_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_sand_100_200_slopes[i] < SD_cs_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_ph_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_0.5)
  SD_cs_ph_0_5_lm_sum <- summary(SD_cs_ph_0_5_lm) #extracting the linear regression information
  SD_cs_ph_0.5_slopes <- c(SD_cs_ph_0.5_slopes, SD_cs_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_ph_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
SD_cs_ph_0_5_lm_real_sum <- summary(SD_cs_ph_0_5_lm_real) #extract the summary 
SD_cs_ph_0_5_lm_real_slope <- SD_cs_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_ph_0.5_slopes[i] < SD_cs_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_ph_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_100.200)
  SD_cs_ph_100_200_lm_sum <- summary(SD_cs_ph_100_200_lm) #extracting the linear regression information
  SD_cs_ph_100_200_slopes <- c(SD_cs_ph_100_200_slopes, SD_cs_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_ph_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
SD_cs_ph_100_200_lm_real_sum <- summary(SD_cs_ph_100_200_lm_real) #extract the summary 
SD_cs_ph_100_200_lm_real_slope <- SD_cs_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_ph_100_200_slopes[i] > SD_cs_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_soc_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  SD_cs_soc_0_5_lm_sum <- summary(SD_cs_soc_0_5_lm) #extracting the linear regression information
  SD_cs_soc_0.5_slopes <- c(SD_cs_soc_0.5_slopes, SD_cs_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_soc_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
SD_cs_soc_0_5_lm_real_sum <- summary(SD_cs_soc_0_5_lm_real) #extract the summary 
SD_cs_soc_0_5_lm_real_slope <- SD_cs_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_soc_0.5_slopes[i] < SD_cs_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_soc_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  SD_cs_soc_100_200_lm_sum <- summary(SD_cs_soc_100_200_lm) #extracting the linear regression information
  SD_cs_soc_100_200_slopes <- c(SD_cs_soc_100_200_slopes, SD_cs_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_soc_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
SD_cs_soc_100_200_lm_real_sum <- summary(SD_cs_soc_100_200_lm_real) #extract the summary 
SD_cs_soc_100_200_lm_real_slope <- SD_cs_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_soc_100_200_slopes[i] < SD_cs_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_vol_10_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  SD_cs_vol_10_0_5_lm_sum <- summary(SD_cs_vol_10_0_5_lm) #extracting the linear regression information
  SD_cs_vol_10_0.5_slopes <- c(SD_cs_vol_10_0.5_slopes, SD_cs_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_vol_10_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
SD_cs_vol_10_0_5_lm_real_sum <- summary(SD_cs_vol_10_0_5_lm_real) #extract the summary 
SD_cs_vol_10_0_5_lm_real_slope <- SD_cs_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_vol_10_0.5_slopes[i] > SD_cs_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_vol_10_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  SD_cs_vol_10_100_200_lm_sum <- summary(SD_cs_vol_10_100_200_lm) #extracting the linear regression information
  SD_cs_vol_10_100_200_slopes <- c(SD_cs_vol_10_100_200_slopes, SD_cs_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_vol_10_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
SD_cs_vol_10_100_200_lm_real_sum <- summary(SD_cs_vol_10_100_200_lm_real) #extract the summary 
SD_cs_vol_10_100_200_lm_real_slope <- SD_cs_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_vol_10_100_200_slopes[i] < SD_cs_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_vol_1500_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  SD_cs_vol_1500_0_5_lm_sum <- summary(SD_cs_vol_1500_0_5_lm) #extracting the linear regression information
  SD_cs_vol_1500_0.5_slopes <- c(SD_cs_vol_1500_0.5_slopes, SD_cs_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_vol_1500_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
SD_cs_vol_1500_0_5_lm_real_sum <- summary(SD_cs_vol_1500_0_5_lm_real) #extract the summary 
SD_cs_vol_1500_0_5_lm_real_slope <- SD_cs_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_vol_1500_0.5_slopes[i] < SD_cs_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_vol_1500_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  SD_cs_vol_1500_100_200_lm_sum <- summary(SD_cs_vol_1500_100_200_lm) #extracting the linear regression information
  SD_cs_vol_1500_100_200_slopes <- c(SD_cs_vol_1500_100_200_slopes, SD_cs_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_vol_1500_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
SD_cs_vol_1500_100_200_lm_real_sum <- summary(SD_cs_vol_1500_100_200_lm_real) #extract the summary 
SD_cs_vol_1500_100_200_lm_real_slope <- SD_cs_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_vol_1500_100_200_slopes[i] > SD_cs_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_nitrogen_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  SD_cs_nitrogen_0_5_lm_sum <- summary(SD_cs_nitrogen_0_5_lm) #extracting the linear regression information
  SD_cs_nitrogen_0.5_slopes <- c(SD_cs_nitrogen_0.5_slopes, SD_cs_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_nitrogen_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
SD_cs_nitrogen_0_5_lm_real_sum <- summary(SD_cs_nitrogen_0_5_lm_real) #extract the summary 
SD_cs_nitrogen_0_5_lm_real_slope <- SD_cs_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_nitrogen_0.5_slopes[i] > SD_cs_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_cs_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, Crown_spread.shuffled = sample(Crown_spread)) #create a data frame with a shuffled 
  SD_cs_nitrogen_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$Crown_spread.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  SD_cs_nitrogen_100_200_lm_sum <- summary(SD_cs_nitrogen_100_200_lm) #extracting the linear regression information
  SD_cs_nitrogen_100_200_slopes <- c(SD_cs_nitrogen_100_200_slopes, SD_cs_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized cs values to the list of stored slopes
}

#extracting the slope of our points
SD_cs_nitrogen_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$Crown_spread~SD_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
SD_cs_nitrogen_100_200_lm_real_sum <- summary(SD_cs_nitrogen_100_200_lm_real) #extract the summary 
SD_cs_nitrogen_100_200_lm_real_slope <- SD_cs_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_cs_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_cs_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled cs vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_cs_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_cs_nitrogen_100_200_slopes[i] < SD_cs_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_cs_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN



#DBH_ag


# Clay Content 0-5 cm

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_clay_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_clay_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.0.5)
  SD_dbh_clay_0_5_lm_sum <- summary(SD_dbh_clay_0_5_lm) #extracting the linear regression information
  SD_dbh_clay_0.5_slopes <- c(SD_dbh_clay_0.5_slopes, SD_dbh_clay_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_clay_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$clay.content.0.5) #creating the linear regression
SD_dbh_clay_0_5_lm_real_sum <- summary(SD_dbh_clay_0_5_lm_real) #extract the summary 
SD_dbh_clay_0_5_lm_real_slope <- SD_dbh_clay_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_clay_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_clay_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_clay_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_clay_0.5_slopes[i] > SD_dbh_clay_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_clay_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

# Clay Content 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_clay_100_200_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_clay_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$clay.content.100.200)
  SD_dbh_clay_100_200_lm_sum <- summary(SD_dbh_clay_100_200_lm) #extracting the linear regression information
  SD_dbh_clay_100_200_slopes <- c(SD_dbh_clay_100_200_slopes, SD_dbh_clay_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_clay_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$clay.content.100.200) #creating the linear regression
SD_dbh_clay_100_200_lm_real_sum <- summary(SD_dbh_clay_100_200_lm_real) #extract the summary 
SD_dbh_clay_100_200_lm_real_slope <- SD_dbh_clay_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_clay_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_clay_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Clay Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_clay_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_clay_100_200_slopes[i] > SD_dbh_clay_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_clay_100_200_slopes)) #the proportion of random ANNs that are less than our ANN

#silt 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_silt_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_silt_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.0.5)
  SD_dbh_silt_0_5_lm_sum <- summary(SD_dbh_silt_0_5_lm) #extracting the linear regression information
  SD_dbh_silt_0.5_slopes <- c(SD_dbh_silt_0.5_slopes, SD_dbh_silt_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_silt_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$silt.0.5) #creating the linear regression
SD_dbh_silt_0_5_lm_real_sum <- summary(SD_dbh_silt_0_5_lm_real) #extract the summary 
SD_dbh_silt_0_5_lm_real_slope <- SD_dbh_silt_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_silt_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_silt_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_silt_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_silt_0.5_slopes[i] > SD_dbh_silt_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_silt_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#silt 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_silt_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_silt_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$silt.100.200)
  SD_dbh_silt_100_200_lm_sum <- summary(SD_dbh_silt_100_200_lm) #extracting the linear regression information
  SD_dbh_silt_100_200_slopes <- c(SD_dbh_silt_100_200_slopes, SD_dbh_silt_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_silt_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$silt.100.200) #creating the linear regression
SD_dbh_silt_100_200_lm_real_sum <- summary(SD_dbh_silt_100_200_lm_real) #extract the summary 
SD_dbh_silt_100_200_lm_real_slope <- SD_dbh_silt_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_silt_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_silt_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Silt Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_silt_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_silt_100_200_slopes[i] > SD_dbh_silt_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_silt_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#sand 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_sand_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_sand_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.0.5)
  SD_dbh_sand_0_5_lm_sum <- summary(SD_dbh_sand_0_5_lm) #extracting the linear regression information
  SD_dbh_sand_0.5_slopes <- c(SD_dbh_sand_0.5_slopes, SD_dbh_sand_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_sand_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$sand.0.5) #creating the linear regression
SD_dbh_sand_0_5_lm_real_sum <- summary(SD_dbh_sand_0_5_lm_real) #extract the summary 
SD_dbh_sand_0_5_lm_real_slope <- SD_dbh_sand_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_sand_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_sand_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_sand_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_sand_0.5_slopes[i] < SD_dbh_sand_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_sand_0.5_slopes)) #the proportion of random ANNs that are less than our ANN

#sand 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_sand_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_sand_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$sand.100.200)
  SD_dbh_sand_100_200_lm_sum <- summary(SD_dbh_sand_100_200_lm) #extracting the linear regression information
  SD_dbh_sand_100_200_slopes <- c(SD_dbh_sand_100_200_slopes, SD_dbh_sand_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_sand_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$sand.100.200) #creating the linear regression
SD_dbh_sand_100_200_lm_real_sum <- summary(SD_dbh_sand_100_200_lm_real) #extract the summary 
SD_dbh_sand_100_200_lm_real_slope <- SD_dbh_sand_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_sand_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_sand_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Sand Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_sand_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_sand_100_200_slopes[i] < SD_dbh_sand_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_sand_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


##ph 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_ph_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_ph_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_0.5)
  SD_dbh_ph_0_5_lm_sum <- summary(SD_dbh_ph_0_5_lm) #extracting the linear regression information
  SD_dbh_ph_0.5_slopes <- c(SD_dbh_ph_0.5_slopes, SD_dbh_ph_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_ph_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$ph_0.5) #creating the linear regression
SD_dbh_ph_0_5_lm_real_sum <- summary(SD_dbh_ph_0_5_lm_real) #extract the summary 
SD_dbh_ph_0_5_lm_real_slope <- SD_dbh_ph_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_ph_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_ph_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_ph_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_ph_0.5_slopes[i] > SD_dbh_ph_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_ph_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#ph 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_ph_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils, DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_ph_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$ph_100.200)
  SD_dbh_ph_100_200_lm_sum <- summary(SD_dbh_ph_100_200_lm) #extracting the linear regression information
  SD_dbh_ph_100_200_slopes <- c(SD_dbh_ph_100_200_slopes, SD_dbh_ph_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_ph_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$ph_100.200) #creating the linear regression
SD_dbh_ph_100_200_lm_real_sum <- summary(SD_dbh_ph_100_200_lm_real) #extract the summary 
SD_dbh_ph_100_200_lm_real_slope <- SD_dbh_ph_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_ph_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_ph_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. ph Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_ph_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_ph_100_200_slopes[i] > SD_dbh_ph_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_ph_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_soc_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_soc_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.0.5)
  SD_dbh_soc_0_5_lm_sum <- summary(SD_dbh_soc_0_5_lm) #extracting the linear regression information
  SD_dbh_soc_0.5_slopes <- c(SD_dbh_soc_0.5_slopes, SD_dbh_soc_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_soc_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$SOC.0.5) #creating the linear regression
SD_dbh_soc_0_5_lm_real_sum <- summary(SD_dbh_soc_0_5_lm_real) #extract the summary 
SD_dbh_soc_0_5_lm_real_slope <- SD_dbh_soc_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_soc_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_soc_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_soc_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_soc_0.5_slopes[i] > SD_dbh_soc_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_soc_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#soil organic carbon 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_soc_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_soc_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$SOC.100.200)
  SD_dbh_soc_100_200_lm_sum <- summary(SD_dbh_soc_100_200_lm) #extracting the linear regression information
  SD_dbh_soc_100_200_slopes <- c(SD_dbh_soc_100_200_slopes, SD_dbh_soc_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_soc_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$SOC.100.200) #creating the linear regression
SD_dbh_soc_100_200_lm_real_sum <- summary(SD_dbh_soc_100_200_lm_real) #extract the summary 
SD_dbh_soc_100_200_lm_real_slope <- SD_dbh_soc_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_soc_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_soc_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Soil Organic Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_soc_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_soc_100_200_slopes[i] > SD_dbh_soc_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_soc_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 0-5

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_vol_10_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_vol_10_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_0.5)
  SD_dbh_vol_10_0_5_lm_sum <- summary(SD_dbh_vol_10_0_5_lm) #extracting the linear regression information
  SD_dbh_vol_10_0.5_slopes <- c(SD_dbh_vol_10_0.5_slopes, SD_dbh_vol_10_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_vol_10_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$vol_water_.10_0.5) #creating the linear regression
SD_dbh_vol_10_0_5_lm_real_sum <- summary(SD_dbh_vol_10_0_5_lm_real) #extract the summary 
SD_dbh_vol_10_0_5_lm_real_slope <- SD_dbh_vol_10_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_vol_10_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_vol_10_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_vol_10_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_vol_10_0.5_slopes[i] > SD_dbh_vol_10_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_vol_10_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -10 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_vol_10_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_vol_10_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.10_100.200)
  SD_dbh_vol_10_100_200_lm_sum <- summary(SD_dbh_vol_10_100_200_lm) #extracting the linear regression information
  SD_dbh_vol_10_100_200_slopes <- c(SD_dbh_vol_10_100_200_slopes, SD_dbh_vol_10_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_vol_10_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$vol_water_.10_100.200) #creating the linear regression
SD_dbh_vol_10_100_200_lm_real_sum <- summary(SD_dbh_vol_10_100_200_lm_real) #extract the summary 
SD_dbh_vol_10_100_200_lm_real_slope <- SD_dbh_vol_10_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_vol_10_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_vol_10_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -10 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_vol_10_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_vol_10_100_200_slopes[i] < SD_dbh_vol_10_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_vol_10_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 0-5


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_vol_1500_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_vol_1500_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500kPa_0.5)
  SD_dbh_vol_1500_0_5_lm_sum <- summary(SD_dbh_vol_1500_0_5_lm) #extracting the linear regression information
  SD_dbh_vol_1500_0.5_slopes <- c(SD_dbh_vol_1500_0.5_slopes, SD_dbh_vol_1500_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_vol_1500_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$vol_water_.1500kPa_0.5) #creating the linear regression
SD_dbh_vol_1500_0_5_lm_real_sum <- summary(SD_dbh_vol_1500_0_5_lm_real) #extract the summary 
SD_dbh_vol_1500_0_5_lm_real_slope <- SD_dbh_vol_1500_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_vol_1500_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_vol_1500_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_vol_1500_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_vol_1500_0.5_slopes[i] < SD_dbh_vol_1500_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_vol_1500_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#volume of water content at -1500 kpa 100-200


#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_vol_1500_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_vol_1500_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$vol_water_.1500_100.200)
  SD_dbh_vol_1500_100_200_lm_sum <- summary(SD_dbh_vol_1500_100_200_lm) #extracting the linear regression information
  SD_dbh_vol_1500_100_200_slopes <- c(SD_dbh_vol_1500_100_200_slopes, SD_dbh_vol_1500_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_vol_1500_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$vol_water_.1500_100.200) #creating the linear regression
SD_dbh_vol_1500_100_200_lm_real_sum <- summary(SD_dbh_vol_1500_100_200_lm_real) #extract the summary 
SD_dbh_vol_1500_100_200_lm_real_slope <- SD_dbh_vol_1500_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_vol_1500_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_vol_1500_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Volume of Water at -1500 kPa 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_vol_1500_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_vol_1500_100_200_slopes[i] > SD_dbh_vol_1500_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_vol_1500_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 05

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_nitrogen_0.5_slopes <- c() #creating empty list to collect p values 

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_nitrogen_0_5_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.0.5)
  SD_dbh_nitrogen_0_5_lm_sum <- summary(SD_dbh_nitrogen_0_5_lm) #extracting the linear regression information
  SD_dbh_nitrogen_0.5_slopes <- c(SD_dbh_nitrogen_0.5_slopes, SD_dbh_nitrogen_0_5_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_nitrogen_0_5_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$nitrogen.0.5) #creating the linear regression
SD_dbh_nitrogen_0_5_lm_real_sum <- summary(SD_dbh_nitrogen_0_5_lm_real) #extract the summary 
SD_dbh_nitrogen_0_5_lm_real_slope <- SD_dbh_nitrogen_0_5_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_nitrogen_0.5_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_nitrogen_0_5_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 0-5 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_nitrogen_0.5_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_nitrogen_0.5_slopes[i] > SD_dbh_nitrogen_0_5_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_nitrogen_0.5_slopes)) #the proportion of random ANNs that are less than our ANN


#nitrogen 100-200

#extracting slopes from comparing soil values with randomized shape/size values with linear regressions
SD_dbh_nitrogen_100_200_slopes <- c() #creating empty list to collect p values

set.seed(21)
for (i in 1:1000){ #for 1000 permutations
  SD_fixed_field_data_processed_soils_shuffled <- transform(SD_fixed_field_data_processed_soils,DBH_ag.shuffled = sample(DBH_ag)) #create a data frame with a shuffled 
  SD_dbh_nitrogen_100_200_lm <- lm(SD_fixed_field_data_processed_soils_shuffled$DBH_ag.shuffled~SD_fixed_field_data_processed_soils_shuffled$nitrogen.100.200)
  SD_dbh_nitrogen_100_200_lm_sum <- summary(SD_dbh_nitrogen_100_200_lm) #extracting the linear regression information
  SD_dbh_nitrogen_100_200_slopes <- c(SD_dbh_nitrogen_100_200_slopes, SD_dbh_nitrogen_100_200_lm_sum$coefficients[2]) #add the current p-value from the randomized dbh values to the list of stored slopes
}

#extracting the slope of our points
SD_dbh_nitrogen_100_200_lm_real <- lm(SD_fixed_field_data_processed_soils$DBH_ag~SD_fixed_field_data_processed_soils$nitrogen.100.200) #creating the linear regression
SD_dbh_nitrogen_100_200_lm_real_sum <- summary(SD_dbh_nitrogen_100_200_lm_real) #extract the summary 
SD_dbh_nitrogen_100_200_lm_real_slope <- SD_dbh_nitrogen_100_200_lm_real_sum$coefficients[2] #storing the slope

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=SD_dbh_nitrogen_100_200_slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=SD_dbh_nitrogen_100_200_lm_real_slope, col = "red")+ #line of our real slope
  xlab("Slopes of Shuffled dbh vs. Nitrogen Content 100-200 cm")+
  theme_classic()


#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(SD_dbh_nitrogen_100_200_slopes)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (SD_dbh_nitrogen_100_200_slopes[i] > SD_dbh_nitrogen_100_200_lm_real_slope){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(SD_dbh_nitrogen_100_200_slopes)) #the proportion of random ANNs that are less than our ANN


#Summarizing results across population, size/shape, and soil variables

#loading in all of the results from comparing the real slopes to randomized slopes 
Size_Slope_Results.df <- read.csv("./data/Hypothesis_4_Size_Shape_Test_Results_Wanger_REU_2024.csv") #imports the csv created from analyzing_morpho_data_cleaned.R
View(Size_Slope_Results.df)

#summarizing the data frame 
summary(Size_Slope_Results.df)

#number of yes or no by population 
table(Size_Slope_Results.df$Population, Size_Slope_Results.df$Sig.)

# number of yes or nos for each size variables 

#LM
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="LM"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LM"])

#LC
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="LC"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LC"])

#SD
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="SD"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="SD"])

# number of yes or nos for each soil variables 

#LM
table(Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="LM"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LM"])

#LC
table(Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="LC"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LC"])

#SD
table(Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="SD"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="SD"])


## number of yes or nos for each size variables and each soil variable 

#LM 
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="LM"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LM"], Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="LM"])

#LC 
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="LC"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="LC"], Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="LC"])

#SD 
table(Size_Slope_Results.df$Shape.Size.Variable[Size_Slope_Results.df$Population=="SD"], Size_Slope_Results.df$Sig.[Size_Slope_Results.df$Population=="SD"], Size_Slope_Results.df$Variable[Size_Slope_Results.df$Population=="SD"])


#plots to visualize the results

#plot of population, yes/no, and the interaciton between size/shape and soil values
ggplot(Size_Slope_Results.df, aes(x = Sig., fill = interaction(Shape.Size.Variable,Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and  size/shape 
ggplot(Size_Slope_Results.df, aes(x = Sig., fill = interaction(Shape.Size.Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and soil values
ggplot(Size_Slope_Results.df, aes(x = Sig., fill = interaction(Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()



### PART 3: Comparing average soil values from inside populations to outside populations ###

#downloading the locations of the 16 known populations
all_pop_locations.df <- read.csv(file = "./data/All QUBR sites- all sites.csv")
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

plot(BCS_polygon_UTM)

#cropping the BCS polygon to just be the southern region 
BCS_polygon_box <- st_bbox(BCS_polygon_UTM)


#creating a cropped bbox 
BCS_polygon_box_sf <- BCS_polygon_box %>% #turning the bbox into polygon
  st_as_sfc()
BCS_polygon_box_spatial <- as(BCS_polygon_box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
BCS_polygon_box_spatial_cropped <- raster::crop(BCS_polygon_box_spatial, extent((BCS_polygon_box[1]+411199), (BCS_polygon_box[3]), (BCS_polygon_box[2]), (BCS_polygon_box[4]-432999))) #cropping the xmin, xmax, ymin, and ymax by 20 m inside
BCS_polygon_box_sf_cropped <-  BCS_polygon_box_spatial_cropped %>% #turning the spatial polygon into a polygon
  st_as_sfc()

#cropping the points by the cropped box
BCS_polygon_UTM_cropped<- st_crop(BCS_polygon_UTM, BCS_polygon_box_sf_cropped) #cropping the BCS polygon with the cropped bbox extent

BCS_polygon_UTM_cropped <- BCS_polygon_UTM_cropped %>% #turning the cropped BCS polygon into a sf file
  st_as_sf

#plotting the original extent and polygon and the cropped extent
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)
  
#creating a bbox around the cropped BCS file
BCS_polygon_crop_bbox <- st_bbox(BCS_polygon_UTM_cropped) 


#Load in extra environmental rasters because the lower portion of BCS needs to be covered by the 250 m soil rasters


#loading in soil textures from CONABIO, theses are too larger, about 1 km^2 I believe
clay_05_2 <- raster(paste0("./data/Soil Grid/clay content/clay content 0-5 part 2.tif"))
clay_200_2 <- raster(paste0("./data/Soil Grid/clay content/clay contentn 100-200 part 2.tif"))
silt_05_2 <- raster(paste0("./data/Soil Grid/silt/silt 0-5 part 2.tif"))
silt_200_2 <-raster(paste0("./data/Soil Grid/silt/silt 100-200 part 2.tif"))
sand_05_2 <- raster(paste0("./data/Soil Grid/sand/sand 0-5 part 2.tif"))
sand_200_2 <- raster(paste0("./data/Soil Grid/sand/sand 100-200 part 2.tif"))
ph_05_2 <- raster(paste0("./data/Soil Grid/pH/ph_0-5 part 2.tif")) #0-5 cm ph
ph_200_2 <- raster(paste0("./data/Soil Grid/pH/ph_100-200 part 2.tif")) #100-200 ph
ocd_05_2 <- raster(paste0("./data/Soil Grid/organic carbon density/OCD_0-5 part 2.tif")) #0-5cm organic carbon density
ocd_200_2 <- raster(paste0("./data/Soil Grid/organic carbon density/OCD_100-200 part 2.tif")) #100-200cm organic carbon density
coarse_frag_05_2 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_0-5_part_2.tif")) #0-5 cm coarse fragments
coarse_frag_200_2 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_100-200_part_2.tif")) #100-200 cm coarse fragments
cat_ex_cap_05_2 <-raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_0-5_part_2.tif")) #0-5 cm cation exchange capacity
cat_ex_cap_200_2 <- raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_100-200_part_2.tif")) #100-200 cm cation exchange capacity
bulk_dens_05_2 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_0-5_part_2.tif")) #0-5 cm bulk density
bulk_dens_200_2 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_100_200_part_2.tif")) #100-200 cm bulk density
vol_wat_10kpa_05_2 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_0-5  part 2.tif"))  #0-5 cm -10 kpa volumn water content
vol_wat_10kpa_200_2 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_100-200 part 2.tif"))  #100-200 cm -10 kpa volumn water content
vol_wat_1500kpa_05_2 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500_0-5 part 2.tif"))  #0-5 cm -1500 kpa volumn water content
vol_wat_1500kpa_200_2 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500_100-200 part 2.tif")) #100-200 cm -1500 kpa volumn water content
nitrogen_05_2 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 0-5 part 2.tif"))
nitrogen_200_2 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 100-200 part 2.tif"))
Soil_Organic_Carbon_05_2 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 0-5 part 2.tif"))
Soil_Organic_Carbon_200_2 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 100-200 part 2.tif"))



#project rasters to equal area projection (UTM 12N), uses meters as distance measurement 
clay_05_utm_2 <- projectRaster(clay_05_2, crs=26912) #converting the 0-5 cm clay raster to utm 12
clay_200_utm_2 <- projectRaster(clay_200_2, crs=26912) #converting the 90-200 cm clay raster to utm 12
silt_05_utm_2 <- projectRaster(silt_05_2, crs=26912)
silt_200_utm_2 <- projectRaster(silt_200_2, crs=26912)
sand_05_utm_2 <- projectRaster(sand_05_2, crs=26912)
sand_200_utm_2 <- projectRaster(sand_200_2, crs=26912)
ph_05_utm_2 <- projectRaster(ph_05_2, crs=26912) 
ph_200_utm_2 <- projectRaster(ph_200_2, crs=26912) 
ocd_05_utm_2 <- projectRaster(ocd_05_2, crs=26912)
ocd_200_utm_2 <- projectRaster(ocd_200_2, crs=26912)
coarse_frag_05_utm_2 <- projectRaster(coarse_frag_05_2, crs=26912)
coarse_frag_200_utm_2 <- projectRaster(coarse_frag_200_2, crs=26912)
cat_ex_cap_05_utm_2 <- projectRaster(cat_ex_cap_05_2, crs=26912)
cat_ex_cap_200_utm_2 <- projectRaster(cat_ex_cap_200_2, crs=26912)
bulk_dens_05_utm_2 <- projectRaster(bulk_dens_05_2, crs=26912)
bulk_dens_200_utm_2 <- projectRaster(bulk_dens_200_2, crs=26912)
vol_wat_10kpa_05_utm_2 <- projectRaster(vol_wat_10kpa_05_2, crs=26912)
vol_wat_10kpa_200_utm_2 <- projectRaster(vol_wat_10kpa_200_2, crs=26912)
vol_wat_1500kpa_05_utm_2 <- projectRaster(vol_wat_1500kpa_05_2, crs=26912)
vol_wat_1500kpa_200_utm_2 <- projectRaster(vol_wat_1500kpa_200_2, crs=26912)
nitrogen_05_utm_2 <- projectRaster(nitrogen_05_2, crs=26912)
nitrogen_200_utm_2 <- projectRaster(nitrogen_200_2, crs=26912)
Soil_Organic_Carbon_05_utm_2 <- projectRaster(Soil_Organic_Carbon_05_2, crs=26912)
Soil_Organic_Carbon_200_utm_2 <- projectRaster(Soil_Organic_Carbon_200_2, crs=26912)


#cropping the new rasters so they do not overlap with the older ones

clay_05_utm_2_bbox <- clay_05_utm_2 %>%
  st_bbox()
clay_05_utm_2_crop <- crop(clay_05_utm_2, extent((clay_05_utm_2_bbox[1]), (clay_05_utm_2_bbox[3]), (clay_05_utm_2_bbox[2]), (2538644)))

clay_200_utm_2_bbox <- clay_200_utm_2 %>%
  st_bbox()
clay_200_utm_2_crop <- crop(clay_200_utm_2, extent((clay_200_utm_2_bbox[1]), (clay_200_utm_2_bbox[3]), (clay_200_utm_2_bbox[2]), (2538644)))

silt_05_utm_2_bbox <- silt_05_utm_2 %>%
  st_bbox()
silt_05_utm_2_crop <- crop(silt_05_utm_2, extent((silt_05_utm_2_bbox[1]), (silt_05_utm_2_bbox[3]), (silt_05_utm_2_bbox[2]), (2538644)))

silt_200_utm_2_bbox <- silt_200_utm_2 %>%
  st_bbox()
silt_200_utm_2_crop <- crop(silt_200_utm_2, extent((silt_200_utm_2_bbox[1]), (silt_200_utm_2_bbox[3]), (silt_200_utm_2_bbox[2]), (2538644)))

sand_05_utm_2_bbox <- sand_05_utm_2 %>%
  st_bbox()
sand_05_utm_2_crop <- crop(sand_05_utm_2, extent((sand_05_utm_2_bbox[1]), (sand_05_utm_2_bbox[3]), (sand_05_utm_2_bbox[2]), (2538644)))

sand_200_utm_2_bbox <- sand_200_utm_2 %>%
  st_bbox()
sand_200_utm_2_crop <- crop(sand_200_utm_2, extent((sand_200_utm_2_bbox[1]), (sand_200_utm_2_bbox[3]), (sand_200_utm_2_bbox[2]), (2538644)))

ph_05_utm_2_bbox <- ph_05_utm_2 %>%
  st_bbox()
ph_05_utm_2_crop <- crop(ph_05_utm_2, extent((ph_05_utm_2_bbox[1]), (ph_05_utm_2_bbox[3]), (ph_05_utm_2_bbox[2]), (2538644)))

ph_200_utm_2_bbox <- ph_200_utm_2 %>%
  st_bbox()
ph_200_utm_2_crop <- crop(ph_200_utm_2, extent((ph_200_utm_2_bbox[1]), (ph_200_utm_2_bbox[3]), (ph_200_utm_2_bbox[2]), (2538644)))

ocd_05_utm_2_bbox <- ocd_05_utm_2 %>%
  st_bbox()
ocd_05_utm_2_crop <- crop(ocd_05_utm_2, extent((ocd_05_utm_2_bbox[1]), (ocd_05_utm_2_bbox[3]), (ocd_05_utm_2_bbox[2]), (2538644)))

ocd_200_utm_2_bbox <- ocd_200_utm_2 %>%
  st_bbox()
ocd_200_utm_2_crop <- crop(ocd_200_utm_2, extent((ocd_200_utm_2_bbox[1]), (ocd_200_utm_2_bbox[3]), (ocd_200_utm_2_bbox[2]), (2538644)))

coarse_frag_05_utm_2_bbox <- coarse_frag_05_utm_2 %>%
  st_bbox()
coarse_frag_05_utm_2_crop <- crop(coarse_frag_05_utm_2, extent((coarse_frag_05_utm_2_bbox[1]), (coarse_frag_05_utm_2_bbox[3]), (coarse_frag_05_utm_2_bbox[2]), (2538644)))

coarse_frag_200_utm_2_bbox <- coarse_frag_200_utm_2 %>%
  st_bbox()
coarse_frag_200_utm_2_crop <- crop(coarse_frag_200_utm_2, extent((coarse_frag_200_utm_2_bbox[1]), (coarse_frag_200_utm_2_bbox[3]), (coarse_frag_200_utm_2_bbox[2]), (2538644)))

cat_ex_cap_05_utm_2_bbox <- cat_ex_cap_05_utm_2 %>%
  st_bbox()
cat_ex_cap_05_utm_2_crop <- crop(cat_ex_cap_05_utm_2, extent((cat_ex_cap_05_utm_2_bbox[1]), (cat_ex_cap_05_utm_2_bbox[3]), (cat_ex_cap_05_utm_2_bbox[2]), (2538644)))

cat_ex_cap_200_utm_2_bbox <- cat_ex_cap_05_utm_2 %>%
  st_bbox()
cat_ex_cap_200_utm_2_crop <- crop(cat_ex_cap_05_utm_2, extent((cat_ex_cap_05_utm_2_bbox[1]), (cat_ex_cap_05_utm_2_bbox[3]), (cat_ex_cap_05_utm_2_bbox[2]), (2538644)))

bulk_dens_05_utm_2_bbox <- bulk_dens_05_utm_2 %>%
  st_bbox()
bulk_dens_05_utm_2_crop <- crop(bulk_dens_05_utm_2, extent((bulk_dens_05_utm_2_bbox[1]), (bulk_dens_05_utm_2_bbox[3]), (bulk_dens_05_utm_2_bbox[2]), (2538644)))

bulk_dens_200_utm_2_bbox <- bulk_dens_200_utm_2 %>%
  st_bbox()
bulk_dens_200_utm_2_crop <- crop(bulk_dens_200_utm_2, extent((bulk_dens_200_utm_2_bbox[1]), (bulk_dens_200_utm_2_bbox[3]), (bulk_dens_200_utm_2_bbox[2]), (2538644)))

vol_wat_10kpa_05_utm_2_bbox <- vol_wat_10kpa_05_utm_2 %>%
  st_bbox()
vol_wat_10kpa_05_utm_2_crop <- crop(vol_wat_10kpa_05_utm_2, extent((vol_wat_10kpa_05_utm_2_bbox[1]), (vol_wat_10kpa_05_utm_2_bbox[3]), (vol_wat_10kpa_05_utm_2_bbox[2]), (2538644)))

vol_wat_10kpa_200_utm_2_bbox <- vol_wat_10kpa_200_utm_2 %>%
  st_bbox()
vol_wat_10kpa_200_utm_2_crop <- crop(vol_wat_10kpa_200_utm_2, extent((vol_wat_10kpa_200_utm_2_bbox[1]), (vol_wat_10kpa_200_utm_2_bbox[3]), (vol_wat_10kpa_200_utm_2_bbox[2]), (2538644)))

vol_wat_1500kpa_05_utm_2_bbox <- vol_wat_1500kpa_05_utm_2 %>%
  st_bbox()
vol_wat_1500kpa_05_utm_2_crop <- crop(vol_wat_1500kpa_05_utm_2, extent((vol_wat_1500kpa_05_utm_2_bbox[1]), (vol_wat_1500kpa_05_utm_2_bbox[3]), (vol_wat_1500kpa_05_utm_2_bbox[2]), (2538644)))

vol_wat_1500kpa_200_utm_2_bbox <- vol_wat_1500kpa_200_utm_2 %>%
  st_bbox()
vol_wat_1500kpa_200_utm_2_crop <- crop(vol_wat_1500kpa_200_utm_2, extent((vol_wat_1500kpa_200_utm_2_bbox[1]), (vol_wat_1500kpa_200_utm_2_bbox[3]), (vol_wat_1500kpa_200_utm_2_bbox[2]), (2538644)))

nitrogen_05_utm_2_bbox <- nitrogen_05_utm_2 %>%
  st_bbox()
nitrogen_05_utm_2_crop <- crop(nitrogen_05_utm_2, extent((nitrogen_05_utm_2_bbox[1]), (nitrogen_05_utm_2_bbox[3]), (nitrogen_05_utm_2_bbox[2]), (2538644)))

nitrogen_200_utm_2_bbox <- nitrogen_200_utm_2 %>%
  st_bbox()
nitrogen_200_utm_2_crop <- crop(nitrogen_200_utm_2, extent((nitrogen_200_utm_2_bbox[1]), (nitrogen_200_utm_2_bbox[3]), (nitrogen_200_utm_2_bbox[2]), (2538644)))

Soil_Organic_Carbon_05_utm_2_bbox <- Soil_Organic_Carbon_05_utm_2 %>%
  st_bbox()
Soil_Organic_Carbon_05_utm_2_crop <- crop(Soil_Organic_Carbon_05_utm_2, extent((Soil_Organic_Carbon_05_utm_2_bbox[1]), (Soil_Organic_Carbon_05_utm_2_bbox[3]), (Soil_Organic_Carbon_05_utm_2_bbox[2]), (2538644)))

Soil_Organic_Carbon_200_utm_2_bbox <- Soil_Organic_Carbon_200_utm_2 %>%
  st_bbox()
Soil_Organic_Carbon_200_utm_2_crop <- crop(Soil_Organic_Carbon_200_utm_2, extent((Soil_Organic_Carbon_200_utm_2_bbox[1]), (Soil_Organic_Carbon_200_utm_2_bbox[3]), (Soil_Organic_Carbon_200_utm_2_bbox[2]), (2538644)))


#crop the soil raster for the southern portion of baja

#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_all_pop <- crop(clay_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4])) 
clay_200_all_pop <- crop(clay_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
silt_05_all_pop <- crop(silt_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
silt_200_all_pop <- crop(silt_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
sand_05_all_pop <- crop(sand_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
sand_200_all_pop <- crop(sand_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ph_05_all_pop <- crop(ph_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ph_200_all_pop <- crop(ph_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ocd_05_all_pop <- crop(ocd_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ocd_200_all_pop <- crop(ocd_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
coarse_frag_05_all_pop <- crop(coarse_frag_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
coarse_frag_200_all_pop <- crop(coarse_frag_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
cat_ex_cap_05_all_pop <- crop(cat_ex_cap_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
cat_ex_cap_200_all_pop <- crop(cat_ex_cap_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
bulk_dens_05_all_pop <- crop(bulk_dens_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
bulk_dens_200_all_pop <- crop(bulk_dens_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_10kpa_05_all_pop <- crop(vol_wat_10kpa_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_10kpa_200_all_pop <- crop(vol_wat_10kpa_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_33kpa_05_all_pop <- crop(vol_wat_33kpa_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_33kpa_200_all_pop <- crop(vol_wat_33kpa_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_1500kpa_05_all_pop <- crop(vol_wat_1500kpa_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_1500kpa_200_all_pop <- crop(vol_wat_1500kpa_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
nitrogen_05_all_pop <-  crop(nitrogen_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
nitrogen_200_all_pop <- crop(nitrogen_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
Soil_Organic_Carbon_05_all_pop <- crop(Soil_Organic_Carbon_05_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
Soil_Organic_Carbon_200_all_pop <- crop(Soil_Organic_Carbon_200_utm, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))



#crop the soil raster for the southern portion of baja, for the southern most area

#using the extent of the box around the rivers to crop the raster for each soil texture layer
clay_05_all_pop_2 <- crop(clay_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4])) 
clay_200_all_pop_2 <- crop(clay_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
silt_05_all_pop_2 <- crop(silt_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
silt_200_all_pop_2 <- crop(silt_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
sand_05_all_pop_2 <- crop(sand_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
sand_200_all_pop_2 <- crop(sand_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ph_05_all_pop_2 <- crop(ph_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ph_200_all_pop_2 <- crop(ph_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ocd_05_all_pop_2 <- crop(ocd_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
ocd_200_all_pop_2 <- crop(ocd_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
coarse_frag_05_all_pop_2 <- crop(coarse_frag_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
coarse_frag_200_all_pop_2 <- crop(coarse_frag_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
cat_ex_cap_05_all_pop_2 <- crop(cat_ex_cap_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
cat_ex_cap_200_all_pop_2 <- crop(cat_ex_cap_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
bulk_dens_05_all_pop_2 <- crop(bulk_dens_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
bulk_dens_200_all_pop_2 <- crop(bulk_dens_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_10kpa_05_all_pop_2 <- crop(vol_wat_10kpa_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_10kpa_200_all_pop_2 <- crop(vol_wat_10kpa_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_1500kpa_05_all_pop_2 <- crop(vol_wat_1500kpa_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
vol_wat_1500kpa_200_all_pop_2 <- crop(vol_wat_1500kpa_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
nitrogen_05_all_pop_2 <-  crop(nitrogen_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
nitrogen_200_all_pop_2 <- crop(nitrogen_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
Soil_Organic_Carbon_05_all_pop_2 <- crop(Soil_Organic_Carbon_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
Soil_Organic_Carbon_200_all_pop_2 <- crop(Soil_Organic_Carbon_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))



#attempt of using ggplot to plot clay layer with river shapefile
ggplot()+
  geom_raster(data = as.data.frame(clay_05_utm, xy=T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_raster(data = as.data.frame(clay_05_all_pop, xy=T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = BCS_polygon_UTM_cropped)

#attempt of using ggplot to plot clay layer with river shapefile
ggplot()+
  geom_raster(data = as.data.frame(clay_05_utm_2, xy=T), aes(x=x, y=y, fill = clay.content.0.5.part.2))+
  geom_sf(data = BCS_polygon_UTM_cropped)+
  geom_raster(data = as.data.frame(clay_05_all_pop_2, xy=T), aes(x=x, y=y, fill = clay.content.0.5.part.2))+
  geom_sf(data = river_LM_trans)
  

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
soil_stack_vol_wat_1500kpa <- stack(vol_wat_1500kpa_05_all_pop, vol_wat_1500kpa_200_all_pop) #stacked volume water content at 1500 kpa
soil_stack_nitrogen <- stack(nitrogen_05_all_pop, nitrogen_200_all_pop)  #stacked volume water content at 10 kpa
soil_stack_soc <- stack(Soil_Organic_Carbon_05_all_pop, Soil_Organic_Carbon_200_all_pop)  #stacked volume water content at 10 kpa



#creating a stack of the raster layers for the southern raster

soil_stack_clay_2 <- stack(clay_05_all_pop_2, clay_200_all_pop_2)  #stacked clay
soil_stack_silt_2 <- stack(silt_05_all_pop_2, silt_200_all_pop_2) #stacked silt
soil_stack_sand_2 <- stack(sand_05_all_pop_2, sand_200_all_pop_2) #stacked sand
soil_stack_ph_2 <- stack(ph_05_all_pop_2, ph_200_all_pop_2) #stacked ph
soil_stack_ocd_2 <- stack(ocd_05_all_pop_2, ocd_200_all_pop_2) #stacked ocd
soil_stack_coarse_frag_2 <- stack(coarse_frag_05_all_pop_2, coarse_frag_200_all_pop_2)#stacked coarse fragment
soil_stack_cat_ex_2 <- stack(cat_ex_cap_05_all_pop_2, cat_ex_cap_200_all_pop_2)  #stacked cation exchange capacity
soil_stack_bulk_dens_2 <- stack(bulk_dens_05_all_pop_2, bulk_dens_200_all_pop_2) #stacked bulk density
soil_stack_vol_wat_10kpa_2 <- stack(vol_wat_10kpa_05_all_pop_2, vol_wat_10kpa_200_all_pop_2) #stacked volume water content at 10 kpa
soil_stack_vol_wat_1500kpa_2 <- stack(Soil_Organic_Carbon_05_all_pop_2, vol_wat_1500kpa_200_all_pop_2) #stacked volume water content at 1500 kpa
soil_stack_nitrogen_2 <- stack(nitrogen_05_all_pop_2, nitrogen_200_all_pop_2)  #stacked volume water content at 10 kpa
soil_stack_soc_2 <- stack(Soil_Organic_Carbon_05_all_pop_2, Soil_Organic_Carbon_200_all_pop_2)  #stacked volume water content at 10 kpa


soil_stack_LM.df <- as.data.frame(getValues(soil_stack_LM))

#plotting the stacked rasters, example with clay
plot(soil_stack_clay) #version with soil textures
plot(soil_stack_clay, zlim = c(0, 400)) #version where the plots have the same scale


#extracting the soil data for each point of the known 16 points 

all_known_pop_soil_clay <- extract(soil_stack_clay, all_pop_locations.df_sf_trans_coordinates) #extracting soil textures for each point value
all_known_pop_soil_silt <- extract(soil_stack_silt, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_sand <- extract(soil_stack_sand, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_ph <- extract(soil_stack_ph, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_ocd <- extract(soil_stack_ocd, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_coarse_frag <- extract(soil_stack_coarse_frag, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_cat_ex <- extract(soil_stack_cat_ex, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_bulk_dens <- extract(soil_stack_bulk_dens, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_vol_wat_10kpa <- extract(soil_stack_vol_wat_10kpa, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_vol_wat_1500kpa <- extract(soil_stack_vol_wat_1500kpa, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_nitrogen <- extract(soil_stack_nitrogen, all_pop_locations.df_sf_trans_coordinates)
all_known_pop_soil_soc <- extract(soil_stack_soc, all_pop_locations.df_sf_trans_coordinates)

all_known_pop_soils <- cbind(all_pop_locations.df_sf_trans_coordinates, all_known_pop_soil_clay) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_silt) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_sand) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_ph) #bind the soil textures data for each point to the all point point dataframe
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_ocd)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_cat_ex)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_bulk_dens)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_vol_wat_10kpa)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_vol_wat_1500kpa)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_nitrogen)
all_known_pop_soils <- cbind(all_known_pop_soils, all_known_pop_soil_soc)

View(all_known_pop_soils)

#clay

#extracting means from randomly selected 16 points 

random_clay_0.5_means <- c() #creating empty list to collect means
random_clay_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations

  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_clay <- extract(soil_stack_clay, random_16)
  random_16_pop_soil_clay_2 <- extract(soil_stack_clay_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_clay)){ #adding in the missing clay values from the second raster
    if (is.na(random_16_pop_soil_clay[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_clay[i] = random_16_pop_soil_clay_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_clay_0.5_mean <- mean(random_16_pop_soil_clay[,1]) #storing the mean of the 0-5 value
  random_clay_100.200_mean <- mean(random_16_pop_soil_clay[,2]) #storing the mean of the 100-200 value
  
  random_clay_0.5_means <- c(random_clay_0.5_means, random_clay_0.5_mean) #adding the 0-5 mean to the list of means
  random_clay_100.200_means <- c(random_clay_100.200_means, random_clay_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_clay_0.5_mean <- mean(all_known_pop_soils$clay.content.0.5)
all_known_clay_100.200_mean <- mean(all_known_pop_soils$clay.content.100.200)

#for Clay 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean Clay 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_clay_0.5_means <- na.omit(random_clay_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_clay_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_0.5_means[i] > all_known_clay_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_clay_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for Clay 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_clay_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_clay_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean Clay 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_clay_100.200_means <- na.omit(random_clay_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_clay_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_clay_100.200_means[i] > all_known_clay_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_clay_100.200_means)) #the proportion of random ANNs that are less than our ANN




#silt

#extracting means from randomly selected 16 points 

random_silt_0.5_means <- c() #creating empty list to collect means
random_silt_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_silt <- extract(soil_stack_silt, random_16)
  random_16_pop_soil_silt_2 <- extract(soil_stack_silt_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_silt)){ #adding in the missing silt values from the second raster
    if (is.na(random_16_pop_soil_silt[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_silt[i] = random_16_pop_soil_silt_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_silt_0.5_mean <- mean(random_16_pop_soil_silt[,1]) #storing the mean of the 0-5 value
  random_silt_100.200_mean <- mean(random_16_pop_soil_silt[,2]) #storing the mean of the 100-200 value
  
  random_silt_0.5_means <- c(random_silt_0.5_means, random_silt_0.5_mean) #adding the 0-5 mean to the list of means
  random_silt_100.200_means <- c(random_silt_100.200_means, random_silt_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_silt_0.5_mean <- mean(all_known_pop_soils$silt.0.5)
all_known_silt_100.200_mean <- mean(all_known_pop_soils$silt.100.200)

#for silt 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_silt_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_silt_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean silt 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_silt_0.5_means <- na.omit(random_silt_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_silt_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_silt_0.5_means[i] < all_known_silt_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_silt_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for silt 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_silt_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_silt_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean silt 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_silt_100.200_means <- na.omit(random_silt_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_silt_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_silt_100.200_means[i] > all_known_silt_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_silt_100.200_means)) #the proportion of random ANNs that are less than our ANN



#sand

#extracting means from randomly selected 16 points 

random_sand_0.5_means <- c() #creating empty list to collect means
random_sand_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_sand <- extract(soil_stack_sand, random_16)
  random_16_pop_soil_sand_2 <- extract(soil_stack_sand_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_sand)){ #adding in the missing sand values from the second raster
    if (is.na(random_16_pop_soil_sand[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_sand[i] = random_16_pop_soil_sand_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_sand_0.5_mean <- mean(random_16_pop_soil_sand[,1]) #storing the mean of the 0-5 value
  random_sand_100.200_mean <- mean(random_16_pop_soil_sand[,2]) #storing the mean of the 100-200 value
  
  random_sand_0.5_means <- c(random_sand_0.5_means, random_sand_0.5_mean) #adding the 0-5 mean to the list of means
  random_sand_100.200_means <- c(random_sand_100.200_means, random_sand_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

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
(total / length(random_sand_0.5_means)) #the proportion of random ANNs that are less than our ANN

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
  if (random_sand_100.200_means[i] < all_known_sand_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_sand_100.200_means)) #the proportion of random ANNs that are less than our ANN


#ph

#extracting means from randomly selected 16 points 

random_ph_0.5_means <- c() #creating empty list to collect means
random_ph_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_ph <- extract(soil_stack_ph, random_16)
  random_16_pop_soil_ph_2 <- extract(soil_stack_ph_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_ph)){ #adding in the missing ph values from the second raster
    if (is.na(random_16_pop_soil_ph[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_ph[i] = random_16_pop_soil_ph_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_ph_0.5_mean <- mean(random_16_pop_soil_ph[,1]) #storing the mean of the 0-5 value
  random_ph_100.200_mean <- mean(random_16_pop_soil_ph[,2]) #storing the mean of the 100-200 value
  
  random_ph_0.5_means <- c(random_ph_0.5_means, random_ph_0.5_mean) #adding the 0-5 mean to the list of means
  random_ph_100.200_means <- c(random_ph_100.200_means, random_ph_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

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
total = 0  #set empty vaue
for (i in 1:length(random_ph_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_ph_0.5_means[i] < all_known_ph_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_ph_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for ph 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_ph_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_ph_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean ph 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_ph_100.200_means <- na.omit(random_ph_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_ph_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_ph_100.200_means[i] < all_known_ph_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_ph_100.200_means)) #the proportion of random ANNs that are less than our ANN

Soil_Organic_Carbon_05_all_pop_2 <- crop(Soil_Organic_Carbon_05_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))
Soil_Organic_Carbon_200_all_pop_2 <- crop(Soil_Organic_Carbon_200_utm_2_crop, extent(BCS_polygon_crop_bbox[1], BCS_polygon_crop_bbox[3], BCS_polygon_crop_bbox[2], BCS_polygon_crop_bbox[4]))



#soc

#extracting means from randomly selected 16 points 

random_soc_0.5_means <- c() #creating empty list to collect means
random_soc_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_soc <- extract(soil_stack_soc, random_16)
  random_16_pop_soil_soc_2 <- extract(soil_stack_soc_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_soc)){ #adding in the missing soc values from the second raster
    if (is.na(random_16_pop_soil_soc[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_soc[i] = random_16_pop_soil_soc_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_soc_0.5_mean <- mean(random_16_pop_soil_soc[,1]) #storing the mean of the 0-5 value
  random_soc_100.200_mean <- mean(random_16_pop_soil_soc[,2]) #storing the mean of the 100-200 value
  
  random_soc_0.5_means <- c(random_soc_0.5_means, random_soc_0.5_mean) #adding the 0-5 mean to the list of means
  random_soc_100.200_means <- c(random_soc_100.200_means, random_soc_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_soc_0.5_mean <- mean(all_known_pop_soils$SOC.0.5)
all_known_soc_100.200_mean <- mean(all_known_pop_soils$SOC.100.200)

#for soc 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_soc_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_soc_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean soc 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_soc_0.5_means <- na.omit(random_soc_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_soc_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_soc_0.5_means[i] > all_known_soc_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_soc_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for soc 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_soc_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_soc_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean soc 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_soc_100.200_means <- na.omit(random_soc_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_soc_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_soc_100.200_means[i] > all_known_soc_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_soc_100.200_means)) #the proportion of random ANNs that are less than our ANN



#vol_wat_10kpa

#extracting means from randomly selected 16 points 

random_vol_wat_10kpa_0.5_means <- c() #creating empty list to collect means
random_vol_wat_10kpa_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_vol_wat_10kpa <- extract(soil_stack_vol_wat_10kpa, random_16)
  random_16_pop_soil_vol_wat_10kpa_2 <- extract(soil_stack_vol_wat_10kpa_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_vol_wat_10kpa)){ #adding in the missing vol_wat_10kpa values from the second raster
    if (is.na(random_16_pop_soil_vol_wat_10kpa[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_vol_wat_10kpa[i] = random_16_pop_soil_vol_wat_10kpa_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_vol_wat_10kpa_0.5_mean <- mean(random_16_pop_soil_vol_wat_10kpa[,1]) #storing the mean of the 0-5 value
  random_vol_wat_10kpa_100.200_mean <- mean(random_16_pop_soil_vol_wat_10kpa[,2]) #storing the mean of the 100-200 value
  
  random_vol_wat_10kpa_0.5_means <- c(random_vol_wat_10kpa_0.5_means, random_vol_wat_10kpa_0.5_mean) #adding the 0-5 mean to the list of means
  random_vol_wat_10kpa_100.200_means <- c(random_vol_wat_10kpa_100.200_means, random_vol_wat_10kpa_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_vol_wat_10kpa_0.5_mean <- mean(all_known_pop_soils$vol_water_.10_0.5)
all_known_vol_wat_10kpa_100.200_mean <- mean(all_known_pop_soils$vol_water_.10_100.200)

#for vol_wat_10kpa 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_10kpa_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_10kpa_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_10kpa 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_vol_wat_10kpa_0.5_means <- na.omit(random_vol_wat_10kpa_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_10kpa_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_10kpa_0.5_means[i] > all_known_vol_wat_10kpa_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_vol_wat_10kpa_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for vol_wat_10kpa 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_10kpa_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_10kpa_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_10kpa 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_vol_wat_10kpa_100.200_means <- na.omit(random_vol_wat_10kpa_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_10kpa_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_10kpa_100.200_means[i] > all_known_vol_wat_10kpa_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_vol_wat_10kpa_100.200_means)) #the proportion of random ANNs that are less than our ANN


#vol_wat_1500kpa

#extracting means from randomly selected 16 points 

random_vol_wat_1500kpa_0.5_means <- c() #creating empty list to collect means
random_vol_wat_1500kpa_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_vol_wat_1500kpa <- extract(soil_stack_vol_wat_1500kpa, random_16)
  random_16_pop_soil_vol_wat_1500kpa_2 <- extract(soil_stack_vol_wat_1500kpa_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_vol_wat_1500kpa)){ #adding in the missing vol_wat_1500kpa values from the second raster
    if (is.na(random_16_pop_soil_vol_wat_1500kpa[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_vol_wat_1500kpa[i] = random_16_pop_soil_vol_wat_1500kpa_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_vol_wat_1500kpa_0.5_mean <- mean(random_16_pop_soil_vol_wat_1500kpa[,1]) #storing the mean of the 0-5 value
  random_vol_wat_1500kpa_100.200_mean <- mean(random_16_pop_soil_vol_wat_1500kpa[,2]) #storing the mean of the 100-200 value
  
  random_vol_wat_1500kpa_0.5_means <- c(random_vol_wat_1500kpa_0.5_means, random_vol_wat_1500kpa_0.5_mean) #adding the 0-5 mean to the list of means
  random_vol_wat_1500kpa_100.200_means <- c(random_vol_wat_1500kpa_100.200_means, random_vol_wat_1500kpa_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_vol_wat_1500kpa_0.5_mean <- mean(all_known_pop_soils$vol_water_.10_0.5)
all_known_vol_wat_1500kpa_100.200_mean <- mean(all_known_pop_soils$vol_water_.10_100.200)

#for vol_wat_1500kpa 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_1500kpa_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_1500kpa_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_1500kpa 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_vol_wat_1500kpa_0.5_means <- na.omit(random_vol_wat_1500kpa_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_1500kpa_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_1500kpa_0.5_means[i] > all_known_vol_wat_1500kpa_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_vol_wat_1500kpa_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for vol_wat_1500kpa 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_vol_wat_1500kpa_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_vol_wat_1500kpa_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean vol_wat_1500kpa 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_vol_wat_1500kpa_100.200_means <- na.omit(random_vol_wat_1500kpa_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_vol_wat_1500kpa_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_vol_wat_1500kpa_100.200_means[i] > all_known_vol_wat_1500kpa_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_vol_wat_1500kpa_100.200_means)) #the proportion of random ANNs that are less than our ANN


#nitrogen

#extracting means from randomly selected 16 points 

random_nitrogen_0.5_means <- c() #creating empty list to collect means
random_nitrogen_100.200_means <- c() #creating empty list to collect means

set.seed(20)
for (i in 1:1000){ #for 1000 permutations
  
  random_16 <- st_sample(BCS_polygon_UTM_cropped, 16) #select rando 16 points within the cropped BCS polygon
  random_16 <- random_16 %>%
    st_as_sf()
  random_16_pop_soil_nitrogen <- extract(soil_stack_nitrogen, random_16)
  random_16_pop_soil_nitrogen_2 <- extract(soil_stack_nitrogen_2, random_16)
  
  for (i in 1:length(random_16_pop_soil_nitrogen)){ #adding in the missing nitrogen values from the second raster
    if (is.na(random_16_pop_soil_nitrogen[i]) == T) { #checking if there are NAs in the list of values
      random_16_pop_soil_nitrogen[i] = random_16_pop_soil_nitrogen_2[i] #if there are NAs replace iwth value from the other raster
    }
  }
  
  random_nitrogen_0.5_mean <- mean(random_16_pop_soil_nitrogen[,1]) #storing the mean of the 0-5 value
  random_nitrogen_100.200_mean <- mean(random_16_pop_soil_nitrogen[,2]) #storing the mean of the 100-200 value
  
  random_nitrogen_0.5_means <- c(random_nitrogen_0.5_means, random_nitrogen_0.5_mean) #adding the 0-5 mean to the list of means
  random_nitrogen_100.200_means <- c(random_nitrogen_100.200_means, random_nitrogen_100.200_mean) #adding the 100-200 mean to the list of means
  
}

#plotting the randomly selected points
ggplot()+
  geom_sf(data=BCS_polygon_box_sf)+
  geom_sf(data=BCS_polygon_box_sf_cropped)+
  geom_sf(data=BCS_polygon_UTM)+
  geom_sf(data=BCS_polygon_UTM_cropped, color = "red")+
  geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
  geom_sf(data=random_16, color ="blue")

#storing the real means
all_known_nitrogen_0.5_mean <- mean(all_known_pop_soils$nitrogen.0.5)
all_known_nitrogen_100.200_mean <- mean(all_known_pop_soils$nitrogen.100.200)

#for nitrogen 0-5

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_nitrogen_0.5_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_nitrogen_0.5_mean, col = "red")+ #line of our real slope
  xlab("Mean nitrogen 0-5 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_nitrogen_0.5_means <- na.omit(random_nitrogen_0.5_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_nitrogen_0.5_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_nitrogen_0.5_means[i] > all_known_nitrogen_0.5_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_nitrogen_0.5_means)) #the proportion of random ANNs that are less than our ANN

#for nitrogen 100-200

#plotting the histogram of the randomly distributed p-values and our real slope
ggplot()+
  geom_histogram(aes(x=random_nitrogen_100.200_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
  geom_vline(xintercept=all_known_nitrogen_100.200_mean, col = "red")+ #line of our real slope
  xlab("Mean nitrogen 100-200 of Random Populations vs.Known Populations (n=16)")+
  theme_classic()

random_nitrogen_100.200_means <- na.omit(random_nitrogen_100.200_means) #remove NAs

#calculating pseudo p-value for 
total = 0  #set empty vaue
for (i in 1:length(random_nitrogen_100.200_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
  if (random_nitrogen_100.200_means[i] > all_known_nitrogen_100.200_mean){
    total = total + 1
  }
} #add number of values of in the random set of ANN values that are less than our mean ANN
(total / length(random_nitrogen_100.200_means)) #the proportion of random ANNs that are less than our ANN



# creating a heat map of all of the p-values

LM_soil_values_y <- 
LM_size_variables_x <-
LM_p_values <- 
data <- expand.grid(X=x, Y= LM_soil_values_y)



