# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Quercus brandegeei tree shape are Influenced by Their Soil Metrics%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to determine whether there is a significant relationship between the 
# size (short canopy axis, DBH, etc.) of Q. brandegeei trees at each population and their soil
# characteristics (clay content, pH, etc.).
# A significant relationship could indicate which soil characteristics may influence the size and shape of the 
# trees and potentially help explain their distribution. 

# To test this, we used Simple Linear Regressions for each population, using the soil metrics as the explanatory 
# variables and each size metric as the response variable. We used a permutation test, comparing the slope of the 
# real linear regressions to a distribution of slopes for the linear regression if the size and shape metrics 
# were randomly shuffled, to support if the relationship is significant without assumption a normal distribution.

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
          # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
          # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
#boxs around the  rivers, stacking the rasters for each population, and processing them into one dataframe for all and each population
          # Creating the four new soil metrics: Sandy available Water (0-5 and 100-200 cm) and Clay/Loam Available Water (0-5 and 100-200 cm)
# 2) Making a function for calculating the real and randomized single linear regressions/slopes comparing tree 
# size to soil metrics for each population and calculating the p-values. 
# 3) Running the function and storing the output, graphing the histograms/slopes, turning the results into a dataframe
# 4) Using tables and heat maps to summarize the results (the slopes and their significances)

#### Loading libraries and Relevant Data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) #ggplot extention
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

# loading in the processed tree data 
source("./analyses/Data_Processing_Script.R")

# # loading in the tree data (size, elevation, lat/lon, ID, size/shape)
# 
# fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R
# 
# #adding a sequential column, "X," to number each tree
# 
# fixed_field_data_processed <- fixed_field_data_processed %>%
#   mutate(X = row_number())
# 
# #creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
# #sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
# fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
#                                           coords = c("long", "lat"), crs = 4326)
# 
# #creating a transformed point shapefile with UTM 12 N an equal area projection
# fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) 
# 
# #create dataframe with X and Y UTM coordinates
# 
# fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with seperate x and y columns from the UTM 12N transformation
# fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
#   cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe
# 
# #add average nearest neighbor for each individual column
# fixed_field_data_processed_NN_UTM <- fixed_field_data_processed_sf_trans_coordinates %>%  #creates a dataframe with the ANN of the closest 5 individual trees for each individual
#   mutate(dist1 = nndist(X = X.1, Y= Y, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
#   mutate(dist2 = nndist(X = X.1, Y= Y, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
#   mutate(dist3 = nndist(X = X.1, Y= Y, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
#   mutate(dist4 = nndist(X = X.1, Y= Y, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
#   mutate(dist5 = nndist(X = X.1, Y= Y, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
#   rowwise()%>% #so that in the next part we take the averages across rows
#   mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5)))  %>% #creates a column of the average distances (1-5) of each individual
#   dplyr::select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances
# 
# soil_stack_clay <- stack(clay_05_all_pop, clay_200_all_pop) #stacked clay
# soil_stack_silt <- stack(silt_05_all_pop, silt_200_all_pop) #stacked silt
# soil_stack_sand <- stack(sand_05_all_pop, sand_200_all_pop) #stacked sand
# soil_stack_ph <- stack(ph_05_all_pop, ph_200_all_pop) #stacked ph
# soil_stack_ocd <- stack(ocd_05_all_pop, ocd_200_all_pop)  #stacked ocd
# soil_stack_coarse_frag <- stack(coarse_frag_05_all_pop, coarse_frag_200_all_pop)  #stacked coarse fragment
# soil_stack_cat_ex <- stack(cat_ex_cap_05_all_pop, cat_ex_cap_200_all_pop)  #stacked cation exchange capacity
# soil_stack_bulk_dens <- stack(bulk_dens_05_all_pop, bulk_dens_200_all_pop)  #stacked bulk density
# soil_stack_vol_wat_10kpa <- stack(vol_wat_10kpa_05_all_pop, vol_wat_10kpa_200_all_pop)  #stacked volume water content at -10 kpa
# soil_stack_vol_wat_33kpa <- stack(vol_wat_33kpa_05_all_pop, vol_wat_33kpa_200_all_pop)  #stacked volume water content at -33 kpa
# soil_stack_vol_wat_1500kpa <- stack(vol_wat_1500kpa_05_all_pop, vol_wat_1500kpa_200_all_pop) #stacked volume water content at -1500 kpa
# soil_stack_nitrogen <- stack(nitrogen_05_all_pop, nitrogen_200_all_pop)  #stacked nitrogen
# soil_stack_soc <- stack(Soil_Organic_Carbon_05_all_pop, Soil_Organic_Carbon_200_all_pop)  #stacked soil organic carbon
# soil_stack_sandy_water <- stack(sandy_avail_water_0.5_all_pop, sandy_avail_water_100.200_all_pop) #stacked sandy available water
# soil_stack_clay_loam_water <- stack(clay_loam_avail_water_0.5_all_pop, clay_loam_avail_water_100.200_all_pop) #stacked clay/loam available water
# 
# 
# #creating shapefiles for each population, turning sf of all points into sfc
# 
# LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LM") %>%
#   st_as_sfc()
# 
# LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LC") %>%
#   st_as_sfc()
# 
# SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "SD") %>%
#   st_as_sfc()
# 
# # Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns
# 
# LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
#   filter(Locality == "LM")
# 
# LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
#   filter(Locality == "LC")
# 
# SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
#   filter(Locality == "SD")
# 
# #Upload ArcGIS river shapefile and filter out polygons for each population
# 
# river_LM <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LM River/LM_Rivers_Final.shp")
# river_LM  <- river_LM$geometry[1]
# plot(river_LM)
# 
# river_LC  <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LC River/LC_Rivers_Final.shp")
# river_LC  <- river_LC$geometry[1]
# plot(river_LC)
# 
# river_SD <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/SD River/SD_Rivers_Final.shp")
# river_SD <- river_SD$geometry[1]
# plot(river_SD)
# 
# #changing the coordinate reference system of the river polygons to be equal area projection (UTM 12N), uses meters as distance measurement 
# river_LM_trans <- st_as_sf(st_transform(river_LM, crs = 26912))
# river_LC_trans <- st_as_sf(st_transform(river_LC, crs = 26912))
# river_SD_trans <- st_as_sf(st_transform(river_SD, crs = 26912))
# 
# #creating bounding boxes for each population
# 
# #creating a boundry box of LM with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
# LM_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LM") %>%
#   st_bbox %>%
#   st_as_sfc()
# 
# #creating a boundry box of LC with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
# LC_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LC") %>%
#   st_bbox %>%
#   st_as_sfc()
# 
# #creating a boundry box of SD with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
# SD_fixed_field_data_processed_box <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "SD") %>%
#   st_bbox %>%
#   st_as_sfc()
# 
# #creating bboxs for all of the river shapefiles for each population
# LM_box <- st_bbox(river_LM_trans)
# LC_box <- st_bbox(river_LC_trans)
# SD_box <- st_bbox(river_SD_trans)
# 
# ## Load in Soil Metric Data ##
# 
# # loading in soil textures raster tifs from CONABIO
# 
# clay_05 <- raster(paste0("./data/Soil Grid/clay content/clay content 0-5.tif"))
# clay_200 <- raster(paste0("./data/Soil Grid/clay content/clay content 100-200.tif"))
# silt_05 <- raster(paste0("./data/Soil Grid/silt/silt 0-5.tif"))
# silt_200 <-raster(paste0("./data/Soil Grid/silt/silt 100-200.tif"))
# sand_05 <- raster(paste0("./data/Soil Grid/sand/sand 0-5.tif"))
# sand_200 <- raster(paste0("./data/Soil Grid/sand/sand 100-200.tif"))
# ph_05 <- raster(paste0("./data/Soil Grid/pH/ph_0-5.tif")) #0-5 cm ph
# ph_200 <- raster(paste0("./data/Soil Grid/pH/ph_100-200.tif")) #100-200 ph
# ocd_05 <- raster(paste0("./data/Soil Grid/organic carbon density/OCD_0-5.tif")) #0-5cm organic carbon density
# ocd_200 <- raster(paste0("./data/Soil Grid/organic carbon density/OCR_100-200.tif")) #100-200cm organic carbon density
# coarse_frag_05 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_0-5.tif")) #0-5 cm coarse fragments
# coarse_frag_200 <- raster(paste0("./data/Soil Grid/coarse fragments/coarse_fragments_100-200.tif")) #100-200 cm coarse fragments
# cat_ex_cap_05 <-raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_0-5.tif")) #0-5 cm cation exchange capacity
# cat_ex_cap_200 <- raster(paste0("./data/Soil Grid/cation exchange capacity/Cat_exc_100-200.tif")) #100-200 cm cation exchange capacity
# bulk_dens_05 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_0-5.tif")) #0-5 cm bulk density
# bulk_dens_200 <- raster(paste0("./data/Soil Grid/bulk density/bulk_density_100-200.tif")) #100-200 cm bulk density
# vol_wat_10kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_0-5.tif"))  #0-5 cm -10 kpa volumn water content
# vol_wat_10kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -10 kPa/vol_water_-10_100-200.tif"))  #100-200 cm -10 kpa volumn water content
# vol_wat_33kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -33 kPa /vol_water_0-5.tif")) #0-5 cm -33 kpa volumn water content
# vol_wat_33kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -33 kPa /vol_water_100-200.tif")) #100-200 cm -33 kpa volumn water content
# vol_wat_1500kpa_05 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500kPa_0-5.tif"))  #0-5 cm -1500 kpa volumn water content
# vol_wat_1500kpa_200 <- raster(paste0("./data/Soil Grid/vol. water content at -1500 kPa/vol_water_-1500_100-200.tif")) #100-200 cm -1500 kpa volumn water content
# nitrogen_05 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 0-5.tif")) #0-5 cm nitrogen content 
# nitrogen_200 <- raster(paste0("./data/Soil Grid/Nitrogen/nitrogen 100-200.tif")) #100-200 cm nitrogen content 
# Soil_Organic_Carbon_05 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 0-5.tif")) #0-5 cm soil organic carbon
# Soil_Organic_Carbon_200 <- raster(paste0("./data/Soil Grid/Soil Organic Carbon/SOC 100-200.tif")) #100-200 cm soil organic carbon
# 
# #project the soil metric rasters to equal area projection (UTM 12N) which uses meters as distance measurement 
# 
# clay_05_utm <- projectRaster(clay_05, crs=26912) #converting the 0-5 cm clay raster to utm 12
# clay_200_utm <- projectRaster(clay_200, crs=26912) #converting the 90-200 cm clay raster to utm 12
# silt_05_utm <- projectRaster(silt_05, crs=26912)
# silt_200_utm <- projectRaster(silt_200, crs=26912)
# sand_05_utm <- projectRaster(sand_05, crs=26912)
# sand_200_utm <- projectRaster(sand_200, crs=26912)
# ph_05_utm <- projectRaster(ph_05, crs=26912) 
# ph_200_utm <- projectRaster(ph_200, crs=26912) 
# ocd_05_utm <- projectRaster(ocd_05, crs=26912)
# ocd_200_utm <- projectRaster(ocd_200, crs=26912)
# coarse_frag_05_utm <- projectRaster(coarse_frag_05, crs=26912)
# coarse_frag_200_utm <- projectRaster(coarse_frag_200, crs=26912)
# cat_ex_cap_05_utm <- projectRaster(cat_ex_cap_05, crs=26912)
# cat_ex_cap_200_utm <- projectRaster(cat_ex_cap_200, crs=26912)
# bulk_dens_05_utm <- projectRaster(bulk_dens_05, crs=26912)
# bulk_dens_200_utm <- projectRaster(bulk_dens_200, crs=26912)
# vol_wat_10kpa_05_utm <- projectRaster(vol_wat_10kpa_05, crs=26912)
# vol_wat_10kpa_200_utm <- projectRaster(vol_wat_10kpa_200, crs=26912)
# vol_wat_33kpa_05_utm <- projectRaster(vol_wat_33kpa_05, crs=26912)
# vol_wat_33kpa_200_utm <- projectRaster(vol_wat_33kpa_200, crs=26912)
# vol_wat_1500kpa_05_utm <- projectRaster(vol_wat_1500kpa_05, crs=26912)
# vol_wat_1500kpa_200_utm <- projectRaster(vol_wat_1500kpa_200, crs=26912)
# nitrogen_05_utm <- projectRaster(nitrogen_05, crs=26912)
# nitrogen_200_utm <- projectRaster(nitrogen_200, crs=26912)
# Soil_Organic_Carbon_05_utm <- projectRaster(Soil_Organic_Carbon_05, crs=26912)
# Soil_Organic_Carbon_200_utm <- projectRaster(Soil_Organic_Carbon_200, crs=26912)
# 
# ## cropping the soil metrics to bounding boxes around each population ##
# 
# # LM
# 
# # using the extent of the box around the rivers to crop the raster for each soil texture layer
# clay_05_LM <- crop(clay_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100)) 
# clay_200_LM <- crop(clay_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# silt_05_LM <- crop(silt_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# silt_200_LM <- crop(silt_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# sand_05_LM <- crop(sand_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# sand_200_LM <- crop(sand_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# ph_05_LM <- crop(ph_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# ph_200_LM <- crop(ph_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# ocd_05_LM <- crop(ocd_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# ocd_200_LM <- crop(ocd_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# coarse_frag_05_LM <- crop(coarse_frag_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# coarse_frag_200_LM <- crop(coarse_frag_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# cat_ex_cap_05_LM <- crop(cat_ex_cap_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# cat_ex_cap_200_LM <- crop(cat_ex_cap_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# bulk_dens_05_LM <- crop(bulk_dens_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# bulk_dens_200_LM <- crop(bulk_dens_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_10kpa_05_LM <- crop(vol_wat_10kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_10kpa_200_LM <- crop(vol_wat_10kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_33kpa_05_LM <- crop(vol_wat_33kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_33kpa_200_LM <- crop(vol_wat_33kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_1500kpa_05_LM <- crop(vol_wat_1500kpa_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# vol_wat_1500kpa_200_LM <- crop(vol_wat_1500kpa_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# nitrogen_05_LM <-  crop(nitrogen_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# nitrogen_200_LM <- crop(nitrogen_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# Soil_Organic_Carbon_05_LM <- crop(Soil_Organic_Carbon_05_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# Soil_Organic_Carbon_200_LM <- crop(Soil_Organic_Carbon_200_utm, extent(LM_box[1]-100, LM_box[3]+100, LM_box[2]-100, LM_box[4]+100))
# 
# 
# # plotting an example soil metric raster with river shapefile and tree points  
# ggplot()+
#   geom_raster(data = as.data.frame(Soil_Organic_Carbon_05_LM, xy=T), aes(x=x, y=y, fill = SOC.0.5))+
#   geom_sf(data = river_LM_trans)+
#   geom_sf(data = LM_fixed_field_data_processed)
# 
# #creating stacks of the soil raster layers (using multiple for visualization purposes)
# soil_stack_LM_soil_text <- stack(clay_05_LM, clay_200_LM, silt_05_LM, silt_200_LM, sand_05_LM, sand_200_LM) #the stack of all of the soil texture rasters
# soil_stack_LM_other <- stack(ph_05_LM, ph_200_LM, ocd_05_LM, ocd_200_LM, coarse_frag_05_LM, coarse_frag_200_LM, #the stack of all of the other soil variables, with different extents than the soil texture rasters
#                              cat_ex_cap_05_LM, cat_ex_cap_200_LM, bulk_dens_05_LM, bulk_dens_200_LM, vol_wat_10kpa_05_LM,
#                              vol_wat_10kpa_200_LM, vol_wat_33kpa_05_LM, vol_wat_33kpa_200_LM, vol_wat_1500kpa_05_LM, 
#                              vol_wat_1500kpa_200_LM) 
# soil_stack_LM_extra <- stack(nitrogen_05_LM, nitrogen_200_LM, Soil_Organic_Carbon_05_LM, Soil_Organic_Carbon_200_LM)
# 
# 
# 
# #plotting the stacked rasters
# plot(soil_stack_LM_soil_text) #version with soil textures
# plot(soil_stack_LM_soil_text, zlim = c(100, 710)) #version where the plots have the same scale
# plot(soil_stack_LM_other) #version with other variables
# plot(soil_stack_LM_other, zlim = c(30, 360)) #version where the plots have the same scale
# plot(soil_stack_LM_extra) #version with other variables
# plot(soil_stack_LM_extra, zlim = c(30, 360)) #version where the plots have the same scale
# 
# 
# #LC
# 
# # using the extent of the box around the rivers to crop the raster for each soil texture layer
# clay_05_LC <- crop(clay_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100)) 
# clay_200_LC <- crop(clay_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# silt_05_LC <- crop(silt_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# silt_200_LC <- crop(silt_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# sand_05_LC <- crop(sand_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# sand_200_LC <- crop(sand_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# ph_05_LC <- crop(ph_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# ph_200_LC <- crop(ph_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# ocd_05_LC <- crop(ocd_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# ocd_200_LC <- crop(ocd_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# coarse_frag_05_LC <- crop(coarse_frag_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# coarse_frag_200_LC <- crop(coarse_frag_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# cat_ex_cap_05_LC <- crop(cat_ex_cap_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# cat_ex_cap_200_LC <- crop(cat_ex_cap_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# bulk_dens_05_LC <- crop(bulk_dens_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# bulk_dens_200_LC <- crop(bulk_dens_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_10kpa_05_LC <- crop(vol_wat_10kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_10kpa_200_LC <- crop(vol_wat_10kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_33kpa_05_LC <- crop(vol_wat_33kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_33kpa_200_LC <- crop(vol_wat_33kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_1500kpa_05_LC <- crop(vol_wat_1500kpa_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# vol_wat_1500kpa_200_LC <- crop(vol_wat_1500kpa_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# nitrogen_05_LC <-  crop(nitrogen_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# nitrogen_200_LC <- crop(nitrogen_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# Soil_Organic_Carbon_05_LC <- crop(Soil_Organic_Carbon_05_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# Soil_Organic_Carbon_200_LC <- crop(Soil_Organic_Carbon_200_utm, extent(LC_box[1]-100, LC_box[3]+100, LC_box[2]-100, LC_box[4]+100))
# 
# #creating stacks of the soil raster layers (using multiple for visualization purposes)
# soil_stack_LC_soil_text <- stack(clay_05_LC, clay_200_LC, silt_05_LC, silt_200_LC, sand_05_LC, sand_200_LC) #the stack of all of the soil texture rasters
# soil_stack_LC_other <- stack(ph_05_LC, ph_200_LC, ocd_05_LC, ocd_200_LC, coarse_frag_05_LC, coarse_frag_200_LC, #the stack of all of the other soil variables, with different extents than the soil texture rasters
#                              cat_ex_cap_05_LC, cat_ex_cap_200_LC, bulk_dens_05_LC, bulk_dens_200_LC, vol_wat_10kpa_05_LC,
#                              vol_wat_10kpa_200_LC, vol_wat_33kpa_05_LC, vol_wat_33kpa_200_LC, vol_wat_1500kpa_05_LC, 
#                              vol_wat_1500kpa_200_LC) 
# soil_stack_LC_extra <- stack(nitrogen_05_LC, nitrogen_200_LC, Soil_Organic_Carbon_05_LC, Soil_Organic_Carbon_200_LC)
# 
# #plotting the stacked rasters
# plot(soil_stack_LC_soil_text) #version with soil textures
# plot(soil_stack_LC_soil_text, zlim = c(100, 710)) #version where the plots have the same scale
# plot(soil_stack_LC_other) #version with other variables
# plot(soil_stack_LC_other, zlim = c(30, 360)) #version where the plots have the same scale
# plot(soil_stack_LC_extra) #version with other variables
# plot(soil_stack_LC_extra, zlim = c(30, 180)) #version where the plots have the same scale
# 
# 
# #SD
# 
# #using the extent of the box around the rivers to crop the raster for each soil texture layer
# clay_05_SD <- crop(clay_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100)) 
# clay_200_SD <- crop(clay_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# silt_05_SD <- crop(silt_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# silt_200_SD <- crop(silt_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# sand_05_SD <- crop(sand_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# sand_200_SD <- crop(sand_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# ph_05_SD <- crop(ph_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# ph_200_SD <- crop(ph_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# ocd_05_SD <- crop(ocd_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# ocd_200_SD <- crop(ocd_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# coarse_frag_05_SD <- crop(coarse_frag_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# coarse_frag_200_SD <- crop(coarse_frag_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# cat_ex_cap_05_SD <- crop(cat_ex_cap_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# cat_ex_cap_200_SD <- crop(cat_ex_cap_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# bulk_dens_05_SD <- crop(bulk_dens_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# bulk_dens_200_SD <- crop(bulk_dens_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_10kpa_05_SD <- crop(vol_wat_10kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_10kpa_200_SD <- crop(vol_wat_10kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_33kpa_05_SD <- crop(vol_wat_33kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_33kpa_200_SD <- crop(vol_wat_33kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_1500kpa_05_SD <- crop(vol_wat_1500kpa_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# vol_wat_1500kpa_200_SD <- crop(vol_wat_1500kpa_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# nitrogen_05_SD <-  crop(nitrogen_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# nitrogen_200_SD <- crop(nitrogen_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# Soil_Organic_Carbon_05_SD <- crop(Soil_Organic_Carbon_05_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# Soil_Organic_Carbon_200_SD <- crop(Soil_Organic_Carbon_200_utm, extent(SD_box[1]-100, SD_box[3]+100, SD_box[2]-100, SD_box[4]+100))
# 
# #creating stacks of the soil raster layers (using multiple for visualization purposes)
# soil_stack_SD_soil_text <- stack(clay_05_SD, clay_200_SD, silt_05_SD, silt_200_SD, sand_05_SD, sand_200_SD) #the stack of all of the soil texture rasters
# soil_stack_SD_other <- stack(ph_05_SD, ph_200_SD, ocd_05_SD, ocd_200_SD, coarse_frag_05_SD, coarse_frag_200_SD, #the stack of all of the other soil variables, with different extents than the soil texture rasters
#                              cat_ex_cap_05_SD, cat_ex_cap_200_SD, bulk_dens_05_SD, bulk_dens_200_SD, vol_wat_10kpa_05_SD,
#                              vol_wat_10kpa_200_SD, vol_wat_33kpa_05_SD, vol_wat_33kpa_200_SD, vol_wat_1500kpa_05_SD, 
#                              vol_wat_1500kpa_200_SD) 
# soil_stack_SD_extra <- stack(nitrogen_05_SD, nitrogen_200_SD,Soil_Organic_Carbon_05_SD,  Soil_Organic_Carbon_200_SD)
# 
# #plotting the stacked rasters
# plot(soil_stack_SD_soil_text)
# plot(soil_stack_SD_soil_text, zlim = c(130, 710)) #version where the plots have the same scale
# plot(soil_stack_SD_other)
# plot(soil_stack_SD_other, zlim = c(45, 360)) #version where the plots have the same scale
# plot(soil_stack_SD_extra)
# plot(soil_stack_SD_extra, zlim = c(25, 340)) #version where the plots have the same scale
# 
# #creating X sequential columns in LC and SD point data which will make it easier to select random points from each grid later
# 
# #creating an x_sequential column that is 1 through the number of LM points, which will make it easier to randomly choose one point
# LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
#   mutate(X_sequential = 1:nrow(LM_fixed_field_data_processed))
# 
# #creating an x_sequential column that is 1 through the number of LC points, which will make it easier to randomly choose one point
# LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
#   mutate(X_sequential = 1:nrow(LC_fixed_field_data_processed))
# 
# #creating an x_sequential column that is 1 through the number of SD points, which will make it easier to randomly choose one point
# SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
#   mutate(X_sequential = 1:nrow(SD_fixed_field_data_processed))
# 
# ## Extracting the soil data to the tree points 
# 
# #LM
# LM_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LM_soil_text, LM_fixed_field_data_processed) #extracting soil textures for each point value
# LM_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LM_other, LM_fixed_field_data_processed) #extracting the other soil variables for each point value
# LM_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LM_extra, LM_fixed_field_data_processed) #extracting the extra soil variables for each point value
# LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed, LM_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LM point dataframe
# LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LM point dataframe
# LM_fixed_field_data_processed_soils <- cbind(LM_fixed_field_data_processed_soils, LM_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LM point dataframe
# 
# #LC
# LC_soil_text_raster_250_data_pts <- raster::extract(soil_stack_LC_soil_text, LC_fixed_field_data_processed) #extracting soil textures for each point value
# LC_soil_other_raster_250_data_pts <- raster::extract(soil_stack_LC_other, LC_fixed_field_data_processed) #extracting the other soil variables for each point value
# LC_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_LC_extra, LC_fixed_field_data_processed) #extracting the extra soil variables for each point value
# LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed, LC_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
# LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
# LC_fixed_field_data_processed_soils <- cbind(LC_fixed_field_data_processed_soils, LC_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe
# 
# #SD
# SD_soil_text_raster_250_data_pts <- raster::extract(soil_stack_SD_soil_text, SD_fixed_field_data_processed) #extracting soil textures for each point value
# SD_soil_other_raster_250_data_pts <- raster::extract(soil_stack_SD_other, SD_fixed_field_data_processed) #extracting the other soil variables for each point value
# SD_soil_extra_raster_250_data_pts <- raster::extract(soil_stack_SD_extra, SD_fixed_field_data_processed) #extracting the extra soil variables for each point value
# SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed, SD_soil_text_raster_250_data_pts) #bind the soil textures data for each point to the LC point dataframe
# SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_other_raster_250_data_pts) #bind the other soil variable data for each point to the LC point dataframe
# SD_fixed_field_data_processed_soils <- cbind(SD_fixed_field_data_processed_soils, SD_soil_extra_raster_250_data_pts) #bind the extra soil variable data for each point to the LC point dataframe
# 
# 
# ## Creating Sandy and Clay/Loamy Available Water Columns ##
# 
# # Making four new soil metric variables for soil available water, which equals the field capacity - permanent wilting point
# 
# # LM
# 
# LM_fixed_field_data_processed_soils <- LM_fixed_field_data_processed_soils %>%
#   mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
#   mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
#   mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
#   mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm
# 
# # LC
# 
# LC_fixed_field_data_processed_soils <- LC_fixed_field_data_processed_soils %>%
#   mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
#   mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
#   mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
#   mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200)  # Clay/Loam Available Water 100-200 cm
# 
# # SD
# 
# SD_fixed_field_data_processed_soils <- SD_fixed_field_data_processed_soils %>%
#   mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
#   mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
#   mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
#   mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm

#### Creating the Function Comparing the soil vs. size values ####

# This function that generates real and randomly generated slopes from single linear regressions comparing tree 
# sizes to soil metrics for each population to see if there is a significant relationship between tree size and soil characteristics. 

# The function components are 
    # 1) generating an empty list for the slopes and the p-values from comparing the real and randomized slopes. 
    # 2) a loop to iterate over each population, over each soil metric, and over each size characteristics 
          # in which we generate the list of randomized slopes (1,000 permutations, shuffling the size metrics 
              # across the trees) and real slope, 
          # plot the histogram of randomized slopes and real slope,
          # generate the p-values (using significantly different),
          # and store the real slopes and p-values in their respective arrays

# The function returns the slope array, p-value array, and stored histogram plot list.

slopes_simulations <- function () {
  
  # Creating the slope and p-value arrays 
  
  Populations <- c("LM", "LC", "SD") #creating a list of population names
  Size.metrics <- c("SCA", "LCA", "CA", "CS", "DBH") #creating a list of size metric names
  Soil.metrics <- c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                     "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                     "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                    "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                    "Volume of water content -1500 kpa 100-200", 
                    "Nitrogen 0-5", "Nitrogen 100-200", 
                    "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                    "Sand Available Water 0-5", "Sand Available Water 100-200",
                    "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200") #creating a list of soil metric names
  #creating the slopes and p-value arrays [soil metrics x size metrics x populations]
  slopes_array <- array(NA, dim = c(length(Soil.metrics), length(Size.metrics), length(Populations)), 
                        list(Soil.metrics, Size.metrics, Populations)) #real slopes array
  p_values_array <- array(NA, dim = c(length(Soil.metrics), length(Size.metrics), length(Populations)), 
                        list(Soil.metrics, Size.metrics, Populations)) # p-value array

  #to store generated histogram plots
  plot_list <- list()   
  
  #starting a loop iterating over populations
  for (i in 1:length(Populations)) { 
   
     #assigning the population based on the current tree
     if (Populations[i] == "LM"){  #LM trees
       Population = "LM"
       dataframe_soils = LM_fixed_field_data_processed_soils #assigning the LM tree/soils dataframe
     } else if (Populations[i] == "LC"){ #LC trees
       Population = "LC"
       dataframe_soils = LC_fixed_field_data_processed_soils #assigning the LC tree/soils dataframe
     } else if (Populations[i] == "SD"){ #SD trees
       Population = "SD"
       dataframe_soils = SD_fixed_field_data_processed_soils #assigning the SD tree/soils dataframe
     }

    #creating a dataframe of just the soil characteristics
    fixed_field_data_processed_soils.condensed <-  st_drop_geometry(dataframe_soils) #dropping the spatial geometry to easily remove columns
    fixed_field_data_processed_soils.condensed <- fixed_field_data_processed_soils.condensed[, c(c("clay.content.0.5", "clay.content.100.200", "silt.0.5", "silt.100.200", "sand.0.5", "sand.100.200",
                                                                                                   "ph_0.5", "ph_100.200",  "vol_water_.10_0.5",
                                                                                                   "vol_water_.10_100.200", "vol_water_0.5",
                                                                                                   "vol_water_100.200", "vol_water_.1500kPa_0.5", 
                                                                                                   "vol_water_.1500_100.200", 
                                                                                                   "nitrogen.0.5", "nitrogen.100.200", 
                                                                                                   "SOC.0.5", "SOC.100.200",
                                                                                                   "sandy_avail_water_0.5", "sandy_avail_water_100.200",
                                                                                                   "clay_loam_avail_water_0.5", "clay_loam_avail_water_100.200"))] #dropping all columns except the relevant soil metrics
   
    
    #creating a dataframe of just the size variables
    fixed_field_data_processed_size_variables <- st_drop_geometry(dataframe_soils) #dropping the spatial geometry to easily remove columns
    fixed_field_data_processed_size_variables <- fixed_field_data_processed_size_variables[, c("Canopy_short", "Canopy_long",
                                                                                               "Canopy_area","Crown_spread", 
                                                                                               "DBH_ag")] #dropping all columns except the relevant size metrics
    
    #iterating over the number of soil variables
    for (j in 1:(ncol(fixed_field_data_processed_soils.condensed))){ 
      
      #iterating over the number of soil variables
      for (k in 1:(ncol(fixed_field_data_processed_size_variables))) {
        
        #extracting slopes from comparing soil values with randomized shape/size values with linear regressions
        
        #creating empty list to collect randomized slope values
        slopes <- c() 
        
        #creating the list of randomized slopes
        set.seed(21) #setting the seed
        for (y in 1:1000){ #for 1000 permutations
          fixed_field_data_processed_soils_shuffled <- transform(dataframe_soils, variable.shuffled = sample(fixed_field_data_processed_size_variables[,k])) #create a data frame with shuffled soil metrics across trees for the respective population
          lm <- lm(fixed_field_data_processed_soils_shuffled$variable.shuffled ~ fixed_field_data_processed_soils.condensed[,j]) #creating the single linear regression of the shuffle size metric vs. soil metrics for the population's trees
          lm_sum <- summary(lm) #extracting the linear regression information
          slopes <- c(slopes, lm_sum$coefficients[2]) #add the current p-value from the randomized single linear regressions to the list of stored slopes
        }
        
        #extracting the real slope from our points
        lm_real <- lm(fixed_field_data_processed_size_variables[,k] ~ fixed_field_data_processed_soils.condensed[,j]) #creating the linear regression of the size metric vs. soil metrics for the respective population
        lm_real_sum <- summary(lm_real) #extract the summary 
        lm_real_slope <- lm_real_sum$coefficients[2] #storing the slope from the lienar regression summary
        
        #creating the histogram plot of the randomly distributed p-values and our real slope
        Size_Variable = colnames(fixed_field_data_processed_soils.condensed)[j] #extracting the size variable name
        plot_out <- ggplot()+
          geom_histogram(aes(x=slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
          geom_vline(xintercept=lm_real_slope, col = "red") + #line of our real slope
          xlab(paste(Population, "Slopes of Shuffled", Size_Variable, "vs. our", Size_Variable)) +
          ylab("Frequency") +
          theme_classic()
        
        # store the histogram plot list with a descriptive name 
        plot_name <- paste(Populations[i], #population name
                           colnames(fixed_field_data_processed_soils.condensed)[j], #soil metric name
                           colnames(fixed_field_data_processed_size_variables)[k], #size metric name
                           sep = "_")
        plot_list[[plot_name]] <- plot_out
        
        #calculating pseudo p-value based on where real slope overlaps with histogram of randomized slopes
        
        # if using greater than hypothesis (which we did not use in this case)
        
          # total = 0  #set empty value
          # for (p in 1:length(slopes)){ #loop that adds 1 to the value total if the simulated slope value is less than our average value for our trees
          #   if (slopes[p] > lm_real_slope){
          #     total = total + 1
          #   }
          # } #add number of values of in the random set of slope values that are less than our mean slope
          # p_value <- (total / length(slopes)) #the proportion of random ANNs that are less than our slope
        
        # using the significantly different alternative hypothesis 
        
        p_value_greater_than <- sum(slopes >= lm_real_slope)/length(slopes)   # proportion of simulated slopes higher than our real slope
        p_value_less_than <- sum(slopes <= lm_real_slope)/length(slopes)   # proportion of simulated slopes lower than our real slope
        p_value <- min(1, 2 * min(p_value_greater_than, p_value_less_than)) # find twice the smaller tail (the "more extreme" one) and if the value is over 1, use 1
        
        #storing the results
        
        slopes_array[j, k, i] = lm_real_slope #assigning the each index in the array with the new slope value
        p_values_array[j, k, i] = p_value #assigning the each index in the array with the new p value

       # print(paste("Updating:", j, k, i)) #soil metric, size metric, population
       # print(paste("Real slope:", lm_real_slope))
       # print(paste("P-value:", p_value))
        
      }
    }
    
  }
  
  #returning the slope array, p-value array, and histogram plot list
  return(list(slopes_array = slopes_array, p_values_array = p_values_array, plot_list = plot_list))
  
}

#### Running the Function and Storing the Results ####

#running the simulation function
slope_simultations <- slopes_simulations()

#Example of extracting one of the histograms comparing the slopes for our original soil vs. size metrics to the shuffled ones
slope_simultations$plot_list$LM_clay.content.0.5_Canopy_short
plot <- slope_simultations$plot_list$LM_clay.content.0.5_Canopy_short #storing the example plot
plot +  #example of changing the x-axis label
  xlab("LM Slopes of Linear Regressions with Shuffled Clay Content (0-5 cm) vs Known Clay Content (0-5 cm)")

#if you want to see all of the plots at once run: 
    #slope_simultations$plot_list

#turning the simulation results into a data frame

slopes_df <- as.data.frame.table(slope_simultations$slopes_array, 
                                 responseName = "slope") #creating a dataframe from the slope array
pvals_df <- as.data.frame.table(slope_simultations$p_values_array, 
                                responseName = "p_value") #creating a dataframe from the p value array

# merging the slopes and p-values dataframes 
size.pop.slopes.df <- merge(slopes_df, pvals_df, by = c("Var3", "Var2", "Var1")) #merging the two dataframes into one
names(size.pop.slopes.df) <- c("Population", "Size.Variable", "Soil.Metric", "Slope", "P.value") #re-naming the columns to be more appropriate

# Bonferroni correction of the p-values because of multiple testing
size.pop.slopes.df <- size.pop.slopes.df %>%
  mutate(p_bonf_corrected = p.adjust(size.pop.slopes.df$P.value, method = "bonferroni"))

#creating a column for whether the p values are significant or not
size.pop.slopes.df <- size.pop.slopes.df %>%
  mutate(Significance.reg = case_when(P.value < 0.05 ~ "Y",
                                  P.value >= 0.05 ~ "N")) %>% #using uncorrected p-values
  mutate(Significance.bonf = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                       p_bonf_corrected >= 0.05 ~ "N")) #using Bonferonni corrected p-values

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
  
#### Summarizing the results ####

#summarizing the results data frame 
summary(size.pop.slopes.df)

## Tables

#table of number of yes or no by population 
table(size.pop.slopes.df$Population, size.pop.slopes.df$Significance.bonf)

# table of number of yes or nos for each size variables 

#LM
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"])

# table of number of yes or nos for each soil variables 

#LM
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"])


# table of number of yes or nos for each size variables and each soil variable 

#LM 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LM"])

#LC 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LC"])

#SD 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="SD"])

## plots to visualize the results

#plot of population, yes/no, and the interaction between size/shape and soil values
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Size.Variable,Soil.Metric)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, andsize/shape 
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Size.Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and soil values
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Soil.Metric)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()


## Heat Maps

#setting population to a factor and ensuring the appropriate levels
size.pop.slopes.df$Population <- as.factor(size.pop.slopes.df$Population)
levels(size.pop.slopes.df$Population)

#changing the number of digits included to be 3
options(digits=3) 

#across all populations heat map of Bonferonni corrected P-values 
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Values",
       title = "All Sites Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) 
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  #coord_flip()

#across LM, heat map of Bonferonni corrected P-values 
size.pop.slopes.df.LM <- size.pop.slopes.df[size.pop.slopes.df$Population == "LM", ] #isolating the LM p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P.value < 0.05, P.value, NA))) 

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, round(Slope, 3), NA)), color = "black")



#across LC, heat map of Bonferonni corrected P-values 
size.pop.slopes.df.LC <- size.pop.slopes.df[size.pop.slopes.df$Population == "LC", ] #isolating the LC p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)))

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, round(Slope, 3), NA)))


#across SD
size.pop.slopes.df.SD <- size.pop.slopes.df[size.pop.slopes.df$Population == "SD", ] #isolating the SD p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P Values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)))

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(P.value < 0.05, P.value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P.value < 0.05, round(Slope, 3), NA))) 

