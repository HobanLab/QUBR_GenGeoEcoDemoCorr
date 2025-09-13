# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if the Quercus brandegeei populations have different mean soil metrics%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is determine whether there are significantly different mean
# soil metrics (clay content, pH, etc.) between the three Quercus brandegeei populations.
# A lack of differences could indicate which characteristics may explain why the trees are found at
# these distinct sites versus other locations. 
# Observed differences could may explain other differences in tree growth/survival between the sites, 
# as well as relate to environmental/topographic characteristics. 

# To test this, we used difference in means tests based on which conditions were met for each soil metric.
    # If the residuals were not normal, we used a Kruskal-Wallis test and Post-Hoc Wilcoxon Rank Sum Test
    # If the residuals were normal but the variance was NOT equal, we used a Welch's ANOVA test and Post-Hoc Tamhane's T2 Test
    # If the residuals were normal and the variance was equal, we used a Traditional ANOVA test and Post-Hoc Pairwise T-test

# The script is broken into sections of 
# 1) loading and processing the packages and spatial/size/shape data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations and loading in the river outline shapefiles, 
# 2) processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
#boxs around the  rivers, stacking the rasters for each population, and processing them into one dataframe for all and each population
# 3) Creating the four new soil metrics: Sandy available Water (0-5 and 100-200 cm) and Clay/Loam Available Water (0-5 and 100-200 cm)
# 4) Choosing random trees per grid cell for each population to avoid independence issues during tests. 
# 5) Making the function for checking conditions and running the appropriate difference in means test
# 6) Running the function and storing the outputs for each soil metric

#### Loading libraries and relevant data ####

library(tidyverse) # for graphing and data organization
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) # ggplot extension
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

# loading in the tree data (size, elevation, lat/lon, ID, size/shape)

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#creating a transformed point shapefile with UTM 12 N an equal area projection
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) 

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

# Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")


#Upload ArcGIS river shapefile and filter out polygons for each population

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
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed)

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
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(LM_fixed_field_data_processed))

#creating an x_sequential column that is 1 through the number of LC points, which will make it easier to randomly choose one point
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(LC_fixed_field_data_processed))

#creating an x_sequential column that is 1 through the number of SD points, which will make it easier to randomly choose one point
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(X_sequential = 1:nrow(SD_fixed_field_data_processed))

## Extracting the soil data to the tree points 

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

# Making four new soil metric variables for soil available water, which equals the field capacity - permanent wilting point

# LM

LM_fixed_field_data_processed_soils <- LM_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm

# LC

LC_fixed_field_data_processed_soils <- LC_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>%
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200)  # Clay/Loam Available Water 100-200 cm

# SD

SD_fixed_field_data_processed_soils <- SD_fixed_field_data_processed_soils %>%
  mutate(sandy_avail_water_0.5 = vol_water_0.5 - vol_water_.1500kPa_0.5) %>% # Sand Available Water 0-5 cm
  mutate(sandy_avail_water_100.200 = vol_water_100.200 - vol_water_.1500_100.200) %>% # Sand Available Water 100-200 cm
  mutate(clay_loam_avail_water_0.5 = vol_water_.10_0.5 - vol_water_.1500kPa_0.5) %>% # Clay/Loam Available Water 0-5 cm
  mutate(clay_loam_avail_water_100.200 = vol_water_.10_100.200 - vol_water_.1500_100.200) # Clay/Loam Available Water 100-200 cm

#### Choosing a Random Tree per Grid Cell ####

# For each population, randomly selecting one tree from each grid cell to avoid issues because 
# the rasters are only 250 m resolution and certain populations might have more points in specific 
# cells, skewing the mean results away from the potentially real mean

#LM

#creating a grid over the soil cells
LM_tree_grid_cropped <- st_make_grid(soil_stack_LM_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster to ensure they were made properly
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_LM_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = LM_tree_grid_cropped, fill = NA)

#selecting a point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, LM_fixed_field_data_processed_sf, sparse =T) #making sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
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
LM_list_grids_and_point_trees_df <- as.data.frame(unlist(LM_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(LM_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
#filters out grid cells that do not have trees within them
LM_list_grids_and_trees_fixed <- LM_list_grids_and_point_trees_df %>% 
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
LM_fixed_field_data_processed_trees_soils <- LM_fixed_field_data_processed_soils %>%
  filter(X %in% LM_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with the row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data= LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_fixed_field_data_processed_trees_soils, color = "red")


#LC

#creating a grid over the soil cells
LC_tree_grid_cropped <- st_make_grid(soil_stack_LC_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster to ensure they were made properly
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




#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
LC_list_grids_and_point_trees_df <- as.data.frame(unlist(LC_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(LC_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
#filters out grid cells that do not have trees within them

#### CHECK THAT THE TREE WE THINK WE ARE GRABBING IS ACTUALLY WHAT WE ARE GRABBING ####

LC_list_grids_and_trees_fixed <- LC_list_grids_and_point_trees_df %>% 
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  mutate(data_row = LC_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them


#filtering out point data to be just the focal points
LC_fixed_field_data_processed_trees_soils <- LC_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% LC_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with the row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data= LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_fixed_field_data_processed_trees_soils, color = "red")

#SD

#creating a grid over the soil cells
SD_tree_grid_cropped <- st_make_grid(soil_stack_SD_soil_text, cellsize = c(230, 265))

#plotting the grid over an example soil raster to ensure they were made properly
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


#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
SD_list_grids_and_point_trees_df <- as.data.frame(unlist(SD_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(SD_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
#filters out grid cells that do not have trees within them
SD_list_grids_and_trees_fixed <- SD_list_grids_and_point_trees_df %>% 
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
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

## Finalizing the tree soil metric dataframe

#combining the LM, LC, and SD tree dataframes with the soil metrics and randomly chosen points within each grid cell

fixed_field_data_processed_trees_soils <- rbind(LM_fixed_field_data_processed_trees_soils, LC_fixed_field_data_processed_trees_soils) #combining the LM and LC soil and randomly chosen tree data
fixed_field_data_processed_trees_soils <- rbind(fixed_field_data_processed_trees_soils, SD_fixed_field_data_processed_trees_soils) #combining the SD tree point data to the LM and LC soil and randomly chosen tree point data

#creating a column/variable with locality as a factor to be able to use it in the Tamhane's T2 Test later

fixed_field_data_processed_trees_soils$Locality_Factor <- as.factor(fixed_field_data_processed_trees_soils$Locality)

#### Making Function for Differences in Means ####

# Function that checks the conditions and then runs the appropriate difference in means test
    # If the residuals were not normal, we used a Kruskal-Wallis test and Post-Hoc Wilcoxon Rank Sum Test.
        # We also performed this test for ever soil metric because it is a non-parametric test allowing for comparisons across the soil metrics. 
    # If the residuals were normal but the variance was NOT equal, we used a Welch's ANOVA test and Post-Hoc Tamhane's T2 Test
    # If the residuals were normal and the variance was equal, we used a Traditional ANOVA test and Post-Hoc Pairwise T-test

# For checking the residuals were normal, we used a Shapiro-Wilks Test
# For checking equal variance, we check that test outcomes were both true, 1) rule of thumb (< 2 means equal variance) test and 
    # if the residuals were normal, 
          #we used the 2) Bartlett's test (p>0.05 means equal variance) to test if there was equal variance because it works better with normal residuals
    # if the residuals were NOT normal,
          #we used the 2) Fligner-Killeen test to test (p>0.05 means equal variance) if there was equal variance because it works better with non-normal residuals
# We also computed the Levene's Test, as an additional, less robust test for additional evidence, but we did not factor it into our difference in means test decision.

# The function returns:
  # a) an ANOVA model, 
  # b) the Shapiro-Wilks Test results, 
  # c) Fligner-Killeen Test Results,
  # d) Bartlett's Test results,
  # e) Levene's Test results,
  # f) Rule of Thumb Test results,
  # g) the chosen Difference in Means Test results,
  # h) the chosen Difference in Means Test Post-Hoc Test results,
  # i) a print of the difference in means tests chosen,
  # j) Kruskal-Wallis Test results,
  # k) andPost-Hoc Wilcoxon Rank Sum Test results.

mean_soil_function <- function(soil_group, data = fixed_field_data_processed_trees_soils, Populations = "Locality") {
  
  # Building the formula for the difference in means tests
  formula <- as.formula(paste(soil_group, "~", Populations)) #creating the formula of the soil group as the y variable and the locality as the x variable for expediting the following tests
  
  # Initial ANOVA test 
  anova_model <- aov(formula, data = data) 
  
  # checking to see if the residuals are normal
  
  shapiro_test <- shapiro.test(anova_model$residuals) #Shapiro-Wilks Test
  
  #Equal variance tests
  
  #Fligner-Killeen, more useful when data is not normal or there are outliers 
  fligner_test <- fligner.test(formula, data = fixed_field_data_processed_trees_soils)
  #bartlett's test for equal variances when data is normal, which in this case it is
  bartlett_test <- bartlett.test(formula, data = fixed_field_data_processed_trees_soils)
  #levene's test, not super robust to strong differences to normality
  levenes_test <- car::leveneTest(formula, data = fixed_field_data_processed_trees_soils)
  #rule of thumb test
  thumb_test <- tapply(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality, sd)
  thumb_test_results <- max(thumb_test, na.rm = T) / min(thumb_test, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass
  
  # checking conditions to choose with difference in means test to use
  if (shapiro_test$p.value < 0.05) { #if the residuals are NOT normally distributed
    #kruskal-Wallis test because the data has normally distributed residuals and equal variance of residuals
    test <- kruskal.test(formula, data = fixed_field_data_processed_trees_soils)
    #post-hoc Wilcoxon rank sum tests
    post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality,
                         p.adjust.method = "fdr") #p value adjusted using false discovery rate method
    #storing the tests used
    test_type = "kruskal-Wallis + Wilcoxon Rank Sum Test"
    #printing out which test was used
    print(paste("kruskal-Wallis Test with a Post-Hoc Wilcoxon Rank Sum Test"))
    
  } else if (shapiro_test$p.value > 0.05) { #if the residuals are normally distributed
    if (bartlett_test$p.value < 0.05 & thumb_test_results > 2) { #if the equal variance of residuals condition is NOT met
      
      #Welch's ANOVA, does not assume equal variances, but does meet normality condition
      test <- oneway.test(formula, data = fixed_field_data_processed_trees_soils, var.equal = F)
      #post hoc Welch's ANOVA test: Tamhane's T2 Test
      post_hoc <- tamhaneT2Test(fixed_field_data_processed_trees_soils[[soil_group]]~fixed_field_data_processed_trees_soils$Locality_Factor, data = fixed_field_data_processed_trees_soils)
      #storing the tests used
      test_type = "Welch's ANOVA + Tamhane's Test"
      #printing out which test was used
      print(paste("Welch's ANOVA with a Post-Hoc Tamhane's Test"))
      
    } else if (bartlett_test$p.value > 0.05 & thumb_test_results < 2) { #if the equal variance of residuals condition IS met
      #traditional ANOVA because equal variance and normality were met
      test <- anova(anova_model)
      #post-hoc Tukey's HSD
      post_hoc <- TukeyHSD(anova_model)
      #storing the tests used
      test_type = "ANOVA + Tukey's HSD"
      #printing out which test was used
      print(paste("ANOVA Test with a Post-Hoc Pairwise T-Test Test"))
    }
  }
  
  #kruskal-Wallis test because it is non-parametric
  kruskal_test <- kruskal.test(formula, data = fixed_field_data_processed_trees_soils)
  #post-hoc Wilcoxon rank sum tests
  kruskal_post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality,
                                   p.adjust.method = "fdr") #p value adjusted using false discovery rate method
  
  
  return(list(
    anova_model = summary(anova_model),
    shapiro_test = shapiro_test,
    fligner_test = fligner_test,
    bartlett_test = bartlett_test,
    levenes_test = levenes_test,
    thumb_test_results = thumb_test_results,
    final_test = test,
    posthoc = post_hoc,
    test_type = test_type,
    kruskal_test = kruskal_test, 
    kruskal_post_hoc = kruskal_post_hoc
  ))
}

#### Running the Function Results ####
  

##clay 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_0.5 <- mean_soil_function("clay.content.0.5")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_clay_0.5$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_clay_0.5$posthoc

#Kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_0.5_mean_p.value <- mean_soil_function_clay_0.5$kruskal_test$p.value

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_0_5 <- aov(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content vs. Population")


##clay 100-200 

#running the Difference in Means Analysis function
mean_soil_function_clay_100.200 <- mean_soil_function("clay.content.100.200")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_clay_100.200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_100.200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_100.200_mean_p.value <- mean_soil_function_clay_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_100_200 <- aov(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")


#silt 0-5

#running the Difference in Means Analysis function
mean_soil_function_silt_0.5 <- mean_soil_function("silt.0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_silt_0.5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_silt_0.5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_silt_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
silt_0.5_mean_p.value <- mean_soil_function_silt_0.5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_silt_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_silt_0_5 <- aov(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_silt_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

##silt 100-200

#running the Difference in Means Analysis function
mean_soil_function_silt_100.200 <- mean_soil_function("silt.100.200")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_silt_100.200$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_silt_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_silt_100.200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
silt_100.200_mean_p.value <- mean_soil_function_silt_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_silt_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_silt_100_200 <- aov(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_silt_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

##sand  0-5 

#running the Difference in Means Analysis function
mean_soil_function_sand_0.5 <- mean_soil_function("sand.0.5")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_sand_0.5$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_sand_0.5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sand_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sand_0.5_mean_p.value <- mean_soil_function_sand_0.5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sand_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_sand_0_5 <- aov(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sand_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

## sand 100-200

#running the Difference in Means Analysis function
mean_soil_function_sand_100.200 <- mean_soil_function("sand.100.200")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_sand_100.200$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_sand_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sand_100.200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sand_100.200_mean_p.value <- mean_soil_function_sand_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sand_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_sand_100_200 <- aov(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sand_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

## ph 0-5

#running the Difference in Means Analysis function
mean_soil_function_ph_0_5 <- mean_soil_function("ph_0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_ph_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_ph_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
ph_0.5_mean_p.value <- mean_soil_function_ph_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_ph_0_5 <- aov(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_ph_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for pH at 0-5 cm vs. Population")

##ph 100-200

#running the Difference in Means Analysis function
mean_soil_function_ph_100_200 <- mean_soil_function("ph_0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_ph_100_200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_ph_100_200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
ph_100.200_mean_p.value <- mean_soil_function_ph_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_ph_100_200 <- aov(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_ph_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for pH at 100-200 cm vs. Population")

##soil organic carbon 0-5

#running the Difference in Means Analysis function
mean_soil_function_SOC_0_5 <- mean_soil_function("SOC.0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_SOC_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_SOC_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
soc_0.5_mean_p.value <- mean_soil_function_SOC_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_soc_0_5 <- aov(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_soc_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Organic Carbon at 0-5 cm vs. Population")

#soil organic carbon 100-200

#running the Difference in Means Analysis function
mean_soil_function_SOC_100_200 <- mean_soil_function("SOC.100.200")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_SOC_100_200$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_SOC_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_SOC_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
soc_100.200_mean_p.value <- mean_soil_function_SOC_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_soc_100_200 <- aov(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_soc_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Organic Carbon at 100-200 cm vs. Population")

#volume of water content at -10 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_10_0_5 <- mean_soil_function("vol_water_.10_0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_10_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_10_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_10kpa_0.5_mean_p.value <- mean_soil_function_vol_water_10_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_10_0.5 <- aov(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_10_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -10 kpa at 0-5 cm vs. Population")

#volume of water content at -10 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_10_100_200 <- mean_soil_function("vol_water_.10_100.200")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_10_100_200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_10_100_200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_10kpa_100.200_mean_p.value <- mean_soil_function_vol_water_10_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_10_100.200 <- aov(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_10_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -10 kpa at 100-200 cm vs. Population")

#volume of water content at -33 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_33_0_5 <- mean_soil_function("vol_water_0.5")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_33_0_5$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_vol_water_33_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_33_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_33kpa_0.5_mean_p.value <- mean_soil_function_vol_water_33_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_33_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_33_0.5 <- aov(vol_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_33_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -33 kpa at 0-5 cm vs. Population")

#volume of water content at -33 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_33_100_200 <- mean_soil_function("vol_water_100.200")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_33_100_200$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_vol_water_33_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_33_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_33kpa_100.200_mean_p.value <- mean_soil_function_vol_water_33_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_33_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_33_100.200 <- aov(vol_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_33_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -33 kpa at 100-200 cm vs. Population")

#volume of water content at -1500 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_1500_0_5 <- mean_soil_function("vol_water_.1500kPa_0.5")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_1500_0_5$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_vol_water_1500_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_1500_0_5$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_1500kpa_0.5_mean_p.value <- mean_soil_function_vol_water_1500_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500kPa_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_1500_0.5 <- aov(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_1500_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -1500 kpa at 0-5 cm vs. Population")

#volume of water content at -1500 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_1500_100_200 <- mean_soil_function("vol_water_.1500_100.200")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_1500_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_1500_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_1500kpa_100.200_mean_p.value <- mean_soil_function_vol_water_1500_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_1500_100.200 <- aov(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_1500_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -1500 kpa at 100-200 cm vs. Population")

#nitrogen 0-5

#running the Difference in Means Analysis function
mean_soil_function_nitrogen_0_5 <- mean_soil_function("nitrogen.0.5")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_nitrogen_0_5$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_nitrogen_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
nitrogen_0.5_mean_p.value <- mean_soil_function_nitrogen_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_nitrogen_0.5 <- aov(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_nitrogen_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Nitrogen Content at 0-5 cm vs. Population")

# nitrogen 100-200

#running the Difference in Means Analysis function
mean_soil_function_nitrogen_100_200 <- mean_soil_function("nitrogen.100.200")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_nitrogen_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_nitrogen_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
nitrogen_100.200_mean_p.value <- mean_soil_function_nitrogen_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_nitrogen_100.200 <- aov(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_nitrogen_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Nitrogen Content at 100-200 cm vs. Population")

# sandy available water 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_sandy_avail_water_0_5 <- mean_soil_function("sandy_avail_water_0.5")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_sandy_avail_water_0_5$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_sandy_avail_water_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sandy_avail_water_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sandy_avail_water_0.5_mean_p.value <- mean_soil_function_sandy_avail_water_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sandy_avail_water_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sandy_avail_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_sandy_avail_water_0.5 <- aov(sandy_avail_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sandy_avail_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Sand Available Water at 0-5 cm vs. Population")

# sandy available water 100-200 cm

#running the Difference in Means Analysis function
mean_soil_function_sandy_avail_water_100_200 <- mean_soil_function("sandy_avail_water_100.200")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_sandy_avail_water_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sandy_avail_water_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sandy_avail_water_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sandy_avail_water_100.200_mean_p.value <- mean_soil_function_sandy_avail_water_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sandy_avail_water_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sandy_avail_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_sandy_avail_water_100.200 <- aov(sandy_avail_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sandy_avail_water_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Sand Available Water at 100-200 cm vs. Population")

# clay loam available water 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_loam_avail_water_0_5 <- mean_soil_function("clay_loam_avail_water_0.5")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_clay_loam_avail_water_0_5$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_loam_avail_water_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_loam_avail_water_0.5_mean_p.value <- mean_soil_function_clay_loam_avail_water_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay_loam_avail_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_loam_avail_water_0.5 <- aov(clay_loam_avail_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_loam_avail_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

# clay loam available water 100-200 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_loam_avail_water_100_200 <- mean_soil_function("clay_loam_avail_water_100.200")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_clay_loam_avail_water_100_200$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_clay_loam_avail_water_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_loam_avail_water_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_loam_avail_water_100.200_mean_p.value <- mean_soil_function_clay_loam_avail_water_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay_loam_avail_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_loam_avail_water_100.200 <- aov(clay_loam_avail_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_loam_avail_water_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 100-200 cm vs. Population")

#### Presenting Results ####

#Heat Map 

#Storing the p-values from the chosen Difference in Means test in a vector
p_value_mean <- c(clay_0.5_mean_p.value, clay_100.200_mean_p.value, silt_0.5_mean_p.value,
                    silt_100.200_mean_p.value,
                    sand_0.5_mean_p.value,
                    sand_100.200_mean_p.value,
                    ph_0.5_mean_p.value,
                    ph_100.200_mean_p.value,
                    soc_0.5_mean_p.value,
                    soc_100.200_mean_p.value,
                    vol_wat_10kpa_0.5_mean_p.value,
                    vol_wat_10kpa_100.200_mean_p.value,
                    vol_wat_33kpa_0.5_mean_p.value,
                    vol_wat_33kpa_100.200_mean_p.value,
                    vol_wat_1500kpa_0.5_mean_p.value,
                    vol_wat_1500kpa_100.200_mean_p.value,
                    nitrogen_0.5_mean_p.value,
                    nitrogen_100.200_mean_p.value,
                    sandy_avail_water_0.5_mean_p.value,
                    sandy_avail_water_100.200_mean_p.value,
                    clay_loam_avail_water_0.5_mean_p.value,
                    clay_loam_avail_water_100.200_mean_p.value
)



# Bonferroni correction of the p-values because of multiple testing
p_bonf_corrected <- p.adjust(p_value_mean, method = "bonferroni")

#creating empty dataframe for inputting the function into
random_pop.df <- data.frame("Shape.Size" = rep(c("Clay 0-5 cm", "Clay 100-200 cm", "Silt 0-5 cm", "Silt 100-200 cm", "Sand 0-5 cm", "Sand 100-200 cm", #column of the soil metric names
                                                 "pH 0-5 cm", "pH 100-200 cm", "Soil Organic Carbon 0-5 cm", "Soil Organic Carbon 100-200 cm", 
                                                 "Volume of water content -10 kPa 0-5 cm", "Volume of water content -10 kPa 100-200 cm",
                                                 "Volume of water content -33 kPa 0-5 cm", "Volume of water content -33 kPa 100-200 cm",
                                                 "Volume of water content -1500 kPa 0-5 cm", "Volume of water content -1500 kPa 100-200 cm", 
                                                 "Nitrogen 0-5 cm", "Nitrogen 100-200 cm", "Sand Available Water 0-5 cm", "Sand Available Water 100-200 cm",
                                                 "Clay/Loam Available Water 0-5 cm", "Clay/Loam Available Water 100-200 cm")),
                            "P_Value" = p_bonf_corrected, #Bonferonni-corrected p-values
                            "Significance" = c(rep(NA, 22))) #whether the p-values are significant (p<0.05) or not

#creating the significance column based on the significance of the p-values
random_pop.df <- random_pop.df %>%
  mutate(Significance = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                  p_bonf_corrected >= 0.05 ~ "N"))

#Turning off scientific notation format
options(scipen=999)

#Creating a heat map in ggplot of the p-values for each soil metric (y-axis) and whether they are 
#significant or not (x-axis), with significant p-values labeled
ggplot(aes(x = fct_reorder(Shape.Size, P_Value), y = Significance, fill = P_Value), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Difference Between Mean Soil Metrics Between Populations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + #setting color pallete
  geom_text(aes(label = ifelse(P_Value < 0.001, "< 0.001", NA)), col = "white") + #labeling significant p-values if less than 0.001 as "<0.001"
  geom_text(aes(label = ifelse(P_Value < 0.5 & P_Value > 0.001, round(P_Value, 8), NA)), col = "white") + #labeling cells with significant p-values less then 0.05 and greater than 0.001, rounding to 8 decimal places
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))





