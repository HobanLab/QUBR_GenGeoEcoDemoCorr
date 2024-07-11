#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(spdep) # to use morna's I functions like lag.listw
library(ape) # for computing the Moran's I stat
library(raster)


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
View(fixed_field_data_processed_sf_trans_coordinates)

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
View(LM_fixed_field_data_processed)
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

LM_fixed_field_data_processed %>%
  mutate(Elevation..m.FIXED = case_when((Elevation..m. > 700) ~ Elevation..m.*0.3048, 
                                        (Elevation..m. < 700) ~ Elevation..m.)) %>%
  ggplot() +
  geom_sf(aes(color = Elevation..m.FIXED))
  


View(LM_fixed_field_data_processed)

ggplot()+
  geom_sf(data = LM_fixed_field_data_processed, aes(color = Elevation..m.FIXED))



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

LM_box <- st_bbox(river_LM_trans)
LC_box <- st_bbox(river_LC_trans)
SD_box <- st_bbox(river_SD_trans)

##Load in environmental rasters 

#loading in soil textures
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

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)

### Sizes vs. Elevation ###

# LM 

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Short Canopy Axis")

LM_fixed_field_data_processed$Elevation..m.

#creating the linear regression

LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short ~ LM_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_elev, aes(x= LM_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_sca_elev, aes(sample = LM_lm_sca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_sca_elev, aes(x = LM_lm_sca_elev$fitted.values, y = LM_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_sca_elev)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Long Canopy Axis")

#creating the linear regression

LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_lca_elev, aes(x= LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_lca_elev, aes(sample = LM_lm_lca_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_lca_elev, aes(x = LM_lm_lca_elev$fitted.values, y = LM_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression

LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CA_elev, aes(x= LM_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CA_elev, aes(sample = LM_lm_CA_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CA_elev)


#crown spread

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread ~ LM_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_CS_elev, aes(x= LM_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_CS_elev, aes(sample = LM_lm_CS_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_CS_elev, aes(x = LM_lm_CS_elev$fitted.values, y = LM_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_CS_elev)


#DBH

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression

LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_DBH_elev, aes(x= LM_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_lm_DBH_elev, aes(sample = LM_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_lm_DBH_elev, aes(x = LM_lm_DBH_elev$fitted.values, y = LM_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(LM_lm_DBH_elev)

