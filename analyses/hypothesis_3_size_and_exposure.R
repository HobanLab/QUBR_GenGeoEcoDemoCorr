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

#creating the aspect and slope rasters

#load in xyz ASCII from INEGI on elevation
elevation_xyz <- read.table("./data/ASCII Elevation Inegi/conjunto_de_datos/f12b43b4_ms.xyz")
View(elevation_xyz)
plot(elevation_xyz)
elevation_sf <-st_as_sf(elevation_xyz, 
                        coords = c("long", "lat"), crs = 4326)

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
#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)

#transformations of variables (log, square root, ) for linear models

#creating columns with transformations: logged all of the variables
LM_fixed_field_data_processed_log <- LM_fixed_field_data_processed %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
LM_fixed_field_data_processed_sqrt <- LM_fixed_field_data_processed %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))


### Sizes vs. Elevation ###

# LM 

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression

LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_lm_sca_elev, aes(x= LM_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
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

#creating the linear regression

LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

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
ggplot(data = LM_fixed_field_data_processed_sqrt, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_log$Canopy_area_lg ~ LM_fixed_field_data_processed_log$Elevation..m.FIXED)

#linear regression with square root transformation of canopy area
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_sqrt$Canopy_area_sqrt ~ LM_fixed_field_data_processed_sqrt$Elevation..m.FIXED)


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
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

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
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with logged transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed_log$DBH_ag_lg ~ LM_fixed_field_data_processed_log$Elevation..m.FIXED)

#linear regression with square root transformation of aggregated DBH
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed_sqrt$DBH_ag_sqrt ~ LM_fixed_field_data_processed_sqrt$Elevation..m.FIXED)


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

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_sca_elev, aes(x= LC_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_elev, aes(sample = LC_lm_sca_elev$residuals))+
  geom_qq()

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
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression

LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long ~ LC_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_lca_elev, aes(x= LC_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_lca_elev, aes(sample = LC_lm_lca_elev$residuals))+
  geom_qq()

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
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed_log$Canopy_area_lg ~ LC_fixed_field_data_processed_log$Elevation..m.)

#linear regression with square root transformation of canopy area
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed_sqrt$Canopy_area_sqrt ~ LC_fixed_field_data_processed_sqrt$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CA_elev, aes(x= LC_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CA_elev, aes(sample = LC_lm_CA_elev$residuals))+
  geom_qq()

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

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread ~ LC_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_CS_elev, aes(x= LC_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_CS_elev, aes(sample = LC_lm_CS_elev$residuals))+
  geom_qq()

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
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed_log$DBH_ag_lg ~ LC_fixed_field_data_processed_log$Elevation..m.)

#linear regression with square root transformation of aggregated DBH
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed_sqrt$DBH_ag_sqrt ~ LC_fixed_field_data_processed_sqrt$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_lm_DBH_elev, aes(x= LC_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_lm_DBH_elev, aes(sample = LC_lm_DBH_elev$residuals))+
  geom_qq()

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

#short canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")


#creating the linear regression

SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short ~ SD_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_sca_elev, aes(x= SD_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_elev, aes(sample = SD_lm_sca_elev$residuals))+
  geom_qq()

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
cor.test(SD_fixed_field_data_processed$Elevation..m., SD_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#creating the linear regression

SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long ~ SD_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_lca_elev, aes(x= SD_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_lca_elev, aes(sample = SD_lm_lca_elev$residuals))+
  geom_qq()

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
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m., y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#creating the linear regression
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area ~ SD_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed_log$Canopy_area_lg ~ SD_fixed_field_data_processed_log$Elevation..m.)

#linear regression with square root transformation of canopy area
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed_sqrt$Canopy_area_sqrt ~ SD_fixed_field_data_processed_sqrt$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CA_elev, aes(x= SD_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CA_elev, aes(sample = SD_lm_CA_elev$residuals))+
  geom_qq()

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
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#creating the linear regression

SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread ~ SD_fixed_field_data_processed$Elevation..m.)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_CS_elev, aes(x= SD_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_CS_elev, aes(sample = SD_lm_CS_elev$residuals))+
  geom_qq()

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
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#creating the linear regression
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag ~ SD_fixed_field_data_processed$Elevation..m.)

#linear regression with logged transformation of aggregated DBH
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed_log$DBH_ag_lg ~ SD_fixed_field_data_processed_log$Elevation..m.)

#linear regression with square root transformation of aggregated DBH
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed_sqrt$DBH_ag_sqrt ~ SD_fixed_field_data_processed_sqrt$Elevation..m.)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_lm_DBH_elev, aes(x= SD_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(SD_lm_DBH_elev, aes(sample = SD_lm_DBH_elev$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = SD_lm_DBH_elev, aes(x = SD_lm_DBH_elev$fitted.values, y = SD_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test visible in summary of the lm
summary(SD_lm_DBH_elev)


