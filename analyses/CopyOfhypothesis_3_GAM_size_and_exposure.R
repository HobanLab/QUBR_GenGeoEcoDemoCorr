#After having trouble with multiple linear rgressions because of issues 
#with the normality conditions and some nervousness about linearity. 
#We decided to use Generalized Additive Models (GAMs)


#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(car) #to create added variable plots and to run levene's test for checking ANOVA conditions
library(stars) # to convert raster into stars
library(gdalUtilities) #to be able to use gdalwarp
library(mgcv) #needed for the gam function for generalized additive models
library(starsExtra) #to use dist_to_nearest
library(MuMIn) #to use the dredge function
library(rpart) #to use the function rpart to check recurissive binary
library(visreg) #package to be able to plot effects of categorical variables
library(gratia) #using the function smooth estimates
library(gridExtra) #way to arrange ggplots in one plot
library(plotly) #3d plotting
devtools::install_github("AckerDWM/gg3D") #3d plotting
library("gg3D") #3d plotting
library(mgcViz) #3d plotting
library(rgl) #3d plotting


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


#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")

LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(Elevation..m. = as.numeric(Elevation..m.))

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
  mutate(Elevation..m.FIXED = case_when((Elevation..m. == 360) ~ NA, 
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
river_LM_trans <- st_as_sf(st_transform(river_LM, crs = 26912))
river_LC_trans <- st_as_sf(st_transform(river_LC, crs = 26912))
river_SD_trans <- st_as_sf(st_transform(river_SD, crs = 26912))

#creating buffers around the rivers
river_buffer_LM <- st_buffer(river_LM_trans, 100) #100 m buffer
ggplot()+
  geom_sf(data = river_buffer_LM)+
  geom_sf(data = river_LM_trans)+
  geom_sf(data = LM_fixed_field_data_processed_sf)

river_buffer_LC<- st_buffer(river_LC_trans, 100) #230 m buffer
ggplot()+
  geom_sf(data = river_buffer_LC)+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_sf)

river_buffer_SD<- st_buffer(river_SD_trans, 70) #70 m buffer
ggplot()+
  geom_sf(data = river_buffer_SD)+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_sf)


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

#Importing the cropped rasters for LM, LC, and SD
CEM_15_utm_LM <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LM.tif"))
CEM_15_utm_LC <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LC.tif"))
CEM_15_utm_SD <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_SD.tif"))

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


## Extracting the aspect 

#all points 

#extracting the slope in degrees, using the queens method (neighbor = 8)
all_points_aspect_raster_15 <- terra::terrain(CEM_15_utm_all_points, unit = 'degrees', neighbors = 8, 'aspect')

#plot the slopes
ggplot()+
  geom_raster(data= as.data.frame(all_points_aspect_raster_15, xy = T), aes(x=x, y=y, fill = aspect))+
  geom_sf(data = fixed_field_data_processed_sf_trans_coordinates)+
  scale_fill_viridis_c()


#creating dataframes for each population and the slope and aspect data by extracting the slope and aspect data fromk each cell for each point and combining the data into a dataframe


#all points
all_points_aspect_raster_15_data_pts <- extract(all_points_aspect_raster_15, fixed_field_data_processed_sf_trans_coordinates) #extracting aspect for each point value
all_points_slope_raster_15_data_pts <- extract(all_points_slope_raster_15, fixed_field_data_processed_sf_trans_coordinates) #extracting slope for each point value
all_points_fixed_field_data_processed_terrain <- cbind(fixed_field_data_processed_sf_trans_coordinates, all_points_aspect_raster_15_data_pts) #bind the aspect data for each point to the LM point dataframe
all_points_fixed_field_data_processed_terrain <- cbind(all_points_fixed_field_data_processed_terrain, all_points_slope_raster_15_data_pts) #bind the slope data for each point to the LM point dataframe

#recategorizing the aspect data

#setting values of 360 to 0 

#all points
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts = case_when((all_points_aspect_raster_15_data_pts == "360") ~  0,
                                                          (all_points_aspect_raster_15_data_pts != "360")~ all_points_aspect_raster_15_data_pts))

# all points

# North, Northeast, East, Southeast, South, Southwest, West, Northwest

# the directions are a range of 45 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_8_categorical = case_when((all_points_aspect_raster_15_data_pts > 0 & all_points_aspect_raster_15_data_pts < 22.5) ~ "N",  #north is between 337.5 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 337.5 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N", #359.99999
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 67.5) ~ "NE", #northeast is between 22.5 and 67.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 67.5 & all_points_aspect_raster_15_data_pts < 112.5) ~ "E", #east is between 67.5 and 112.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 112.5 & all_points_aspect_raster_15_data_pts < 157.5) ~ "SE", #southeast is between 122.5 and 157.5
                                                                        (all_points_aspect_raster_15_data_pts >= 157.5 & all_points_aspect_raster_15_data_pts < 202.5) ~ "S", #south is between 157.5 and 202.5
                                                                        (all_points_aspect_raster_15_data_pts >= 202.5 & all_points_aspect_raster_15_data_pts < 247.5) ~ "SW", #southwest is between 202.5 and 246.5
                                                                        (all_points_aspect_raster_15_data_pts >= 247.5 & all_points_aspect_raster_15_data_pts < 292.5) ~ "W", #west is between 247.5 and 292.5 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 292.5 & all_points_aspect_raster_15_data_pts < 337.5) ~ "NW")) #northwest is between 292.5 and 337.5 degrees

# North, East, South, West

# the directions are a range of 90 degrees 
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  mutate(all_points_aspect_raster_15_data_pts_4_categorical = case_when((all_points_aspect_raster_15_data_pts >= 0 & all_points_aspect_raster_15_data_pts < 45) ~ "N",  #north is between 315 and 22.5
                                                                        (all_points_aspect_raster_15_data_pts >= 315 & all_points_aspect_raster_15_data_pts < 359.999999999999999) ~ "N",
                                                                        (all_points_aspect_raster_15_data_pts >= 22.5 & all_points_aspect_raster_15_data_pts < 135) ~ "E", #northeast is between 22.5 and 135  degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 135 & all_points_aspect_raster_15_data_pts < 225) ~ "S", #south is between 135 and 225 degrees
                                                                        (all_points_aspect_raster_15_data_pts >= 225 & all_points_aspect_raster_15_data_pts < 315) ~ "W")) #west is between 225 and 315

#descriptive summary for LM


### Generalized Additive Models ###

#using only the 8 categories

# all points 

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
all_points_fixed_field_data_processed_terrain_no_NA <- all_points_fixed_field_data_processed_terrain %>%
  filter(is.na(all_points_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(all_points_aspect_raster_15_data_pts_8_categorical) == F)

#Cook's D
plot(all_points_multiple_lm_SCA)
all_points_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_mlm_SCA_cooks <- cooks.distance(all_points_mlm_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (2 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#SCA


all_points_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA)
all_points_add.gam_SCA_interact <- gam(Canopy_short ~ Elevation..m.FIXED * all_points_slope_raster_15_data_pts * all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC, shows the smoothed model is the best fit
AIC(all_points_add.gam_SCA, all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed_first_term, 
    all_points_add.gam_SCA.smoothed_second_term, all_points_add.gam_SCA_interact)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#using different distributions that don't care about the normal distribution: quasi, poisson, quasi-poisson (in order of complexity)
all_points_add.gam_SCA.smoothed.quasi <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasi())
all_points_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA, family = poisson())
all_points_add.gam_SCA.smoothed.quasipoisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                    data = all_points_fixed_field_data_processed_terrain_no_NA, family = quasipoisson())

#we then used liklihood ratio tests to see which level of complexity fits the models the best
anova(all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed.quasi, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_SCA.smoothed.quasi, all_points_add.gam_SCA.smoothed.poisson, test = "LRT") #quasi vs. poisson
anova(all_points_add.gam_SCA.smoothed.quasi, all_points_add.gam_SCA.smoothed.quasipoisson, test = "LRT")  #quasi vs. quasipoisson
anova(all_points_add.gam_SCA.smoothed.poisson, all_points_add.gam_SCA.smoothed.quasipoisson, test = "LRT") #quasipoisson vs. poisson
#these likelihood ratio tests demonstrate that a poisson model is sufficient and a better fit compared  a quasi and quasipoisson model 

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed.poisson)


#comparing the model's the models GCV summary values to see which is lowest
summary(all_points_add.gam_SCA)
summary(all_points_add.gam_SCA.smoothed)
summary(all_points_add.gam_SCA.smoothed.poisson)

#we do not need to dredge the poisson model, but hear is the 
dredge <- dredge(all_points_add.gam_SCA.smoothed.poisson) #using the dredge model to narro the models down to the best choice
dredge[1,] #extracting the best model
all_points_add.gam_SCA.smoothed.dredged <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA)

#Chosen model: all_points_add.gam_SCA.smoothed.poisson

#updating K values, I did not in this scenario but if the k' and edf were close, we would raise the K 
all_points_add.gam_SCA.smoothed.poisson <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA)
k.check(all_points_add.gam_SCA.smoothed.poisson)
#after attempting to try different K values, the default values appear to work the best

plot(all_points_add.gam_SCA.smoothed.poisson, all.terms = T)
#par(mfrow = c(2,2))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Elevation (m)", ylab = expression(f[1]*'(Elevation)'))
plot.gam(all_points_add.gam_SCA.smoothed, xlab = "Slope (º)", ylab = "f_1 (Slope), 3.38")


# Extract smooth effects for Elevation
elev_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, select = "s(Elevation..m.FIXED)")

# Extract smooth effects for Slope
slope_effects <- smooth_estimates(all_points_add.gam_SCA.smoothed, select = "s(all_points_slope_raster_15_data_pts)")

# Plot Elevation Effect
p1 <- ggplot(elev_effects, aes(x = Elevation..m.FIXED, y = .estimate)) +
  #geom_smooth(se = T) +
  geom_line(color = "blue", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "blue", alpha = 0.2) +
  labs(x = "Elevation (m)", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Elevation") +
  theme_minimal()


# Plot Slope Effect
p2 <- ggplot(slope_effects, aes(x = all_points_slope_raster_15_data_pts, y = .estimate)) +
  geom_line(color = "darkgreen", linewidth = 1) +
  geom_ribbon(aes(ymin = .estimate - se, ymax = .estimate + se), fill = "darkgreen", alpha = 0.2) +
  labs(x = "Slope", y = "Effect on Short Canopy Axis", title = "Smooth Effect of Slope") +
  theme_minimal()

p3 <- visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = TRUE, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# Print the plots
grid.arrange(p1, p2, p3, ncol = 2)

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


#plotting with vis.gam
dev.off() #resetting the plot for a new plot
vis.gam(all_points_add.gam_SCA.smoothed, plot.type = "persp", theta = 25,  xlab = "Aspect", 
        ylab = "Elevation (m)")

#extracting the fitted values for the GAM for plotting the model
fitted_values_all_points_add.gam_SCA <- fitted.values(all_points_add.gam_SCA.smoothed)

devtools::install_github("AckerDWM/gg3D")
library("gg3D")
ggplot(all_points_fixed_field_data_processed_terrain_no_NA, aes(x=Elevation..m.FIXED, y=all_points_slope_raster_15_data_pts, 
                                                                z=Canopy_short, color=all_points_aspect_raster_15_data_pts_8_categorical)) + 
  theme_void() +
  axes_3D() +
  stat_3D() + 
  geom_smooth(method = "gam", formula = all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short ~ 
                all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts + 
                all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

