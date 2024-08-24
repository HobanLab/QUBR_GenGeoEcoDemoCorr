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
library(geostatsp) 
library(tmaptools)

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


###creating the inverse distance

#LM

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LM_trans_outline <- st_cast(river_LM_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #create raster of lake edge points
plot(river_LM_trans_point_raster)

river_LM_buffer_trans_outline <- st_cast(river_buffer_LM, "LINESTRING") #turns the polyline of the river into a multipoint object
river_buffer_LM_point_raster <- st_rasterize(river_LM_buffer_trans_outline) #create raster of lake edge points
plot(river_buffer_LM_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, river_LM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
#dist_near_river_buffer_LM_inverse <- 1/dist_near_river_buffer_LM #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LM) #not using inverse distance


#LC

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_LC_trans_points <- st_cast(river_LC_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_LC_trans_point_raster <- st_rasterize(river_LC_trans_points) #create raster of lake edge points
plot(river_LC_trans_point_raster)

river_buffer_LC_points <- st_cast(river_buffer_LC, "LINESTRING") #turns the polyline of the river buffer into a multipoint object in stars
river_buffer_LC_point_raster <- st_rasterize(river_buffer_LC_points) #create raster of lake edge points
plot(river_buffer_LC_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_LC_point_raster[is.na(river_buffer_LC_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_LC <- dist_to_nearest(river_buffer_LC_point_raster, river_LC_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
#dist_near_river_buffer_LC_inverse <- 1/dist_near_river_buffer_LC #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_LC) #not using inverse distance



#SD

#turning river polygon into multipoints and then into a raster for using them to calculate the distances
river_SD_trans_points <- st_cast(river_SD_trans, "LINESTRING") #turns the polyline of the river into a multipoint object
river_SD_trans_point_raster <- st_rasterize(river_SD_trans_points) #create raster of lake edge points
plot(river_SD_trans_point_raster)

river_buffer_SD_points <- st_cast(river_buffer_SD, "LINESTRING") #turns the polyline of the river buffer into a multipoint object
river_buffer_SD_point_raster <- st_rasterize(river_buffer_SD_points) #create raster of lake edge points
plot(river_buffer_SD_point_raster)

#making a stars object of the distances of each cell in the buffer raster from the river edge points
river_buffer_SD_point_raster[is.na(river_buffer_SD_point_raster[])] <- 0  #making sure the points that are not the river buffer have a 0 value
dist_near_river_buffer_SD <- dist_to_nearest(river_buffer_SD_point_raster, river_SD_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the multipoints on the river polygon, this took an hour to run
#dist_near_river_buffer_SD_inverse <- 1/dist_near_river_buffer_SD #creating the inverse of the distance raster so that the higher values are closer to the river and the values are between 0-1
plot(dist_near_river_buffer_SD) #not using inverse distance


#making it so points within or overlapping with the river are assigned 1 

#Assigning points within and overlapping with the river to have true

#LM

LM_points_intersects_river <- st_intersects(LM_fixed_field_data_processed, river_LM_trans, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
LM_fixed_field_data_processed_intersects_river <- cbind(LM_fixed_field_data_processed, LM_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_LM_trans)+
  geom_sf(data=LM_fixed_field_data_processed)+
  geom_sf(data=LM_fixed_field_data_processed_intersects_river, aes(color = LM_points_intersects_river))

#LC 

LC_points_intersects_river <- st_intersects(LC_fixed_field_data_processed, river_LC_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
LC_fixed_field_data_processed_intersects_river <- cbind(LC_fixed_field_data_processed, LC_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_LC_trans)+
  geom_sf(data=LC_fixed_field_data_processed)+
  geom_sf(data=LC_fixed_field_data_processed_intersects_river, aes(color = LC_points_intersects_river))

#SD

SD_points_intersects_river <- st_intersects(SD_fixed_field_data_processed, river_SD_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
SD_fixed_field_data_processed_intersects_river <- cbind(SD_fixed_field_data_processed, SD_points_intersects_river) #binding the list of true or falses with the point data

ggplot()+
  geom_sf(data=river_SD_trans)+
  geom_sf(data=SD_fixed_field_data_processed)+
  geom_sf(data=SD_fixed_field_data_processed_intersects_river, aes(color = SD_points_intersects_river))



###extract distance to river for each point 

#all points


#LM
LM_distance_data_pts <- st_extract(dist_near_river_buffer_LM, LM_fixed_field_data_processed) #extracting aspect for each point value
LM_fixed_field_data_processed_distance  <- cbind(LM_fixed_field_data_processed, LM_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


#LC
LC_distance_data_pts <- st_extract(dist_near_river_buffer_LC, LC_fixed_field_data_processed) #extracting aspect for each point value
LC_fixed_field_data_processed_distance  <- cbind(LC_fixed_field_data_processed, LC_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


#SD
SD_distance_data_pts <- st_extract(dist_near_river_buffer_SD, SD_fixed_field_data_processed) #extracting aspect for each point value
SD_fixed_field_data_processed_distance  <- cbind(SD_fixed_field_data_processed, SD_distance_data_pts) #bind the aspect data for each point to the LM point dataframe


### Assigning all points within/overlapping river to distances of 0

#LM
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>% 
  mutate(d = case_when((LM_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LM_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 


#LC
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(d = case_when((LC_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (LC_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 


#SD

SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(d = case_when((SD_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
                       (SD_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value


#transformations of variables (log, square root, ) for linear models


#LM

#creating columns with transformations: logged all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

#LC

#creating columns with transformations: logged all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

#SD

#creating columns with transformations: logged all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))

#creating columns with transformations: square root all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
  mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
  mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
  mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
  mutate(DBH_ag_sqrt = sqrt(DBH_ag))

#creating columns with transformations: inverse all of the variables
SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
  mutate(Canopy_short_inv = (1/Canopy_short))%>%
  mutate(Canopy_long_inv = (1/Canopy_long))%>%
  mutate(Canopy_area_inv = (1/Canopy_area)) %>%
  mutate(Crown_spread_inv = (1/Crown_spread))%>%
  mutate(DBH_ag_inv = (1/DBH_ag))

#Simple Linear Regressions

#comparing distance to river's edge to size values

#LM

#SCA

#removing outliers

#Cook's D
LM_slr_SCA <- lm(Canopy_short ~ d, data = LM_fixed_field_data_processed_distance)
LM_slr_SCA_cooks <- cooks.distance(LM_slr_SCA) #calculating the cook.s D for each point
plot(LM_slr_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_SCA_cooks[(LM_slr_SCA_cooks > (2 * mean(LM_slr_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_distance_sca_no_outliers <- LM_fixed_field_data_processed_distance[-c(43, 45, 46, 60, 87, 88,  100, 123, 119, 148, 149, 151, 152, 154, 160, 165, 173, 189, 190, 195, 204, 206, 208, 209, 214),]


#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Short Canopy Axis (m)")

#creating the linear regression
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance$Canopy_short ~ LM_fixed_field_data_processed_distance$d)

#linear regression with transformations

#logged transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance$Canopy_short_lg ~ LM_fixed_field_data_processed_distance$d)

#square root transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance$Canopy_short_sqrt ~ LM_fixed_field_data_processed_distance$d)

#inverse transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance$Canopy_short_inv ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations and removal of outliers

#logged transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance_sca_no_outliers$Canopy_short_lg ~ LM_fixed_field_data_processed_distance_sca_no_outliers$d)

#square root transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance_sca_no_outliers$Canopy_short_sqrt ~ LM_fixed_field_data_processed_distance_sca_no_outliers$d)

#inverse transformations
LM_slr_dist_sca  <- lm(LM_fixed_field_data_processed_distance_sca_no_outliers$Canopy_short_inv ~ LM_fixed_field_data_processed_distance_sca_no_outliers$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_slr_dist_sca, aes(x= LM_slr_dist_sca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_slr_dist_sca, aes(sample = LM_slr_dist_sca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_slr_dist_sca$residuals) #significantly normal when square root transformation used and no removal of outliers


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Sqrt(Short Canopy Axis (m))")

                                                                                                 
#checking equal variance
ggplot(data = LM_slr_dist_sca, aes(x = LM_slr_dist_sca$fitted.values, y = LM_slr_dist_sca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LM_slr_dist_sca)


#LCA

#Cook's D
LM_slr_LCA <- lm(Canopy_long ~ d, data = LM_fixed_field_data_processed_distance)
LM_slr_LCA_cooks <- cooks.distance(LM_slr_LCA) #calculating the cook.s D for each point
plot(LM_slr_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_LCA_cooks[(LM_slr_LCA_cooks > (3 * mean(LM_slr_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_distance_lca_no_outliers <- LM_fixed_field_data_processed_distance[-c(43, 45, 50, 87,  103, 119, 124, 151, 152, 165, 173, 186, 190, 204, 206, 208, 212, 213, 214),]


#checking linearity 

LM_fixed_field_data_processed_distance$d

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Long Canopy Axis (m)")

#creating the linear regression
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance$Canopy_long ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance$Canopy_long_lg ~ LM_fixed_field_data_processed_distance$d)

#square root transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance$Canopy_long_sqrt ~ LM_fixed_field_data_processed_distance$d)

#inverse transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance$Canopy_long_inv ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations and removal of outliers

#logged transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance_lca_no_outliers$Canopy_long_lg ~ LM_fixed_field_data_processed_distance_lca_no_outliers$d)

#square root transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance_lca_no_outliers$Canopy_long_sqrt ~ LM_fixed_field_data_processed_distance_lca_no_outliers$d)

#inverse transformations
LM_slr_dist_lca  <- lm(LM_fixed_field_data_processed_distance_lca_no_outliers$Canopy_long_inv ~ LM_fixed_field_data_processed_distance_lca_no_outliers$d)


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance_lca_no_outliers, (aes(x=d, y=Canopy_long_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Long Canopy Axis (m)")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_slr_dist_lca, aes(x= LM_slr_dist_lca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_slr_dist_lca, aes(sample = LM_slr_dist_lca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_slr_dist_lca$residuals) #significantly normal when square root transformation and outliers are removed



#checking equal variance
ggplot(data = LM_slr_dist_lca, aes(x = LM_slr_dist_lca$fitted.values, y = LM_slr_dist_lca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LM_slr_dist_lca)

#CA


#Cook's D
LM_slr_CA <- lm(Canopy_area ~ d, data = LM_fixed_field_data_processed_distance)
LM_slr_CA_cooks <- cooks.distance(LM_slr_CA) #calculating the cook.s D for each point
plot(LM_slr_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_slr_CA_cooks[(LM_slr_CA_cooks > (1.5 * mean(LM_slr_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_processed_distance_ca_no_outliers <- LM_fixed_field_data_processed_distance[-c(43, 45, 50, 60, 87, 88, 100, 119, 124, 151, 152, 
                                                                                                   159, 165, 173, 175, 186, 189, 190, 195, 204, 206, 208, 212, 214),]


#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#creating the linear regression
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance$Canopy_area ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance$Canopy_area_lg ~ LM_fixed_field_data_processed_distance$d)

#square root transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance$Canopy_area_sqrt ~ LM_fixed_field_data_processed_distance$d)

#inverse transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance$Canopy_area_inv ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations and removal of outliers

#logged transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance_ca_no_outliers$Canopy_area_lg ~ LM_fixed_field_data_processed_distance_ca_no_outliers$d)

#square root transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance_ca_no_outliers$Canopy_area_sqrt ~ LM_fixed_field_data_processed_distance_ca_no_outliers$d)

#inverse transformations
LM_slr_dist_ca  <- lm(LM_fixed_field_data_processed_distance_ca_no_outliers$Canopy_area_inv ~ LM_fixed_field_data_processed_distance_ca_no_outliers$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_slr_dist_ca, aes(x= LM_slr_dist_ca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_slr_dist_ca, aes(sample = LM_slr_dist_ca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_slr_dist_ca, aes(x = LM_slr_dist_ca$fitted.values, y = LM_slr_dist_ca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LM_slr_dist_ca)


#CS

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#creating the linear regression
LM_slr_dist_cs  <- lm(LM_fixed_field_data_processed_distance$Crown_spread ~ LM_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LM_slr_dist_cs  <- lm(LM_fixed_field_data_processed_distance$Crown_spread_lg ~ LM_fixed_field_data_processed_distance$d)

#square root transformations
LM_slr_dist_cs  <- lm(LM_fixed_field_data_processed_distance$Crown_spread_sqrt ~ LM_fixed_field_data_processed_distance$d)

#inverse transformations
LM_slr_dist_cs  <- lm(LM_fixed_field_data_processed_distance$Crown_spread_inv ~ LM_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_slr_dist_cs, aes(x= LM_slr_dist_cs$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_slr_dist_cs, aes(sample = LM_slr_dist_cs$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_slr_dist_cs$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_slr_dist_cs, aes(x = LM_slr_dist_cs$fitted.values, y = LM_slr_dist_cs$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LM_slr_dist_cs)

#DBH

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#creating the linear regression
LM_slr_dist_dbh  <- lm(LM_fixed_field_data_processed_distance$DBH_ag ~ LM_fixed_field_data_processed_distance$d)

#linear regression with transformations

#logged transformations
LM_slr_dist_dbh  <- lm(LM_fixed_field_data_processed_distance$DBH_ag_lg ~ LM_fixed_field_data_processed_distance$d)

#square root transformations
LM_slr_dist_dbh  <- lm(LM_fixed_field_data_processed_distance$DBH_ag_sqrt ~ LM_fixed_field_data_processed_distance$d)

#inverse transformations
LM_slr_dist_dbh  <- lm(LM_fixed_field_data_processed_distance$DBH_ag_inv ~ LM_fixed_field_data_processed_distance$d)

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_slr_dist_dbh, aes(x= LM_slr_dist_dbh$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_slr_dist_dbh, aes(sample = LM_slr_dist_dbh$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LM_slr_dist_dbh$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LM_slr_dist_dbh, aes(x = LM_slr_dist_dbh$fitted.values, y = LM_slr_dist_dbh$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LM_slr_dist_dbh)


#LC

#SCA

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#creating the linear regression
LC_slr_dist_sca  <- lm(LC_fixed_field_data_processed_distance$Canopy_short ~ LC_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LC_slr_dist_sca  <- lm(LC_fixed_field_data_processed_distance$Canopy_short_lg ~ LC_fixed_field_data_processed_distance$d)

#square root transformations
LC_slr_dist_sca  <- lm(LC_fixed_field_data_processed_distance$Canopy_short_sqrt ~ LC_fixed_field_data_processed_distance$d)

#inverse transformations
LC_slr_dist_sca  <- lm(LC_fixed_field_data_processed_distance$Canopy_short_inv ~ LC_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_slr_dist_sca, aes(x= LC_slr_dist_sca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_slr_dist_sca, aes(sample = LC_slr_dist_sca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LC_slr_dist_sca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LC_slr_dist_sca, aes(x = LC_slr_dist_sca$fitted.values, y = LC_slr_dist_sca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LC_slr_dist_sca)


#LCA

#checking linearity 

LC_fixed_field_data_processed_distance$d

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#creating the linear regression
LC_slr_dist_lca  <- lm(LC_fixed_field_data_processed_distance$Canopy_long ~ LC_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LC_slr_dist_lca  <- lm(LC_fixed_field_data_processed_distance$Canopy_long_lg ~ LC_fixed_field_data_processed_distance$d)

#square root transformations
LC_slr_dist_lca  <- lm(LC_fixed_field_data_processed_distance$Canopy_long_sqrt ~ LC_fixed_field_data_processed_distance$d)

#inverse transformations
LC_slr_dist_lca  <- lm(CM_fixed_field_data_processed_distance$Canopy_long_inv ~ LC_fixed_field_data_processed_distance$d)



#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_slr_dist_lca, aes(x= LC_slr_dist_lca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_slr_dist_lca, aes(sample = LC_slr_dist_lca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LC_slr_dist_lca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LC_slr_dist_lca, aes(x = LC_slr_dist_lca$fitted.values, y = LC_slr_dist_lca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LC_slr_dist_lca)

#CA

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#creating the linear regression
LC_slr_dist_ca  <- lm(LC_fixed_field_data_processed_distance$Canopy_area ~ LC_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LC_slr_dist_ca  <- lm(LC_fixed_field_data_processed_distance$Canopy_area_lg ~ LC_fixed_field_data_processed_distance$d)

#square root transformations
LC_slr_dist_ca  <- lm(LC_fixed_field_data_processed_distance$Canopy_area_sqrt ~ LC_fixed_field_data_processed_distance$d)

#inverse transformations
LC_slr_dist_ca  <- lm(CM_fixed_field_data_processed_distance$Canopy_area_inv ~ LC_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_slr_dist_ca, aes(x= LC_slr_dist_ca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_slr_dist_ca, aes(sample = LC_slr_dist_ca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LC_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LC_slr_dist_ca, aes(x = LC_slr_dist_ca$fitted.values, y = LC_slr_dist_ca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LC_slr_dist_ca)


#CS

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#creating the linear regression
LC_slr_dist_cs  <- lm(LC_fixed_field_data_processed_distance$Crown_spread ~ LC_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LC_slr_dist_cs  <- lm(LC_fixed_field_data_processed_distance$Crown_spread_lg ~ LC_fixed_field_data_processed_distance$d)

#square root transformations
LC_slr_dist_cs  <- lm(LC_fixed_field_data_processed_distance$Crown_spread_sqrt ~ LC_fixed_field_data_processed_distance$d)

#inverse transformations
LC_slr_dist_cs  <- lm(CM_fixed_field_data_processed_distance$Crown_spread_inv ~ LC_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_slr_dist_cs, aes(x= LC_slr_dist_cs$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_slr_dist_cs, aes(sample = LC_slr_dist_cs$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LC_slr_dist_cs$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LC_slr_dist_cs, aes(x = LC_slr_dist_cs$fitted.values, y = LC_slr_dist_cs$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LC_slr_dist_cs)

#DBH

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#creating the linear regression
LC_slr_dist_dbh  <- lm(LC_fixed_field_data_processed_distance$DBH_ag ~ LC_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
LC_slr_dist_dbh  <- lm(LC_fixed_field_data_processed_distance$DBH_ag_lg ~ LC_fixed_field_data_processed_distance$d)

#square root transformations
LC_slr_dist_dbh  <- lm(LC_fixed_field_data_processed_distance$DBH_ag_sqrt ~ LC_fixed_field_data_processed_distance$d)

#inverse transformations
LC_slr_dist_dbh  <- lm(CM_fixed_field_data_processed_distance$DBH_ag_inv ~ LC_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_slr_dist_dbh, aes(x= LC_slr_dist_dbh$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_slr_dist_dbh, aes(sample = LC_slr_dist_dbh$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(LC_slr_dist_dbh$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = LC_slr_dist_dbh, aes(x = LC_slr_dist_dbh$fitted.values, y = LC_slr_dist_dbh$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#Slope Test visible in summary of the lm
summary(LC_slr_dist_dbh)


#SD

#SCA

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#creating the linear regression
SD_slr_dist_sca  <- lm(SD_fixed_field_data_processed_distance$Canopy_short ~ SD_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
SD_slr_dist_sca  <- lm(SD_fixed_field_data_processed_distance$Canopy_short_lg ~ SD_fixed_field_data_processed_distance$d)

#square root transformations
SD_slr_dist_sca  <- lm(SD_fixed_field_data_processed_distance$Canopy_short_sqrt ~ SD_fixed_field_data_processed_distance$d)

#inverse transformations
SD_slr_dist_sca  <- lm(SD_fixed_field_data_processed_distance$Canopy_short_inv ~ SD_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_slr_dist_sca, aes(x= SD_slr_dist_sca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_slr_dist_sca, aes(sample = SD_slr_dist_sca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(SD_slr_dist_sca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = SD_slr_dist_sca, aes(x = SD_slr_dist_sca$fitted.values, y = SD_slr_dist_sca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(SD_slr_dist_sca)


#LCA

#checking linearity 

SD_fixed_field_data_processed_distance$d

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#creating the linear regression
SD_slr_dist_lca  <- lm(SD_fixed_field_data_processed_distance$Canopy_long ~ SD_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
SD_slr_dist_lca  <- lm(SD_fixed_field_data_processed_distance$Canopy_long_lg ~ SD_fixed_field_data_processed_distance$d)

#square root transformations
SD_slr_dist_lca  <- lm(SD_fixed_field_data_processed_distance$Canopy_long_sqrt ~ SD_fixed_field_data_processed_distance$d)

#inverse transformations
SD_slr_dist_lca  <- lm(SD_fixed_field_data_processed_distance$Canopy_long_inv ~ SD_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_slr_dist_lca, aes(x= SD_slr_dist_lca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_slr_dist_lca, aes(sample = SD_slr_dist_lca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(SD_slr_dist_lca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = SD_slr_dist_lca, aes(x = SD_slr_dist_lca$fitted.values, y = SD_slr_dist_lca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(SD_slr_dist_lca)

#CA

#checking linearity 


#Cook's D
SD_slr_CA <- lm(Canopy_area ~ d, data = SD_fixed_field_data_processed_distance)
SD_slr_CA_cooks <- cooks.distance(SD_slr_CA) #calculating the cook.s D for each point
plot(SD_slr_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_slr_CA_cooks[(SD_slr_CA_cooks > (3 * mean(SD_slr_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_processed_distance_ca_no_outliers <- SD_fixed_field_data_processed_distance[-c(24, 25, 26, 73, 85, 105, 
                                                                                                   110, 134, 135, 136, 144, 145, 
                                                                                                   146, 147, 152, 199, 206, 211),]
                                                                                                   
#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")


#creating the linear regression
SD_slr_dist_ca  <- lm(SD_fixed_field_data_processed_distance$Canopy_area ~ SD_fixed_field_data_processed_distance$d)

#linear regression with transformations

#logged transformations
SD_slr_dist_ca  <- lm(SD_fixed_field_data_processed_distance$Canopy_area_lg ~ SD_fixed_field_data_processed_distance$d)

#square root transformations
SD_slr_dist_ca  <- lm(SD_fixed_field_data_processed_distance$Canopy_area_sqrt ~ SD_fixed_field_data_processed_distance$d)

#inverse transformations
SD_slr_dist_ca  <- lm(SD_fixed_field_data_processed_distance$Canopy_area_inv ~ SD_fixed_field_data_processed_distance$d)


#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_slr_dist_ca, aes(x= SD_slr_dist_ca$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_slr_dist_ca, aes(sample = SD_slr_dist_ca$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(SD_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = SD_slr_dist_ca, aes(x = SD_slr_dist_ca$fitted.values, y = SD_slr_dist_ca$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#Slope Test visible in summary of the lm
summary(SD_slr_dist_ca)


#CS

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#creating the linear regression
SD_slr_dist_cs  <- lm(SD_fixed_field_data_processed_distance$Crown_spread ~ SD_fixed_field_data_processed_distance$d)

#linear regression with transformations

#logged transformations
SD_slr_dist_cs  <- lm(SD_fixed_field_data_processed_distance$Crown_spread_lg ~ SD_fixed_field_data_processed_distance$d)

#square root transformations
SD_slr_dist_cs  <- lm(SD_fixed_field_data_processed_distance$Crown_spread_sqrt ~ SD_fixed_field_data_processed_distance$d)

#inverse transformations
SD_slr_dist_cs  <- lm(SD_fixed_field_data_processed_distance$Crown_spread_inv ~ SD_fixed_field_data_processed_distance$d)

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_slr_dist_cs, aes(x= SD_slr_dist_cs$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_slr_dist_cs, aes(sample = SD_slr_dist_cs$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(SD_slr_dist_cs$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = SD_slr_dist_cs, aes(x = SD_slr_dist_cs$fitted.values, y = SD_slr_dist_cs$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#Slope Test visible in summary of the lm
summary(SD_slr_dist_cs)

#DBH


#Cook's D
SD_slr_dbh <- lm(DBH_ag ~ d, data = SD_fixed_field_data_processed_distance)
SD_slr_dbh_cooks <- cooks.distance(SD_slr_dbh) #calculating the cook.s D for each point
plot(SD_slr_dbh_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_slr_dbh_cooks[(SD_slr_dbh_cooks > (2 * mean(SD_slr_dbh_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_processed_distance_dbh_no_outliers <- SD_fixed_field_data_processed_distance[-c(24, 30, 73, 85, 88, 103, 105, 
                                                                                                   110, 120, 134, 135, 136, 141, 142, 144,
                                                                                                   146, 147, 204, 199, 206, 211),]

#checking linearity 


#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#creating the linear regression
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance$DBH_ag ~ SD_fixed_field_data_processed_distance$d)


#linear regression with transformations

#logged transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance$DBH_ag_lg ~ SD_fixed_field_data_processed_distance$d)

#square root transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance$DBH_ag_sqrt ~ SD_fixed_field_data_processed_distance$d)

#inverse transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance$DBH_ag_inv ~ SD_fixed_field_data_processed_distance$d)


#linear regression with transformations and removal of outliers

#logged transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance_dbh_no_outliers$DBH_ag_lg ~ SD_fixed_field_data_processed_distance_dbh_no_outliers$d)

#square root transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance_dbh_no_outliers$DBH_ag_sqrt ~ SD_fixed_field_data_processed_distance_dbh_no_outliers$d)

#inverse transformations
SD_slr_dist_dbh  <- lm(SD_fixed_field_data_processed_distance_dbh_no_outliers$DBH_ag_inv ~ SD_fixed_field_data_processed_distance_dbh_no_outliers$d)

#plotting the linear model in ggplot for SCA
ggplot(data = SD_fixed_field_data_processed_distance_dbh_no_outliers, (aes(x=d, y=DBH_ag_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("sqrt(DBH)")

#checking normality

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_slr_dist_dbh, aes(x= SD_slr_dist_dbh$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance (m)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_slr_dist_dbh, aes(sample = SD_slr_dist_dbh$residuals))+
  geom_qq()

#shaprio wilk test
shapiro.test(SD_slr_dist_dbh$residuals) #significantly not normal, except when outliers are removed

#checking equal variance
ggplot(data = SD_slr_dist_dbh, aes(x = SD_slr_dist_dbh$fitted.values, y = SD_slr_dist_dbh$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance (m)")

#Slope Test visible in summary of the lm
summary(SD_slr_dist_dbh)




