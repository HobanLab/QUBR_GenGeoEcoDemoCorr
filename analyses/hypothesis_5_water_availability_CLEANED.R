# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Investigating to see if Quercus brandegeei tree shapes are influenced by Water Availability Proxies%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to determine if the trees' distance to river (a proxy for water availability)
# has a significant relationship with the size/shape of trees (short canopy axis, long canopy axis, canopy area, crown spread, DBH)
# A significant relationship could indicate how distance to river may influence the size and shape of the 
# trees and potentially help explain their distribution. 

# To test this, we used Simple Linear Regressions for each population, using distance to river as the explanatory 
# variables and each size metric as the response variable. We used a slope test to check for a significant relationship
# between distance the size/shape of trees.

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
        # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
                #generating river and point buffers and bounding boxes,
        # Extracting the distance to the river of each tree for each population,
# 3) Making a function for creating the single linear regressions, finding the best models (with transformations and removal of outliers or not) 
# that meet the conditions for SLRs
# 4) Running the function and storing the outputs

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
library(lmtest) #to use the Breuch-Pagan Test

# loading in the processed tree data 
source("./analyses/Data_Processing_Script.R")
# 
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
# #creating columns with transformations: logged all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_lg = log(Canopy_short))%>%
#   mutate(Canopy_long_lg = log(Canopy_long))%>%
#   mutate(Canopy_area_lg = log(Canopy_area))%>%
#   mutate(Crown_spread_lg = log(Crown_spread))%>%
#   mutate(DBH_ag_lg = log(DBH_ag))
# 
# #creating columns with transformations: square root all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
#   mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
#   mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
#   mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
#   mutate(DBH_ag_sqrt = sqrt(DBH_ag))
# 
# #creating columns with transformations: inverse all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_inv = (1/Canopy_short))%>%
#   mutate(Canopy_long_inv = (1/Canopy_long))%>%
#   mutate(Canopy_area_inv = (1/Canopy_area)) %>%
#   mutate(Crown_spread_inv = (1/Crown_spread))%>%
#   mutate(DBH_ag_inv = (1/DBH_ag))
# 
# #create dataframe with X and Y UTM coordinates
# fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with seperate x and y columns from the UTM 12N transformation
# fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
#   cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe
# 
# #export the csv of the UTM 12N points for using the file in ArcGIS to make new shapefiles
# fixed_field_data_processed_sf_trans_coordinates_download <- write.csv(fixed_field_data_processed_sf_trans_coordinates, "./data/fixed_field_data_processed_sf_trans_coordinates.csv", row.names = F)
# 
# #transforming the response variables variables (log, square root, inverse) for the linear models
# 
# #creating columns with transformations: log of all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_lg = log(Canopy_short))%>%
#   mutate(Canopy_long_lg = log(Canopy_long))%>%
#   mutate(Canopy_area_lg = log(Canopy_area))%>%
#   mutate(Crown_spread_lg = log(Crown_spread))%>%
#   mutate(DBH_ag_lg = log(DBH_ag))
# 
# #creating columns with transformations: square root of all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_sqrt = sqrt(Canopy_short))%>%
#   mutate(Canopy_long_sqrt = sqrt(Canopy_long))%>%
#   mutate(Canopy_area_sqrt = sqrt(Canopy_area))%>%
#   mutate(Crown_spread_sqrt = sqrt(Crown_spread))%>%
#   mutate(DBH_ag_sqrt = sqrt(DBH_ag))
# 
# #creating columns with transformations: inverse of all of the variables
# fixed_field_data_processed_sf_transformed <- fixed_field_data_processed_sf_transformed %>%
#   mutate(Canopy_short_inv = (1/Canopy_short))%>%
#   mutate(Canopy_long_inv = (1/Canopy_long))%>%
#   mutate(Canopy_area_inv = (1/Canopy_area)) %>%
#   mutate(Crown_spread_inv = (1/Crown_spread))%>%
#   mutate(DBH_ag_inv = (1/DBH_ag))
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
# #creating buffers around the rivers 
# 
# #LM
# river_buffer_LM <- st_buffer(river_LM_trans, 100) #100 m buffer
# ggplot()+
#   geom_sf(data = river_buffer_LM)+
#   geom_sf(data = river_LM_trans)+
#   geom_sf(data = LM_fixed_field_data_processed_sf)
# 
# #LC
# river_buffer_LC<- st_buffer(river_LC_trans, 100) #230 m buffer
# ggplot()+
#   geom_sf(data = river_buffer_LC)+
#   geom_sf(data = river_LC_trans)+
#   geom_sf(data = LC_fixed_field_data_processed_sf)
# 
# #SD
# river_buffer_SD<- st_buffer(river_SD_trans, 70) #70 m buffer
# ggplot()+
#   geom_sf(data = river_buffer_SD)+
#   geom_sf(data = river_SD_trans)+
#   geom_sf(data = SD_fixed_field_data_processed_sf)
# 
# #creating bounding boxes for each population
# 
# #creating a boundary box of LM with the UTM 12 N min and max lat lon values and then turning it into a simple feature geometry
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
# 
# ## Creating the distance to river columns ##
# 
# ## Generating the distance to river rasters for each population
# 
# #LM
# 
# #turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
# river_LM_trans_points <- st_cast(river_LM_trans, "LINESTRING") #turning the polyline of the river into a linestring object
# river_LM_trans_point_raster <- st_rasterize(river_LM_trans_points) #creating a raster out of the river linestring object
# plot(river_LM_trans_point_raster) #plotting the river linestring object
# 
# #turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
# river_LM_buffer_trans_outline <- st_cast(river_buffer_LM, "LINESTRING") #turning the polygon of the river buffer into a linestring object
# river_buffer_LM_point_raster <- st_rasterize(river_LM_buffer_trans_outline) #creating a raster of river buffer linestring object
# plot(river_buffer_LM_point_raster) #plotting the river buffer linestring object
# 
# #generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
# river_buffer_LM_point_raster[is.na(river_buffer_LM_point_raster[])] <- 0  #making sure the cell that are not the river buffer linestring raster have a 0 value
# dist_near_river_buffer_LM <- dist_to_nearest(river_buffer_LM_point_raster, river_LM_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
# plot(dist_near_river_buffer_LM) #plotting the distance to river raster
# 
# #LC
# 
# #turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
# river_LC_trans_points <- st_cast(river_LC_trans, "LINESTRING") #turning the polyline of the river into a linestring object
# river_LC_trans_point_raster <- st_rasterize(river_LC_trans_points) #creating a raster out of the river linestring object
# plot(river_LC_trans_point_raster)
# 
# #turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
# river_buffer_LC_points <- st_cast(river_buffer_LC, "LINESTRING") #turning the polygon of the river buffer into a linestring object
# river_buffer_LC_point_raster <- st_rasterize(river_buffer_LC_points) #creating a raster of river buffer linestring object
# plot(river_buffer_LC_point_raster) #plotting the river buffer linestring object
# 
# #generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
# river_buffer_LC_point_raster[is.na(river_buffer_LC_point_raster[])] <- 0  #making sure the cells that are not part of the the river buffer raster have a 0 value
# dist_near_river_buffer_LC <- dist_to_nearest(river_buffer_LC_point_raster, river_LC_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
# plot(dist_near_river_buffer_LC) #not using inverse distance
# 
# #SD
# 
# #turning the river polygon into a linestring object and then into a raster, to be able to later calculate the distances
# river_SD_trans_points <- st_cast(river_SD_trans, "LINESTRING") #turning the polyline of the river into a linestring object
# river_SD_trans_point_raster <- st_rasterize(river_SD_trans_points) #creating a raster out of the river linestring object
# plot(river_SD_trans_point_raster)
# 
# #turning the river buffer polygon into a linestring object and then into a raster to be able to later calculate the distances
# river_buffer_SD_points <- st_cast(river_buffer_SD, "LINESTRING") #turning the polygon of the river buffer into a linestring object
# river_buffer_SD_point_raster <- st_rasterize(river_buffer_SD_points) #creating a raster of river buffer linestring object
# plot(river_buffer_SD_point_raster) #plotting the river buffer linestring object
# 
# #generating a distance to river raster with the distances of each cell in the buffer raster from the river edge points, whereby the river raster cells are set to a distance of 0 m
# river_buffer_SD_point_raster[is.na(river_buffer_SD_point_raster[])] <- 0  #making sure the cells that are not part of the the river buffer raster have a 0 value
# dist_near_river_buffer_SD <- dist_to_nearest(river_buffer_SD_point_raster, river_SD_trans_points, progress = T) #creating a raster of the distances of each cell in the buffer raster to the linestring object of the river polygon, this can take a while to run
# plot(dist_near_river_buffer_SD) #plotting the distance to river raster
# 
# 
# ## Making it so the cells in the distance raster within or overlapping with the river raster are assigned 1 
# 
# #LM
# 
# #Assigning points within and overlapping with the river to be "true"
# LM_points_intersects_river <- st_intersects(LM_fixed_field_data_processed, river_LM_trans, sparse = F) #creating a list of true or falses for whether points intersect the river shapefiles
# LM_fixed_field_data_processed_intersects_river <- cbind(LM_fixed_field_data_processed, LM_points_intersects_river) #binding the list of true or falses with the point data
# #printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
# ggplot()+
#   geom_sf(data=river_LM_trans)+
#   geom_sf(data=LM_fixed_field_data_processed)+
#   geom_sf(data=LM_fixed_field_data_processed_intersects_river, aes(color = LM_points_intersects_river))
# 
# #LC 
# 
# #Assigning points within and overlapping with the river to be "true"
# LC_points_intersects_river <- st_intersects(LC_fixed_field_data_processed, river_LC_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
# LC_fixed_field_data_processed_intersects_river <- cbind(LC_fixed_field_data_processed, LC_points_intersects_river) #binding the list of true or falses with the point data
# #printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
# ggplot()+
#   geom_sf(data=river_LC_trans)+
#   geom_sf(data=LC_fixed_field_data_processed)+
#   geom_sf(data=LC_fixed_field_data_processed_intersects_river, aes(color = LC_points_intersects_river))
# 
# #SD
# 
# #Assigning points within and overlapping with the river to be "true"
# SD_points_intersects_river <- st_intersects(SD_fixed_field_data_processed, river_SD_trans, sparse = F) #creating a list of true or falses for whether points intersect the rivershapefiles
# SD_fixed_field_data_processed_intersects_river <- cbind(SD_fixed_field_data_processed, SD_points_intersects_river) #binding the list of true or falses with the point data
# #printing the river polygon and the tree points, colored by whether are or aren't within or overlapping with the river
# ggplot()+
#   geom_sf(data=river_SD_trans)+
#   geom_sf(data=SD_fixed_field_data_processed)+
#   geom_sf(data=SD_fixed_field_data_processed_intersects_river, aes(color = SD_points_intersects_river))
# 
# ## Extracting distance to river for each tree using the distance to river raster
# 
# #LM
# LM_distance_data_pts <- st_extract(dist_near_river_buffer_LM, LM_fixed_field_data_processed) #extracting distance to river for each tree
# LM_fixed_field_data_processed_distance  <- cbind(LM_fixed_field_data_processed, LM_distance_data_pts) #binding the distance to river data for each point to the LM point dataframe
# 
# 
# #LC
# LC_distance_data_pts <- st_extract(dist_near_river_buffer_LC, LC_fixed_field_data_processed) #extracting distance to river for each tree
# LC_fixed_field_data_processed_distance  <- cbind(LC_fixed_field_data_processed, LC_distance_data_pts) #binding the distance to river data for each tree to the LC point dataframe
# 
# 
# #SD
# SD_distance_data_pts <- st_extract(dist_near_river_buffer_SD, SD_fixed_field_data_processed) #extracting distance to river for each tree
# SD_fixed_field_data_processed_distance  <- cbind(SD_fixed_field_data_processed, SD_distance_data_pts) #binding the distance to river data for each point to the SD point dataframe
# 
# ## Assigning all points within/overlapping river to distances of 0
# 
# #LM
# LM_fixed_field_data_processed_distance <- LM_fixed_field_data_processed_distance %>% 
#   mutate(d = case_when((LM_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
#                        (LM_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 
# 
# #LC
# LC_fixed_field_data_processed_distance <- LC_fixed_field_data_processed_distance %>%
#   mutate(d = case_when((LC_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
#                        (LC_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value 
# 
# #SD
# SD_fixed_field_data_processed_distance <- SD_fixed_field_data_processed_distance %>%
#   mutate(d = case_when((SD_points_intersects_river == T) ~ 0,  #assigns 0 to points within river
#                        (SD_points_intersects_river == F) ~ d)) #to points outside of river, it leaves the original distance value



#### Creating the Simple Linear Regression function ####

#comparing distance to river's edge to size values

# For all populations/each population and size/shape metric we created single variable regressions by...
#a) creating the single variable linear regressions with the un-transformed response variable, logged variable, and 
# square root of the variable, respectively, either with or without outliers
#b) testing which model best satisfies the conditions for the analysis: LINES
# Linearity, Independence, Normality of residuals, Equal variance of residuals, and simple random sample 
  # we tested Linearity by looking at the scatterplots,
  # we tested Independence by thinking about the explanatory and response variables across the points,
  # we tested Normality of Residuals using histograms, qq norm plots, and the Shapiro-Wilk's test,
  # we tested Equal Variance of Residuals using a fitted vs. residuals plot,
  # we tested Simple Random Sample by thinkg about the data collection method.
#c) We then looked for significant associations (slopes/correlations)
# 1) If the LINES conditions are met...
# we ran a slope test and a Pearson's correlation test to see if there is a significant association
# 2) If the LINES conditions are not met...
# we ran a Mann-Kendall test (non-parametric test) to look for a significant correlation/tau 



#Function for selecting the Simple Linear Regressions based on meeting the equal variance and normal residuals 
#conditions and then returns the model and slope test results

#The tests for Linearity, Independence, and Simple Random Sampling are not included in the function. Linearity is tested in the 
#"Running the Simple Linear Regressions" section and the other two conditions should be tested by the analyst.

simple_linear_regressions <- function(population, size_variable){ #input the population name and the size variable/response variable
  
  #assigning the population based on the inputted population
  if (population == "LM"){  #LM trees
    dataframe_distance = LM_fixed_field_data_processed_distance #assigning the LM tree/distance dataframe
  } else if (population == "LC"){ #LC trees
    dataframe_distance = LC_fixed_field_data_processed_distance #assigning the LC tree/distance dataframe
  } else if (population == "SD"){ #SD trees
    dataframe_distance = SD_fixed_field_data_processed_distance #assigning the SD tree/distance dataframe
  }
  
  #assigning the size/response variable based on the user input
  if (size_variable == "SCA"){  #Short Canopy Axis
    size_variable_name = "Canopy_short" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_short #assigning the short canopy axis variable to the size metric variable
  } else if (size_variable == "LCA"){ #Long Canopy Axis
    size_variable_name = "Canopy_long" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_long #assigning the long canopy axis variable to the size metric variable
  } else if (size_variable == "CA"){ #Canopy Area
    size_variable_name = "Canopy_area" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_area #assigning the canopy area variable to the size metric variable
  } else if (size_variable == "CS"){ #Crown Spread
    size_variable_name = "Crown_spread" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Crown_spread #assigning the crown spread variable to the size metric variable
  } else if (size_variable == "DBH"){ #DBH
    size_variable_name = "DBH_ag" #storing the name of the size variable we are using
    size_metric = dataframe_distance$DBH_ag #assigning the DBH variable to the size metric variable
  }
  
  #creating a dataframe with influential/outlier points removed 
  
  #using Cook's D to check for highly influential points that may skew the linear model results
  slr <- lm(size_metric ~ dataframe_distance$d) #creating a linear regression to use to calculate the Cook's D
  slr_cooks <- cooks.distance(slr) #creating a linear regression to use to calculate the Cook's D
  plot(slr_cooks, type = 'h') #checking to see which Cook's D are unusually high
  influential <- slr_cooks[(slr_cooks > (3 * mean(slr_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
  influential
  
  #removing points that were deemed too influential on the linear model fit
  dataframe_distance_no_outliers <- dataframe_distance[-c(as.numeric(names(influential))),]
  
  #creating the linear regressions
  
  #creating the base linear regression (no removal of outliers, no transformations)
  slr_dist_base  <- lm(size_metric ~ dataframe_distance$d) #generating the linear regression 

  #linear regression with transformations

  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_base_lg  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression

  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_base_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression 

  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_base_inv  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression 

  #linear regression with transformations and removal of outliers
  
  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_no_out_lg  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 

  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_no_out_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 

  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_no_out_inv  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 
 
  #finding and storing the results of the best performing model
  
  #making a list of all of the potential models
  linear_regressions_models <- list("Base Model" = slr_dist_base, # base model slope test p-value
                                    "Base Model with Log Transformation" = slr_dist_base_lg, 
                                    "Base Model with Square Root Transformation" = slr_dist_base_sqrt, 
                                    "Base Model with Inverse Transformation" = slr_dist_base_inv, #slope test p-values with transformations of the response variable
                                    "No Outliers Model with Log Transformation" = slr_dist_no_out_lg, 
                                    "No Outliers Model with Square Root Transformation" =slr_dist_no_out_sqrt, 
                                    "No Outliers Model with Inverse Transformation" =slr_dist_no_out_inv)  #slope test p-values with no outliers and transformations of the response variable
  
  #creating an empty list to store Shapiro-Wilks Test
  shapiro_wilk_list <- c()
  
  #creating an empty list to store Breusch-Pagan Test
  breusch_pagan_test_list <- c()
  
  #assigning a variable to "unmet" until one model meets the conditions
  conditions <- "unmet"
  
  #iterating over each model to test if it meets the normality and equal variance conditions and stopping when a model meets the conditions
  for (k in 1:length(linear_regressions_models)){
    #assigning the correct model based on the iteration
    test_model <- linear_regressions_models[[k]] 
    
    #checking the conditions, specifically the normality and equal variance conditions, assuming the linearity, independence, and simple random sample conditions are met
    
    #checking normality with Shapiro-Wilk Test
    shapiro_wilk <- shapiro.test(test_model$residuals)
    
    #appending the Shapiro-Wilks results in a list
    shapiro_wilk_list <- c(shapiro_wilk_list, shapiro_wilk$p.value)
    
    #Equal variance tests (Bartlett's Test)
    
    # Getting the data used in the model for the Bartlett's Test
    model_data <- model.frame(test_model)
    
    # Getting the response and explanatory variable names for the Bartlett's Test
    response_var <- names(model_data)[1] #all.vars(formula(chosen_model))[2]
    explanatory_var <- names(model_data)[2] #all.vars(formula(chosen_model))[3]
    
    #Breusch-Pagan Test for equal variances when data is normal, if significant, residuals do not have equal variance
    breusch_pagan_test <- bptest(as.formula(paste(response_var, " ~ ", explanatory_var)), # This will cause an error,
                                 data = model_data)
    
    #appending the Breusch-Pagan Test results in a list
    breusch_pagan_test_list <- c(breusch_pagan_test_list, breusch_pagan_test$p.value)
    
    #return test results of the model that meets the conditions
    
    #if the normal residuals and equal variance of residuals conditions ARE met
    if (shapiro_wilk$p.value > 0.05 & breusch_pagan_test$p.value > 0.05) { 
      
      #printing whether the normality and equal variance conditions were met
      print(paste(names(linear_regressions_models)[k], "Data met the Normality and Equal Variance Conditions."))
      
      #Storing the results of the slope test from the summary of the linear model
      slope_test_results <- test_model$coefficients #slope test results

      #assigning the chosen model index to a variable
      chose_model_index <- k
      
      #assigning a variable to "met" when one model met the conditions
      conditions <- "met"
      
      break #exiting the loop when both conditions have been met
      
    } else if (shapiro_wilk$p.value < 0.05 & breusch_pagan_test$p.value > 0.05) {
      print(paste(names(linear_regressions_models)[k], "Data did not meet the Normality Condition but met the Equal Variance Condition.")) #printing whether the normality and equal variance conditions were met
      next #move on to testing if the next model meets both conditions
    } else if (shapiro_wilk$p.value > 0.05 & breusch_pagan_test$p.value < 0.05) {
      print(paste(names(linear_regressions_models)[k], "Data met the Normality Condition but not the Equal Variance Condition.")) #printing whether the normality and equal variance conditions were met
      next #move on to testing if the next model meets both conditions
    }
  }
  
  #if the conditions are not met (such as the residuals were not normal and the variance not equal) we use the model with the closest to normal residuals
  if (conditions != "met") {
    for (k in 1:length(linear_regressions_models)){
      
      #we assign the test model be the one with the highest Shapriro-Wilks Test p-value, so the model with the residuals closest to normal
      test_model <- linear_regressions_models[[which.max(shapiro_wilk_list)]]

      #Storing the results of the slope test from the summary of the linear model
      slope_test_results <- test_model$coefficients #slope test results
      
      #assigning the chosen model index to a variable
      chose_model_index <- which.max(shapiro_wilk_list)
      
      #extracting the best model based on which one has the smallest p-value
      chosen_model <- test_model 
      
      #printing the chosen model
      chosen_model_statement <- print(paste("The chosen model is ", names(linear_regressions_models)[chose_model_index]))
      chosen_model_statement
      
      #printing that none of the models met the conditions
      lack_conditions_statement <- print(paste("None of the models met the normality and equal variance conditions, so we used the one that met the conditions the best"))
      lack_conditions_statement
      
    }
  } else { #if the condition is met
    
    #extracting the best model based on which one has the smallest p-value
    chosen_model <- test_model 
    
    #printing the chosen model
    chosen_model_statement <- print(paste("The chosen model is ", names(linear_regressions_models)[chose_model_index]))
    chosen_model_statement
    
    #printing that one of the models met the conditions
    lack_conditions_statement <- print(paste("Atleast one of the models met the normality and equal variance conditions"))
    lack_conditions_statement
    
  }
  
  #storing the model summary to extract values
  slope_test_summary <- summary(chosen_model) 
  
  #Storing the slope test t-value from the summary of the linear model
  slope_test_t_value <- slope_test_summary$coefficients["d", "t value"]
  
  #Storing the slope test p-value from the summary of the linear model
  slope_test_p_value <- slope_test_summary$coefficients["d", "Pr(>|t|)"]
  
  print(paste("Slope Test Results:", 
              "t =", slope_test_t_value, 
              "p =", slope_test_p_value))
  
  print(slope_test_summary)
  
  #returns the name of the model that performed the best, normality results, equal variance results, and the slope test results
  return(list("chosen_model" = chosen_model, #returning the chosen model
              "shapiro_wilk" = shapiro_wilk, #returning the Shapiro-Wilks Test for normal residual results
              "breusch_pagan_test" = breusch_pagan_test,#returning the Breusch-Pagan Test for equal variance of residuals results
              "slope_test_summary" = slope_test_summary, #returning the linear regression slope test results
              "slope_test_t_value" = slope_test_t_value, #returning the linear regression slope test results t-value
              "slope_test_p_value" = slope_test_p_value, #returning the linear regression slope test results p-value
              "chosen_model_statement" = chosen_model_statement, #returning a statement of which model was chosen
              "lack_conditions_statement" = lack_conditions_statement)) #returning a statement if none of the models met the conditions
  
}


#### Running the Simple Linear Regressions ####

#LM

#SCA

#running the simple linear regression function
simple_linear_regressions_LM_SCA <- simple_linear_regressions("LM", "SCA")
simple_linear_regressions_LM_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_SCA$chosen_model, aes(x= simple_linear_regressions_LM_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_SCA$chosen_model, aes(sample = simple_linear_regressions_LM_SCA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_SCA$chosen_model, aes(x = simple_linear_regressions_LM_SCA$chosen_model$fitted.values, y = simple_linear_regressions_LM_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LM_LCA <- simple_linear_regressions("LM", "LCA")
simple_linear_regressions_LM_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_LCA$chosen_model, aes(x= simple_linear_regressions_LM_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_LCA$chosen_model, aes(sample = simple_linear_regressions_LM_LCA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_LCA$chosen_model, aes(x = simple_linear_regressions_LM_LCA$chosen_model$fitted.values, y = simple_linear_regressions_LM_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LM_CA <- simple_linear_regressions("LM", "CA")
simple_linear_regressions_LM_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Square Root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CA$chosen_model, aes(x= simple_linear_regressions_LM_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CA$chosen_model, aes(sample = simple_linear_regressions_LM_CA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CA$chosen_model, aes(x = simple_linear_regressions_LM_CA$chosen_model$fitted.values, y = simple_linear_regressions_LM_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LM_CS <- simple_linear_regressions("LM", "CS")
simple_linear_regressions_LM_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CS$chosen_model, aes(x= simple_linear_regressions_LM_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CS$chosen_model, aes(x= simple_linear_regressions_LM_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CS$chosen_model, aes(x = simple_linear_regressions_LM_CS$chosen_model$fitted.values, y = simple_linear_regressions_LM_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LM_DBH <- simple_linear_regressions("LM", "DBH")
simple_linear_regressions_LM_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_DBH$chosen_model, aes(x= simple_linear_regressions_LM_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_DBH$chosen_model, aes(x= simple_linear_regressions_LM_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_DBH$chosen_model, aes(x = simple_linear_regressions_LM_DBH$chosen_model$fitted.values, y = simple_linear_regressions_LM_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#LC

#SCA

#running the simple linear regression function
simple_linear_regressions_LC_SCA <- simple_linear_regressions("LC", "SCA")
simple_linear_regressions_LC_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_SCA$chosen_model, aes(x= simple_linear_regressions_LC_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_SCA$chosen_model, aes(x= simple_linear_regressions_LC_SCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_SCA$chosen_model, aes(x = simple_linear_regressions_LC_SCA$chosen_model$fitted.values, y = simple_linear_regressions_LC_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LC_LCA <- simple_linear_regressions("LC", "LCA")
simple_linear_regressions_LC_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_LCA$chosen_model, aes(x= simple_linear_regressions_LC_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_LCA$chosen_model, aes(x= simple_linear_regressions_LC_LCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_LCA$chosen_model, aes(x = simple_linear_regressions_LC_LCA$chosen_model$fitted.values, y = simple_linear_regressions_LC_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LC_CA <- simple_linear_regressions("LC", "CA")
simple_linear_regressions_LC_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CA$chosen_model, aes(x= simple_linear_regressions_LC_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CA$chosen_model, aes(x= simple_linear_regressions_LC_CA$chosen_model$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CA$chosen_model, aes(x = simple_linear_regressions_LC_CA$chosen_model$fitted.values, y = simple_linear_regressions_LC_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LC_CS <- simple_linear_regressions("LC", "CS")
simple_linear_regressions_LC_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CS$chosen_model, aes(x= simple_linear_regressions_LC_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CS$chosen_model, aes(x= simple_linear_regressions_LC_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CS$chosen_model, aes(x = simple_linear_regressions_LC_CS$chosen_model$fitted.values, y = simple_linear_regressions_LC_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LC_DBH <- simple_linear_regressions("LC", "DBH")
simple_linear_regressions_LC_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_DBH$chosen_model, aes(x= simple_linear_regressions_LC_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_DBH$chosen_model, aes(x= simple_linear_regressions_LC_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_DBH$chosen_model, aes(x = simple_linear_regressions_LC_DBH$chosen_model$fitted.values, y = simple_linear_regressions_LC_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#SD

#SCA

#running the simple linear regression function
simple_linear_regressions_SD_SCA <- simple_linear_regressions("SD", "SCA")
simple_linear_regressions_SD_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_SCA$chosen_model, aes(x= simple_linear_regressions_SD_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_SCA$chosen_model, aes(x= simple_linear_regressions_SD_SCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_SCA$chosen_model, aes(x = simple_linear_regressions_SD_SCA$chosen_model$fitted.values, y = simple_linear_regressions_SD_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_SD_LCA <- simple_linear_regressions("SD", "LCA")
simple_linear_regressions_SD_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_LCA$chosen_model, aes(x= simple_linear_regressions_SD_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_LCA$chosen_model, aes(x= simple_linear_regressions_SD_LCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_LCA$chosen_model, aes(x = simple_linear_regressions_SD_LCA$chosen_model$fitted.values, y = simple_linear_regressions_SD_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_SD_CA <- simple_linear_regressions("SD", "CA")
simple_linear_regressions_SD_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CA$chosen_model, aes(x= simple_linear_regressions_SD_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CA$chosen_model, aes(x= simple_linear_regressions_SD_CA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CA$chosen_model, aes(x = simple_linear_regressions_SD_CA$chosen_model$fitted.values, y = simple_linear_regressions_SD_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_SD_CS <- simple_linear_regressions("SD", "CS")
simple_linear_regressions_SD_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CS$chosen_model, aes(x= simple_linear_regressions_SD_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CS$chosen_model, aes(x= simple_linear_regressions_SD_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CS$chosen_model, aes(x = simple_linear_regressions_SD_CS$chosen_model$fitted.values, y = simple_linear_regressions_SD_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_SD_DBH <- simple_linear_regressions("SD", "DBH")
simple_linear_regressions_SD_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_DBH$chosen_model, aes(x= simple_linear_regressions_SD_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance (m)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_DBH$chosen_model, aes(x= simple_linear_regressions_SD_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_DBH$chosen_model, aes(x = simple_linear_regressions_SD_DBH$chosen_model$fitted.values, y = simple_linear_regressions_SD_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance (m)")


