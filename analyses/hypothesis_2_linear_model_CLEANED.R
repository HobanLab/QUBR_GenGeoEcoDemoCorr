# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei compete or facilitate with one another%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by the distance to other individuals of the same species 
# either due to competition or facilitation. 
# If they are impacted by facilitation, we would expect closer trees would be bigger. 
# If they are impacted by competition, we would expect closer trees to be smaller. 
# To test this, we used also performed a linear regression to see if for focal trees, there was a 
# relationship between how much competition the trees face (based on the 
# size of the neighbors over their distance to the focal trees) and the size of the focal trees.

# The script is broken into sections of 
# 1) loading and processing the packages and spatial/size/shape data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations and loading in the river outline shapefiles, 
# 2) using linear regression to see if tree size seem related to local competition  


#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for casdulating the moments of each variable
library(sf) # for plotting spatial objects
library(spatstat) # to run the nndist function
library(spdep) # to use Moran's I functions like lag.listw
library(ape) # for computing the Moran's I stat
library(raster) #to use point distance
library(nlme) # linear mixed effect models
library(MuMIn) #to be able to use model.sel for fitting linear models with spatial autocorrelation
library(geoR) # to be able to use variograms with the lme, requires XQuartz from 
library(Kendall)# to use the Kendall's Tau test to look for non-parametric correlations in the data

# Make a function that is the opposite of the %in% function
`%notin%` <- Negate(`%in%`) 


# loading in the tree data (size, elevation, lat/lon, ID, size/shape)

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

# creating the point shapefiles of the tree locations for each population in UTM 12 N

#creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
#sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#creating a transformed point shapefile with UTM 12 N an equal area projection
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) 

#storing point shapefiles for the trees by population

LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sf()

LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sf()

SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sf()

#create dataframe with X and Y UTM coordinates

fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with separate x and y columns from the UTM 12N transformation
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
  cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

# creating a dataframe with the 5 average nearest neighbors (ANN) for each individual tree/row
fixed_field_data_processed_NN_UTM <- fixed_field_data_processed_sf_trans_coordinates %>%  #creates a dataframe with the ANN of the closest 5 individual trees for each individual
  mutate(dist1 = nndist(X = X.1, Y= Y, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
  mutate(dist2 = nndist(X = X.1, Y= Y, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
  mutate(dist3 = nndist(X = X.1, Y= Y, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
  mutate(dist4 = nndist(X = X.1, Y= Y, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
  mutate(dist5 = nndist(X = X.1, Y= Y, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
  rowwise()%>% #so that in the next part we take the averages across rows
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5))) # %>% #creates a column of the average distances (1-5) of each individual

# Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns

LM_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "SD")



#### Linear Model ####

# To see if trees that face more competition (they face closer and larger trees) are smaller, for each population,
  # 1) we create a bounding box around the points and then cropped the bounding box of the points by 20 m on all sides to avoid edge effects
  # 2) randomly select a tree (focal tree) from each grid cell (width/lenth of 40*mean DBH) with a tree in it
  # 3) filtering out trees to just the focal points and finding the trees that are within a buffer of the 
        #focal tree with a radius of 40 times the mean population DBH
  # 4) filtering and storing the focal trees without any neighbors
  # 5) iterating over a every focal tree in a loop and calculating the sum of the shape/size metrics of the neighbors for each focal tree (compeititon metric for each tree)
  # 6) calculating and use histograms descriptive statistics for the competition metrics 
  # 7) created generalized linear models look at the response variable vs. the sum of the response variable divided by the distance of the focal tree for each focal tree
      # a) Optional: looked for influential points (Points with Cook's D > 3 times the Cook's D) and remove them 
      # b) create different versions of generalized least squares regression (uncorrelated, exponential, 
            #spatial, spherical, linear, and ratio quadratics spatial correlations) and find the best 
            #predictive model by comparing the Akaike's Information Criterion
      # c) check the conditions (linearity, normal residuals, simple random sample). We used scatterplots to check linearity. 
            #We used a Shapiro test/qqnorm/histogram to look at normality of residuals. 
            #For a GLS, the errors are allowed to be correlated and/or unequal variances (heterodescadisticty)
      # d) use a normalized semi-variogram to check whether we controlled for spatial autocorrelation with our GLS model 
            #should hover around 1 if the model is effective
      # e) if conditions are met and the spatial autocorrelation is relatively controlled for, we can look at the model summary
            #and the slope test to see if there is significant impact of competition/facilitation on the growth of the trees.
            # A significant positive slope indicates facilitation.
            # A significant negative slope indicates competition.
      # f) use a non-parametric Kendall's Tau Test to see if there is a significant correlation between the competition metric
            #and the size of the focal trees for the data with no outliers


focal_function <- function(population){
  
  #assigning the dataframes based on the population
  if (population == "LM") {
    dataframe <- LM_fixed_field_data_processed
    dataframe_sf <- LM_fixed_field_data_processed_sf
  } else if (population == "LC") {
    dataframe <- LC_fixed_field_data_processed
    dataframe_sf <- LC_fixed_field_data_processed_sf
  } else if (population == "SD") {
    dataframe <- SD_fixed_field_data_processed
    dataframe_sf <- SD_fixed_field_data_processed_sf
  }
  
  
  
  # creating a bounding box based on the point locations
  box <- st_bbox(dataframe_sf)
  
  #cropping the tree points down by 20 m on all sides
  
  #creating a cropped bbox 
  box_sf <- box %>% #turning the bbox into polygon
    st_as_sfc()
  box_spatial <- as(box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
  box_spatial_cropped <- raster::crop(box_spatial, extent((box[1]+20), (box[3]-20), (box[2]+20), (box[4]-20))) #cropping the bounding box xmin, xmax, ymin, and ymax by 20 m inside
  box_sf_cropped <-  box_spatial_cropped %>% #turning the spatial polygon into a sfc polygon
    st_as_sfc()
  
  #cropping the points by the cropped box
  dataframe_cropped <- st_crop(dataframe_sf, box_sf_cropped)
  
  #Creating a grid over the cropped tree points 
  tree_grid_cropped <- st_make_grid(dataframe_cropped, cellsize = (((40*mean(dataframe$DBH_ag))*2)*2))
  
  #creating an x_sequential column that is 1 through the number of LC points
  dataframe_sf <- dataframe_sf %>%
    mutate(X_sequential = 1:nrow(dataframe_sf))
  
  
  #randomly selecting a focal point from each grid cell with trees within them
  list_grids_and_points <- st_contains(tree_grid_cropped, dataframe_sf, sparse =T) #find which points are within which grid cells, make sure row number in the data frame of grid cells corresponds to the order of the points dataframe within st_contains
  set.seed(25) #setting the seed
  list_grids_and_focal_trees <- lapply(list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
    if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
      focal_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
    }
    else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
      focal_pt <- cell #set the focal point to be the tree that is within the cell
    } else { # if there are no trees
      focal_pt <- NA # set the focal tree point to be NA
    }
    return(focal_pt)
  })
  
  #creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
  list_grids_and_focal_trees_df <- as.data.frame(unlist(list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
  colnames(list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
  list_grids_and_focal_trees_fixed <- list_grids_and_focal_trees_df %>% #filters out grid cells that do not have trees within them
    mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
    mutate(data_row = dataframe_sf$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
    filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them
  
  #filtering out point data to be just the focal points
  fixed_field_data_processed_focal <- dataframe_sf %>%
    filter(X_sequential %in% list_grids_and_focal_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 
  
  #creating the buffer around the focal points
  focal_tree_buffers <-st_buffer(fixed_field_data_processed_focal$geometry, 40*mean(fixed_field_data_processed_focal$DBH_ag))
  
  #calculating the size/distance for focal trees and neighbors within buffers for buffers with only the focal tree and with more 
  
  #create a tibble with the the number of trees within the buffers that contain trees
  tree_buffers_points_within_0 <- st_contains(focal_tree_buffers, dataframe_sf, sparse =F) %>%
    rowSums() %>% #find how many trees are within each grid
    as_tibble() %>% 
    mutate(row = row_number()) %>% #assign a new column with row numbers 
    filter(value > 0) #filter out any buffers with only the focal tree
  
  #filter out the buffers to only have the buffers that contain trees
  tree_buffer_inside_0 <- focal_tree_buffers %>%
    st_as_sf() %>% 
    mutate(row = row_number()) %>% #create a column with row numbers
    filter(row %in% tree_buffers_points_within_0$row) #only keep polygons that match the row number of the grid cells with trees within them 
  
  #Checking that row number in focal dataset is the same as the buffer dataset
  fixed_field_data_processed_focal_row <- fixed_field_data_processed_focal %>%
    as.data.frame() %>%
    mutate(row = as.factor(row_number())) %>%
    st_as_sf()
  tree_buffer_inside_0 <- mutate(tree_buffer_inside_0, row = as.factor(row)) #making sure the buffers have the same row number as the focal data
  
  #calculating the size/distance for focal trees and neighbors within buffers for buffers with more than just the focal tree
  
  #create a tibble with the the number of trees within the buffers that contain trees
  tree_buffers_points_within <- st_contains(focal_tree_buffers, dataframe_sf, sparse =F) %>%
    rowSums() %>% #find how many trees are within each grid
    as_tibble() %>%
    mutate(row = row_number()) %>% #assign a new column with row numbers
    filter(value > 1) #filter out any buffers with only the focal tree
  
  #filter out the buffers to only have the buffers that contain trees
  tree_buffer_inside <- focal_tree_buffers %>%
    st_as_sf() %>% 
    mutate(row = row_number()) %>% #create a column with row numbers
    filter(row %in% tree_buffers_points_within$row) #only keep buffers that match the row number of the buffers cells with trees within them 

  #Checking that row number in focal dataset is the same as the buffer dataset
  fixed_field_data_processed_focal_row <- fixed_field_data_processed_focal %>%
    as.data.frame() %>%
    mutate(row = as.factor(row_number())) %>%
    st_as_sf()
  tree_buffer_inside <- mutate(tree_buffer_inside, row = as.factor(row)) #making sure the buffers have the same row number as the isoalted focal data
  
  fixed_field_data_all_focal_trees <- tibble() #creating the empty tibble 
  
  #calculating the distances of each tree within the buffer to the focal tree and the competition metric values
  
  # in this loop, it iterates over each focal tree with neighbors and calculates the sum of the size metric 
  #divided by the distance for all of the neighbors for each focal tree (competition values for each tree)
  
  for (i in 1:nrow(fixed_field_data_processed_focal)){ #for the length focal trees 
    row_num = i
    tree_buffer_inside_df <- as.data.frame(tree_buffer_inside_0) #uses data of non-isolated and isolated focal trees 
    tree_buffer_inside_df_i <- tree_buffer_inside_df %>% 
      filter(row == row_num) #isolate a row of the buffer dataframe
    tree_buffer_inside_sf_i <- st_as_sf(tree_buffer_inside_df_i) #set the row as a simple feature
    all_pts_buffer <- st_contains(tree_buffer_inside_sf_i, fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
    possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
    fixed_field_data_processed_trees <- fixed_field_data_processed_sf %>%
      filter(X_sequential %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer.
    
    #correct sequence of focal trees
    correct_focal <- fixed_field_data_processed_focal[i,]$X_sequential
    
    #create a dataframe with only the focal tree 
    fixed_field_data_focal_tree <- fixed_field_data_processed_trees %>%
      filter(X_sequential %in% correct_focal$tree_row_num) 
    
    #isolating the neighbor tree data
    fixed_field_data_neighbor_trees <- fixed_field_data_processed_trees %>%
      filter(X_sequential %notin% fixed_field_data_focal_tree$X_sequential) #filtering out tree data for the neighbor trees 
    
    #if there are no neighbors, its sets the sum of the response variable divided by the distance of the tree to the focal tree to 0
    if(nrow(fixed_field_data_neighbor_trees) == 0){
      sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
      sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
      sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
      sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
      sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
    } else{
      
      # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
      fixed_field_data_neighbor_trees <-  fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
        mutate(focal_distance = as.numeric(st_distance(geometry, fixed_field_data_focal_tree$geometry, by_element = T))) %>% #calculate the distance between the focal tree and each tree that neighbors  by_element = T
        mutate(focal_distance = case_when(focal_distance == 0 ~ 0.0000016, 
                                          focal_distance != 0 ~ focal_distance)) %>% #replace values of 0 (if the coords are the same for multiple trees) with a value an order of magnitude smaller than the smallest distance in our dataset to avoid undefined values from division
        mutate(SCA_over_distance = Canopy_short/focal_distance) %>% #creating a column with the short canopy axis size value divided by the tree's distance from the focal tree
        mutate(LCA_over_distance = Canopy_long/focal_distance) %>%
        mutate(CA_over_distance = Canopy_area/focal_distance) %>%
        mutate(CS_over_distance = Crown_spread/focal_distance) %>%
        mutate(DBH_over_distance = DBH_ag/focal_distance)
      
      #create empty variables for the sum of the response variables over the distance of the trees to the focal trees
      sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
      sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
      sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
      sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
      sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
      
      #adding the size values of each neighbor to a sum total of the neighbors size values
      for (y in 1:nrow(fixed_field_data_neighbor_trees)){ 
        sum_SCA_over_distance = sum_SCA_over_distance + fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
        sum_LCA_over_distance = sum_LCA_over_distance + fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the LCA of each neighbor
        sum_CA_over_distance = sum_CA_over_distance + fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
        sum_CS_over_distance = sum_CS_over_distance + fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
        sum_DBH_over_distance = sum_DBH_over_distance + fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
      }
    }
    
    #creating a tibble with all of the calculated sizes over distances and other tree attributes for each focal tree
    all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
    fixed_field_data_focal_tree <- cbind(fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
    fixed_field_data_all_focal_trees <- rbind(fixed_field_data_all_focal_trees, fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
    
  }
  
  return(list(box_sf, box_sf_cropped, dataframe_cropped,
              tree_grid_cropped, focal_tree_buffers, fixed_field_data_processed_focal,
              tree_buffer_inside_0, fixed_field_data_processed_focal_row, tree_buffer_inside, 
              fixed_field_data_all_focal_trees))
}

#LM

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
focal_results <- focal_function("LM")

#assigning necessary dataframes and objects 
LM_box_sf <- focal_results[[1]] #bounding box
LM_box_sf_cropped <- focal_results[[2]] #bounding box cropped by 20m
LM_fixed_field_data_processed_sf_cropped <- focal_results[[3]] # cropped tree data
LM_tree_grid_cropped <- focal_results[[4]] #grid with 40*mean population DBH as grid size
LM_focal_tree_buffers <- focal_results[[5]] #focal tree buffers
LM_fixed_field_data_processed_focal <- focal_results[[6]] #focal tree data
LM_tree_buffer_inside_0 <- focal_results[[7]] #buffers that contain trees
LM_fixed_field_data_processed_focal_row <- focal_results[[8]] #focal trees
LM_tree_buffer_inside <- focal_results[[9]] #focal trees with trees inside buffer
LM_fixed_field_data_all_focal_trees <- focal_results[[10]] #focal tree dataframe

#plotting the original bounding box box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=LM_box_sf)+ #old box
  geom_sf(data=LM_box_sf_cropped)+ #cropped box
  geom_sf(data=LM_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=LM_fixed_field_data_processed_sf_cropped, color = "red") #old points

#graphing the selected focal trees, the buffers, the grid, colored by sequential ID number
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data=LM_focal_tree_buffers, color = "blue") +
  geom_sf(data= LM_fixed_field_data_processed_focal, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LM_tree_grid_cropped) +
  geom_sf(data=LM_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=LM_fixed_field_data_processed_focal_row, aes(color = row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = LM_focal_tree_buffers)+
  geom_sf(data = LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_fixed_field_data_processed_focal, color = 'blue')

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LM_tree_grid_cropped) +
  geom_sf(data=LM_tree_buffer_inside, aes(color = row))+
  geom_sf(data=LM_fixed_field_data_processed_focal_row, aes(color = row))


#descriptive statistics for the focal tree sum of size/shape metrics over distance

#histograms
ggplot(LM_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
LM_field_data_focal_summarized_focal <- LM_fixed_field_data_all_focal_trees %>%
  dplyr::select(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(LM_field_data_focal_summarized_focal)

#### creating the generalized linear effects model ####

#conditions are lINES: linearity, independence, normal distribution of residuals, equal variance, simple random sample

#creating x and y columns of the UTM 12N 
LM_fixed_field_data_all_focal_trees$X.1 <- st_coordinates(LM_fixed_field_data_all_focal_trees)[,1]
LM_fixed_field_data_all_focal_trees$Y <- st_coordinates(LM_fixed_field_data_all_focal_trees)[,2]

#SCA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
LM_lm_focal_SCA <- lm(Canopy_short ~ sum_SCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_SCA_cooks <- cooks.distance(LM_lm_focal_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential, meaning they change the slope of the linear model too much
LM_fixed_field_data_all_focal_trees_no_SCA_outliers <- LM_fixed_field_data_all_focal_trees[-c(21),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion, lowest AIC is the best predictive model
LM_AIC_test <- model.sel(LM_gls_focal_SCA, LM_gls_focal_SCA_exp, LM_gls_focal_SCA_gaus, LM_gls_focal_SCA_spher, LM_gls_focal_SCA_lin, LM_gls_focal_SCA_ratio)
LM_AIC_test

# While LM_gls_focal_SCA_lin has the lowest AIC, but we have had trouble with it, so we are using the second best option. 

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_SCA$residuals))+
  geom_qq()

#Shapiro test
shapiro.test(LM_gls_focal_SCA$residuals) #not significant so it is normal 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees, aes(x = LM_gls_focal_SCA$fitted, y = LM_gls_focal_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#checking we have appropriately removed the spatial autocorrelation
semivario <- Variogram( LM_gls_focal_SCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_SCA)

#non parametric Kendall's Tau Test for if the the residuals are not normal or if scatterplot seems non-linear 
LM_tau_result_SCA <- cor.test(LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance, 
                              LM_fixed_field_data_all_focal_trees$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_SCA)

# Calculate the trend line
LM_trend_line_SCA <- predict(loess(LM_fixed_field_data_all_focal_trees$Canopy_short ~ LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = (LM_fixed_field_data_all_focal_trees$Canopy_short), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = LM_trend_line_SCA), color = "red") +
  labs(x = "SCA over Distance", y = "Short Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#LCA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
LM_lm_focal_LCA <- lm(Canopy_long ~ sum_LCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_LCA_cooks <- cooks.distance(LM_lm_focal_LCA) #calculating the cook.s D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_LCA_cooks[(LM_lm_focal_LCA_cooks > (3 * mean(LM_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_all_focal_trees_no_LCA_outliers <- LM_fixed_field_data_all_focal_trees[-c(27),]


#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_LCA <- gls(Canopy_long ~ sum_LCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_exp <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_gaus <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_spher <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_lin <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_ratio <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_LCA <- model.sel(LM_gls_focal_LCA, LM_gls_focal_LCA_exp, LM_gls_focal_LCA_gaus, LM_gls_focal_LCA_spher, LM_gls_focal_LCA_lin, LM_gls_focal_LCA_ratio)
LM_AIC_test_LCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plots
ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_LCA$residuals))+
  geom_qq()

#shaprio wilk test, not significant so our residuals are normally distributed
shapiro.test(LM_gls_focal_LCA$residuals) 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_LCA$fitted, y = LM_gls_focal_LCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and LCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LM_gls_focal_LCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_LCA)

#non parametric Kendall's Tau Test for the version without outliers
LM_tau_result_LCA <- cor.test(LM_fixed_field_data_all_focal_trees$sum_LCA_over_distance, 
                              LM_fixed_field_data_all_focal_trees$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_LCA)

# Calculate the trend line
LM_trend_line_LCA <- predict(loess(LM_fixed_field_data_all_focal_trees$Canopy_long ~ LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LM_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LM_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#CA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_LA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
LM_lm_focal_CA <- lm(Canopy_area ~ sum_CA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_CA_cooks <- cooks.distance(LM_lm_focal_CA) #calculating the cook.s D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_CA_cooks[(LM_lm_focal_CA_cooks > (3 * mean(LM_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_all_focal_trees_no_CA_outliers <- LM_fixed_field_data_all_focal_trees[-c(27),]


#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_CA <- gls(Canopy_long ~ sum_CA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_exp <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_gaus <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_spher <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_lin <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_ratio <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CA <- model.sel(LM_gls_focal_CA, LM_gls_focal_CA_exp, LM_gls_focal_CA_gaus, LM_gls_focal_CA_spher, LM_gls_focal_CA_ratio) #LM_gls_focal_CA_lin
LM_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, not sign so the residuals are normally distributed
shapiro.test(LM_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_CA$fitted, y = LM_gls_focal_CA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LM_gls_focal_CA_spher, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_CA)

#non parametric Kendall's Tau Test for the version without outliers
LM_tau_result_CA <- cor.test(LM_fixed_field_data_all_focal_trees$sum_CA_over_distance, LM_fixed_field_data_all_focal_trees$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CA)

# Calculate the trend line
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_all_focal_trees$Canopy_area ~ LM_fixed_field_data_all_focal_trees$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = (LM_fixed_field_data_all_focal_trees$Canopy_area), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = LM_trend_line_CA), color = "red") +
  labs(x = "CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_CS_over_distance, y=Crown_spread)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CS over Distance")+
  ylab("Crown Spread")

#Cook's D
LM_lm_focal_CS <- lm(Crown_spread ~ sum_CS_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_CS_cooks <- cooks.distance(LM_lm_focal_CS) #calculating the cook.s D for each point
plot(LM_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_CS_cooks[(LM_lm_focal_CS_cooks > (3 * mean(LM_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_all_focal_trees_no_CS_outliers <- LM_fixed_field_data_all_focal_trees[-c(27),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_CS <- gls(Canopy_long ~ sum_CS_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_exp <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_gaus <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_spher <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_lin <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_ratio <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CS <- model.sel(LM_gls_focal_CS, LM_gls_focal_CS_lin, LM_gls_focal_CS_exp, LM_gls_focal_CS_gaus, LM_gls_focal_CS_spher, LM_gls_focal_CS_ratio) #without linear correlation
LM_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_CS$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_CS$residuals))+
  geom_qq()

# shapiro-wilk, not sign so residuals are normally distributed
shapiro.test(LM_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_CS$fitted, y = LM_gls_focal_CS$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LM_gls_focal_CS, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_CS)

#non parametric Kendall's Tau Test for the version without outliers
LM_tau_result_CS <- cor.test(LM_fixed_field_data_all_focal_trees$sum_Cs_over_distance, LM_fixed_field_data_all_focal_trees$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CS)

# Calculate the trend line
LM_trend_line_CS <- predict(loess(LM_fixed_field_data_all_focal_trees$Crown_spread ~ LM_fixed_field_data_all_focal_trees$sum_CS_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LM_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LM_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#plotting the linear model in ggplot
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_DBH_over_distance, y=DBH_ag)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("DBH over Distance")+
  ylab("DBH")

#Cook's D
LM_lm_focal_DBH <- lm(DBH_ag ~ sum_DBH_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_lm_focal_DBH_cooks <- cooks.distance(LM_lm_focal_DBH) #calculating the cook.s D for each point
plot(LM_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_DBH_cooks[(LM_lm_focal_DBH_cooks > (3 * mean(LM_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_fixed_field_data_all_focal_trees_no_CS_outliers <- LM_fixed_field_data_all_focal_trees[-c(27),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_DBH <- gls(Canopy_short ~ sum_DBH_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_DBH_exp <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_DBH_gaus <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_DBH_spher <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_DBH_lin <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_DBH_ratio <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_DHB <- model.sel(LM_gls_focal_DBH, LM_gls_focal_DBH_exp, LM_gls_focal_DBH_lin, LM_gls_focal_DBH_gaus, LM_gls_focal_DBH_spher, LM_gls_focal_DBH_ratio)
LM_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_DBH$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_DBH$residuals))+
  geom_qq()

# shapiro-wilk, not sign so the residuals are normally distributed
shapiro.test(LM_gls_focal_DBH$residuals) 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_DBH$fitted, y = LM_gls_focal_DBH$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and DBH over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LM_gls_focal_DBH, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_DBH)

#non parametric Kendall's Tau Test for the version without outliers
LM_tau_result_DBH <- cor.test(LM_fixed_field_data_all_focal_trees$sum_DBH_over_distance, LM_fixed_field_data_all_focal_trees$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_DBH)

# Calculate the trend line
LM_trend_line_DBH <- predict(loess(LM_fixed_field_data_all_focal_trees$DBH_ag ~ LM_fixed_field_data_all_focal_trees$sum_DBH_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_all_focal_trees$sum_DBH_over_distance, y = (LM_fixed_field_data_all_focal_trees$DBH_ag), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_all_focal_trees$sum_DBH_over_distance, y = LM_trend_line_DBH), color = "red") +
  labs(x = "DBH over Distance", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()


# LC

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
LC_focal_results <- focal_function("LC")

#assigning necessary dataframes and objects 
LC_box_sf <- LC_focal_results[[1]] #bounding box
LC_box_sf_cropped <- LC_focal_results[[2]] #bounding box cropped by 20m
LC_fixed_field_data_processed_sf_cropped <- LC_focal_results[[3]] # cropped tree data
LC_tree_grid_cropped <- LC_focal_results[[4]] #grid with 40*mean population DBH as grid size
LC_focal_tree_buffers <- LC_focal_results[[5]] #focal tree buffers
LC_fixed_field_data_processed_focal <- LC_focal_results[[6]] #focal tree data
LC_tree_buffer_inside_0 <- LC_focal_results[[7]] #buffers that contain trees
LC_fixed_field_data_processed_focal_row <- LC_focal_results[[8]] #focal trees
LC_tree_buffer_inside <- LC_focal_results[[9]] #focal trees with trees inside buffer
LC_fixed_field_data_all_focal_trees <- LC_focal_results[[10]] #focal tree dataframe

#plotting the original bounding box box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=LC_box_sf)+ #old box
  geom_sf(data=LC_box_sf_cropped)+ #cropped box
  geom_sf(data=LC_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=LC_fixed_field_data_processed_sf_cropped, color = "red") #old points

#graphing the selected focal trees, the buffers, the grid, colored by sequential ID number
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data=LC_focal_tree_buffers, color = "blue") +
  geom_sf(data= LC_fixed_field_data_processed_focal, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LC_tree_grid_cropped) +
  geom_sf(data=LC_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=LC_fixed_field_data_processed_focal_row, aes(color = row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = LC_focal_tree_buffers)+
  geom_sf(data = LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_fixed_field_data_processed_focal, color = 'blue')

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LC_tree_grid_cropped) +
  geom_sf(data=LC_tree_buffer_inside, aes(color = row))+
  geom_sf(data=lc_fixed_field_data_processed_focal_row, aes(color = row))

#plotting the tree points and their competition metrics
ggplot()+
  geom_sf(data=LC_fixed_field_data_all_focal_trees, aes(color = sum_SCA_over_distance))

#descriptive statistics

#histograms
ggplot(LC_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
LC_field_data_focal_summarized_focal <- LC_fixed_field_data_all_focal_trees %>%
  dplyr::select(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(LC_field_data_focal_summarized_focal)


#creating the generalized linear effects model

#creating x and y columns of the UTM 12N 
LC_fixed_field_data_all_focal_trees$X.1 <- st_coordinates(LC_fixed_field_data_all_focal_trees)[,1]
LC_fixed_field_data_all_focal_trees$Y <- st_coordinates(LC_fixed_field_data_all_focal_trees)[,2]

#SCA

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_all_focal_trees, (aes(x=sum_SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")


#Cook's D
LC_lm_focal_SCA <- lm(Canopy_short ~ sum_LCA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_lm_focal_SCA_cooks <- cooks.distance(LC_lm_focal_SCA) #calculating the cook.s D for each point
plot(LC_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_SCA_cooks[(LC_lm_focal_SCA_cooks > (3 * mean(LC_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_fixed_field_data_all_focal_trees_no_SCA_outliers <- LC_fixed_field_data_all_focal_trees[-c(3),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test <- model.sel(LC_gls_focal_SCA, LC_gls_focal_SCA_exp, LC_gls_focal_SCA_gaus, LC_gls_focal_SCA_spher, LC_gls_focal_SCA_lin, LC_gls_focal_SCA_ratio)
LC_AIC_test

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(x= LC_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_fixed_field_data_all_focal_trees, aes(sample = LC_gls_focal_SCA$residuals))+
  geom_qq()

#shapiro-wilk test, not sign so normal residuals
shapiro.test(LC_gls_focal_SCA$residuals) 

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees , aes(x = LC_gls_focal_SCA$fitted, y = LC_gls_focal_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram( LC_gls_focal_SCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_SCA)

#non parametric Kendall's Tau Test for the version without outliers
LC_tau_result_SCA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance, 
                              LC_fixed_field_data_all_focal_trees$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_SCA)

# Calculate the trend line
LC_trend_line_SCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_short ~ LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_short), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = LC_trend_line_SCA), color = "red") +
  labs(x = "SCA over Distance", y = "Short Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#LCA

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_all_focal_trees, (aes(x=sum_LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("LCA over Distance")+
  ylab("Long Canopy Axis")

#Cook's D
LC_lm_focal_LCA <- lm(Canopy_long ~ sum_LCA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_lm_focal_LCA_cooks <- cooks.distance(LC_lm_focal_LCA) #calculating the cook.s D for each point
plot(LC_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_LCA_cooks[(LC_lm_focal_LCA_cooks > (3 * mean(LC_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_fixed_field_data_all_focal_trees_no_LCA_outliers <- LC_fixed_field_data_all_focal_trees[-c(3),]

#unlogged version of generalized linear model
LC_gls_focal_LCA <- gls(Canopy_long ~ sum_LCA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_LCA_exp <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_LCA_gaus <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_LCA_spher <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_LCA_lin <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_LCA_ratio <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_LCA <- model.sel(LC_gls_focal_LCA, LC_gls_focal_LCA_exp, LC_gls_focal_LCA_gaus, LC_gls_focal_LCA_spher, LC_gls_focal_LCA_ratio) #LC_gls_focal_LCA_lin
LC_AIC_test_LCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(x= LC_gls_focal_LCA_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance",
       subtitle = "Using Gaussian Control for Spatial Autocorrelation")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_fixed_field_data_all_focal_trees, aes(sample = LC_gls_focal_LCA_gaus$residuals))+
  geom_qq()

#shapiro-wilk test, significant so non-normal residuals
shapiro.test(LC_gls_focal_LCA_gaus$residuals) 

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees, aes(x = LC_gls_focal_LCA_gaus$fitted, y = LC_gls_focal_LCA_gaus$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and LCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_LCA_gaus, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_LCA_gaus)

#because the residuals are not normal, we will use the Kendall's Tau correlation non-parametric test to see if the relationship is significant

#non parametric Kendall's Tau Test because the data is non-parametric and has ties 
LC_tau_result_LCA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, 
                              LC_fixed_field_data_all_focal_trees$Canopy_long,  
                              method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_LCA)

# Calculate the trend line
LC_trend_line_LCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_long ~ LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Extract fitted values from the GLS model
fitted_canopy <- fitted(LC_gls_focal_LCA_gaus)

# Create the data frame for plotting
line_df <- data.frame(
  sum_LCA_over_distance = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance,
  fitted_canopy = fitted_canopy
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(data = line_df, aes(x = sum_LCA_over_distance, y = fitted_canopy), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#CA

#plotting the linear model in ggplot
ggplot(data = LC_fixed_field_data_all_focal_trees, (aes(x=sum_CA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
LC_lm_focal_CA <- lm(Canopy_area ~ sum_CA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_lm_focal_CA_cooks <- cooks.distance(LC_lm_focal_LCA) #calculating the cook.s D for each point
plot(LC_lm_focal_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_CA_cooks[(LC_lm_focal_CA_cooks > (3 * mean(LC_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_fixed_field_data_all_focal_trees_no_CA_outliers <- LC_fixed_field_data_all_focal_trees[-c(3),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_CA <- gls(Canopy_area ~ sum_CA_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CA_exp <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CA_gaus <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CA_spher <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CA_lin <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CA_ratio <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CA <- model.sel(LC_gls_focal_CA, LC_gls_focal_CA_exp, LC_gls_focal_CA_gaus, LC_gls_focal_CA_spher, LC_gls_focal_CA_lin, LC_gls_focal_CA_ratio)
LC_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(x= LC_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_fixed_field_data_all_focal_trees, aes(sample = LC_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, sign for when residuals so we are using a Kendall's Tau Correlation Test non-parametric data
shapiro.test(LC_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees , aes(x = LC_gls_focal_CA$fitted, y = LC_gls_focal_CA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_CA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_CA)

#non parametric Kendall's Tau Test Test for the version without outliers
LC_tau_result_CA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, 
                             LC_fixed_field_data_all_focal_trees$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CA)

# Calculate the trend line
LC_trend_line_CA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_area ~ LC_fixed_field_data_all_focal_trees$sum_CA_over_distance))

# Create a trend line plot

# Extract fitted values from the GLS model
fitted_canopy <- fitted(LC_gls_focal_CA)

# Create the data frame for plotting
line_df <- data.frame(
  sum_CA_over_distance = LC_fixed_field_data_all_focal_trees$sum_CA_over_distance,
  fitted_canopy = fitted_canopy
)

#plotting
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_area), color = "red")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = fitted_canopy), color = "red") +
  labs(x = "CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot for SCA
ggplot(data = LC_fixed_field_data_all_focal_trees, (aes(x=sum_CS_over_distance, y=Crown_spread)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CS over Distance")+
  ylab("Crown Spread")

#Cook's D
LC_lm_focal_CS <- lm(Crown_spread ~ sum_CS_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_lm_focal_CS_cooks <- cooks.distance(LC_lm_focal_CS) #calculating the cook.s D for each point
plot(LC_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_CS_cooks[(LC_lm_focal_CS_cooks > (3 * mean(LC_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_fixed_field_data_all_focal_trees_no_CS_outliers <- LC_fixed_field_data_all_focal_trees[-c(3),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_CS <- gls(Crown_spread ~ sum_CS_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CS_exp <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CS_gaus <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CS_spher <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CS_lin <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_CS_ratio <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CS <- model.sel(LC_gls_focal_CS, LC_gls_focal_CS_exp, LC_gls_focal_CS_gaus, LC_gls_focal_CS_ratio) #without linear correlation and without spherical
LC_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(x= LC_gls_focal_CS_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq norm
ggplot(LC_fixed_field_data_all_focal_trees, aes(sample = LC_gls_focal_CS_gaus$residuals))+
  geom_qq()

# shapiro-wilk, n sign for both versions with and without outliers so used Kendall's Tau non-parametric test
shapiro.test(LC_gls_focal_CA$residuals) # shapiro-wilk, n sign for both versions with and without outliers so used Kendall's Tau Test non-parametric test

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees , aes(x = LC_gls_focal_CS_gaus$fitted, y = LC_gls_focal_CS_gaus$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_CS_gaus, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_CS_gaus)

#non parametric Kendall's Tau Test Test for the version without outliers
LC_tau_result_CS <- cor.test(LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, 
                             LC_fixed_field_data_all_focal_trees$Crown_spread,  
                             method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CS)

# Calculate the trend line
LC_trend_line_CS <- predict(loess(LC_fixed_field_data_all_focal_trees$Crown_spread ~ LC_fixed_field_data_all_focal_trees$sum_CS_over_distance))

# Extract fitted values from the GLS model
fitted_crown <- fitted(LC_gls_focal_CS_gaus)

# Create the data frame for plotting
line_df <- data.frame(
  sum_CA_over_distance = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance,
  fitted_crown = fitted_crown
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "red")) +
  geom_line(data = line_df, aes(x = sum_CS_over_distance, y = fitted_crown), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#plotting the linear model in ggplot 
ggplot(data = LC_fixed_field_data_all_focal_trees, (aes(x=sum_DBH_over_distance, y=DBH_ag)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("DBH over Distance")+
  ylab("DBH")

#Cook's D
LC_lm_focal_DBH <- lm(DBH_ag ~ sum_DBH_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_lm_focal_DBH_cooks <- cooks.distance(LC_lm_focal_DBH) #calculating the cook.s D for each point
plot(LC_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_DBH_cooks[(LC_lm_focal_DBH_cooks > (3 * mean(LC_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_fixed_field_data_all_focal_trees_no_DBH_outliers <- LC_fixed_field_data_all_focal_trees[-c(16, 19),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_DBH <- gls(DBH_ag ~ sum_DBH_over_distance, data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_DBH_exp <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_DBH_gaus <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_DBH_spher <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_DBH_lin <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_DBH_ratio <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_DHB <- model.sel(LC_gls_focal_DBH, LC_gls_focal_DBH_exp, LC_gls_focal_DBH_gaus, LC_gls_focal_DBH_spher, LC_gls_focal_DBH_ratio) #without linear correlation
LC_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(x= LC_gls_focal_DBH_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plot
ggplot(LC_fixed_field_data_all_focal_trees, aes(sample = LC_gls_focal_DBH_gaus$residuals))+
  geom_qq()

# shapiro-wilk, not significant so normal residuals
shapiro.test(LC_gls_focal_DBH_gaus$residuals)

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees, aes(x = LC_gls_focal_DBH$fitted, y = LC_gls_focal_DBH$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and DBH over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_DBH_gaus, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_DBH_gaus)

#non parametric Kendall's Tau Test Test for the version without outliers
LC_tau_result_DBH <- cor.test(LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, 
                              LC_fixed_field_data_all_focal_trees$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_DBH)

# Calculate the trend line
LC_trend_line_DBH <- predict(loess(LC_fixed_field_data_all_focal_trees$DBH_ag ~ LC_fixed_field_data_all_focal_trees$sum_DBH_over_distance))

# Extract fitted values from the GLS model
fitted_DBH <- fitted(LC_gls_focal_DBH_gaus)

# Create the data frame for plotting
LC_line_df_DBH <- data.frame(
  sum_DBH_over_distance = LC_fixed_field_data_all_focal_trees$sum_DBH_over_distance,
  fitted_DBH = fitted_DBH
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$DBH_ag), color = "red")) +
  geom_line(data = LC_line_df_DBH, aes(x = sum_DBH_over_distance, y = fitted_DBH), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


# SD

focal_function <- function(population){
  
  #assigning the dataframes based on the population
  if (population == "SD") {
    dataframe <- SD_fixed_field_data_processed
    dataframe_sf <- SD_fixed_field_data_processed_sf
  }
  
  # creating a bounding box based on the point locations
  box <- st_bbox(dataframe_sf)
  
  #cropping the tree points down by 20 m on all sides
  
  #creating a cropped bbox 
  box_sf <- box %>% #turning the bbox into polygon
    st_as_sfc()
  box_spatial <- as(box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
  box_spatial_cropped <- raster::crop(box_spatial, extent((box[1]+20), (box[3]-20), (box[2]+20), (box[4]-20))) #cropping the xmin, xmax, ymin, and ymax by 20 m inside
  box_sf_cropped <-  box_spatial_cropped %>% #turning the spatial polygon into a polygon
    st_as_sfc()
  
  #cropping the points by the cropped box
  fixed_field_data_processed_sf_cropped<- st_crop(dataframe_sf, box_sf_cropped)
  
  #Creating a grid over the cropped tree points 
  tree_grid_cropped <- st_make_grid(fixed_field_data_processed_sf_cropped, cellsize = (((40*mean(dataframe$DBH_ag))*2)*2))
  
  #creating an x_sequential column that is 1 through the number of SD points
  dataframe_sf <- dataframe_sf %>%
    mutate(X_sequential = 1:nrow(dataframe_sf))
  
  #selecting a focal point from each grid cell with trees within them
  list_grids_and_points <- st_contains(tree_grid_cropped, dataframe_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
  set.seed(25) #setting the seed
  list_grids_and_focal_trees <- lapply(list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
    if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
      focal_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
    }
    else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
      focal_pt <- cell #set the focal point to be the tree that is within the cell
    } else { # if there are no trees
      focal_pt <- NA # set the focal tree point to be NA
    }
    return(focal_pt)
  })
  
  #creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
  list_grids_and_focal_trees_df <- as.data.frame(unlist(list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
  colnames(list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
  list_grids_and_focal_trees_fixed <- list_grids_and_focal_trees_df %>% #filters out grid cells that do not have trees within them
    mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
    mutate(data_row = dataframe_sf$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
    filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them
  
  #filtering out point data to be just the focal points
  fixed_field_data_processed_focal <- dataframe_sf %>%
    filter(X_sequential %in% list_grids_and_focal_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 
  
  #creating the buffer around the focal points
  focal_tree_buffers <-st_buffer(fixed_field_data_processed_focal$geometry, 40*mean(fixed_field_data_processed_focal$DBH_ag))

  #calculating the size/distance for focal trees and neighbors within buffers for buffers with only the focal tree and with more 
  
  #create a tibble with the the number of trees within the buffers that contain trees
  tree_buffers_points_within_0 <- st_contains(focal_tree_buffers, dataframe_sf, sparse =F) %>%
    rowSums() %>% #find how many trees are within each grid
    as_tibble() %>% 
    mutate(row = row_number()) %>% #assign a new column with row numbers 
    filter(value > 0) #filter out any buffers with only the focal tree
  
  #filter out the buffers to only have the buffers that contain trees
  tree_buffer_inside_0 <- focal_tree_buffers %>%
    st_as_sf() %>% 
    mutate(row = row_number()) %>% #create a column with row numbers
    filter(row %in% tree_buffers_points_within_0$row) #only keep polygons that match the row number of the grid cells with trees within them 
  
  #Checking that row number in focal dataset is the same as the buffer dataset
  fixed_field_data_processed_focal_row <- fixed_field_data_processed_focal %>%
    as.data.frame() %>%
    mutate(row = as.factor(row_number())) %>%
    st_as_sf()
  tree_buffer_inside_0 <- mutate(tree_buffer_inside_0, row = as.factor(row)) #making sure the buffers have the same row number as the focal data
  

  #calculating the size/distance for focal trees and neighbors within buffers for buffers with more than just the focal tree
  
  #create a tibble with the the number of trees within the buffers that contain trees
  tree_buffers_points_within <- st_contains(focal_tree_buffers, dataframe_sf, sparse =F) %>%
    rowSums() %>% #find how many trees are within each grid
    as_tibble() %>% 
    mutate(row = row_number()) %>% #assign a new column with row numbers 
    filter(value > 1) #filter out any buffers with only the focal tree
  
  #filter out the buffers to only have the buffers that contain trees
  tree_buffer_inside <- focal_tree_buffers %>%
    st_as_sf() %>% 
    mutate(row = row_number()) %>% #create a column with row numbers
    filter(row %in% tree_buffers_points_within$row) #only keep buffers that match the row number of the buffers cells with trees within them 
  
  #Checking that row number in focal dataset is the same as the buffer dataset
  fixed_field_data_processed_focal_row <- fixed_field_data_processed_focal %>%
    as.data.frame() %>%
    mutate(row = as.factor(row_number())) %>%
    st_as_sf()
  tree_buffer_inside <- mutate(tree_buffer_inside, row = as.factor(row)) #making sure the buffers have the same row number as the isolated focal data
  
  #creating the empty tibble 
  fixed_field_data_all_focal_trees <- tibble()
  
  #calculating the distances of each tree within the buffer to the focal tree and the competition metric values
  for (i in 1:nrow(fixed_field_data_processed_focal)){ #for the length of the buffers with trees inside of them #fixed_field_data_processed_focal_row
    row_num = i
    tree_buffer_inside_df <- as.data.frame(tree_buffer_inside_0) #uses data of non-isolated and isolated focal trees 
    tree_buffer_inside_df_i <- tree_buffer_inside_df %>% 
      filter(row == row_num) #isolate a row of the buffer dataframe
    tree_buffer_inside_sf_i <- st_as_sf(tree_buffer_inside_df_i) #set the row as a simple feature
    all_pts_buffer <- st_contains(tree_buffer_inside_sf_i, dataframe_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
    possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
    fixed_field_data_processed_trees <- dataframe_sf %>%
      filter(X_sequential %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer. #before it was X %in% possible_pts_buffer
    
    #create column with correct tree number
    correct_focal <- fixed_field_data_processed_focal[i,]$X_sequential
    
    #create a dataframe with only the focal tree
    fixed_field_data_focal_tree <- fixed_field_data_processed_trees %>%
      filter(X_sequential %in% correct_focal) 
    
    #isolating the neighbor tree data
    fixed_field_data_neighbor_trees <- fixed_field_data_processed_trees %>%
      filter(X_sequential %notin% fixed_field_data_focal_tree$X_sequential) #filtering out tree data for the neighbor trees 
    
    #if there are no neighbors, its sets the sum of the response variable divided by the distance of the tree to the focal tree to 0
    if(nrow(fixed_field_data_neighbor_trees) == 0){
      sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
      sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
      sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
      sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
      sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
    } else{
      
      # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
      fixed_field_data_neighbor_trees <-  fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
        mutate(focal_distance = as.numeric(st_distance(geometry,  fixed_field_data_focal_tree$geometry))) %>% #caSDulate the distance between the focal tree and each tree that neighbors it
        mutate(focal_distance = case_when(focal_distance == 0 ~ 0.0000016, 
                                          focal_distance != 0 ~ focal_distance)) %>% #replace values of 0 (if the coords are the same for multiple trees) with a value an order of magnitude smaller than the smallest distance in our dataset
        mutate(SCA_over_distance = Canopy_short/focal_distance) %>% #creating a column with the short canopy axis size value divided by the tree's distance from the focal tree
        mutate(LCA_over_distance = Canopy_long/focal_distance) %>%
        mutate(CA_over_distance = Canopy_area/focal_distance) %>%
        mutate(CS_over_distance = Crown_spread/focal_distance) %>%
        mutate(DBH_over_distance = DBH_ag/focal_distance)
      
      #create empty variables for the sum of the response variables over the distance of the trees to the focal trees
      sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
      sum_SDA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
      sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
      sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
      sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
      
      #adding the size values of each neighbor to a sum total of the neighbors size values
      for (y in 1:nrow(fixed_field_data_neighbor_trees)){ #adding the size values of each neighbor to a sum total of the neighbors size values
        sum_SCA_over_distance = sum_SCA_over_distance + fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
        sum_LCA_over_distance = sum_LCA_over_distance + fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the SDA of each neighbor
        sum_CA_over_distance = sum_CA_over_distance + fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
        sum_CS_over_distance = sum_CS_over_distance + fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
        sum_DBH_over_distance = sum_DBH_over_distance + fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
      }
    }
    #creating a tibble with all of the calculated sizes over distances and other tree attributes for each focal tree
    all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
    fixed_field_data_focal_tree <- cbind(fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
    fixed_field_data_all_focal_trees <- rbind(fixed_field_data_all_focal_trees, fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
    
    
  }
  
  return(list(box_sf, box_sf_cropped, dataframe_cropped,
              tree_grid_cropped, focal_tree_buffers, fixed_field_data_processed_focal,
              tree_buffer_inside_0, fixed_field_data_processed_focal_row, tree_buffer_inside, 
              fixed_field_data_all_focal_trees))
  
}

#creating the empty tibble 
fixed_field_data_all_focal_trees <- tibble()

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
SD_focal_results <- focal_function("SD")

#assigning necessary dataframes and objects 
SD_box_sf <- SD_focal_results[[1]] #bounding box
SD_box_sf_cropped <- SD_focal_results[[2]] #bounding box cropped by 20m
SD_fixed_field_data_processed_sf_cropped <- SD_focal_results[[3]] # cropped tree data
SD_tree_grid_cropped <- SD_focal_results[[4]] #grid with 40*mean population DBH as grid size
SD_focal_tree_buffers <- SD_focal_results[[5]] #focal tree buffers
SD_fixed_field_data_processed_focal <- SD_focal_results[[6]] #focal tree data
SD_tree_buffer_inside_0 <- SD_focal_results[[7]] #buffers that contain trees
SD_fixed_field_data_processed_focal_row <- SD_focal_results[[8]] #focal trees
SD_tree_buffer_inside <- SD_focal_results[[9]] #focal trees with trees inside buffer
SD_fixed_field_data_all_focal_trees <- SD_focal_results[[10]] #focal tree dataframe

#plotting the original bounding box box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=SD_box_sf)+ #old box
  geom_sf(data=SD_box_sf_cropped)+ #cropped box
  geom_sf(data=SD_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=SD_fixed_field_data_processed_sf_cropped, color = "red") #old points

#graphing the selected focal trees, the buffers, the grid, colored by sequential ID number
ggplot()+
  geom_sf(data = SD_tree_grid_cropped)+
  geom_sf(data=SD_focal_tree_buffers, color = "blue") +
  geom_sf(data= SD_fixed_field_data_processed_focal, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data=SD_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=SD_fixed_field_data_processed_focal_row, aes(color = row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = SD_focal_tree_buffers)+
  geom_sf(data = SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_fixed_field_data_processed_focal, color = 'blue')

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data=SD_tree_buffer_inside, aes(color = row))+
  geom_sf(data=SD_fixed_field_data_processed_focal_row, aes(color = row))

#plotting the tree points and their competition metrics
ggplot()+
  geom_sf(data=SD_fixed_field_data_all_focal_trees, aes(color = sum_SCA_over_distance))


#descriptive statistics

#histograms
ggplot(SD_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(SD_fixed_field_data_all_focal_trees) + # Generate the base plot
  geom_histogram(aes(x = sum_DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
SD_field_data_focal_summarized_focal <- SD_fixed_field_data_all_focal_trees %>%
  dplyr::select(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(SD_field_data_focal_summarized_focal)


#creating x and y columns of the UTM 12N 
SD_fixed_field_data_all_focal_trees$X.1 <- st_coordinates(SD_fixed_field_data_all_focal_trees)[,1]
SD_fixed_field_data_all_focal_trees$Y <- st_coordinates(SD_fixed_field_data_all_focal_trees)[,2]

#SCA

#plotting the linear model in ggplot 
ggplot(data = SD_fixed_field_data_all_focal_trees, (aes(x=sum_SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
SD_lm_focal_SCA <- lm(Canopy_short ~ sum_SCA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_SCA_cooks <- cooks.distance(SD_lm_focal_SCA) #calculating the cook.s D for each point
plot(SD_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_SCA_cooks[(SD_lm_focal_SCA_cooks > (3 * mean(SD_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_SCA_outliers <- SD_fixed_field_data_all_focal_trees[-c(3, 24),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_SCA <- model.sel(SD_gls_focal_SCA, SD_gls_focal_SCA_lin, SD_gls_focal_SCA_exp, SD_gls_focal_SCA_gaus, SD_gls_focal_SCA_spher, SD_gls_focal_SCA_ratio) 
SD_AIC_test_SCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq nrom
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_SCA$residuals))+
  geom_qq()

#shapiro-welk test, not significant so normal residuals
shapiro.test(SD_gls_focal_SCA$residuals)

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees , aes(x = SD_gls_focal_SCA$fitted, y = SD_gls_focal_SCA_gaus$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_SCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_SCA)

#non parametric Kendall's Tau Test
SD_tau_result_SCA <- cor.test(SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance, 
                              SD_fixed_field_data_all_focal_trees_no_SCA_outliers$Canopy_short,  
                              method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_SCA)

# Calculate the trend line
SD_trend_line_SCA <- predict(loess(SD_fixed_field_data_all_focal_trees$Canopy_short ~ SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

# Extract fitted values from the GLS model
fitted_SCA <- fitted(SD_gls_focal_SCA)

# Create the data frame for plotting
SD_line_df_SCA <- data.frame(
  sum_SCA_over_distance = SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance,
  fitted_SCA = fitted_SCA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = (SD_fixed_field_data_all_focal_trees$Canopy_short), color = "blue")) +
  geom_line(data = SD_line_df_SCA, aes(x = sum_SCA_over_distance, y = fitted_SCA), color = "red") +
  labs(x = "Sum of SCA over Distance", y = "Short Canopy Axis ", title = "Trend Line Plot") +
  theme_minimal()


#LCA

#plotting the linear model in ggplot for LCA
ggplot(data = SD_fixed_field_data_all_focal_trees, (aes(x=sum_LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Long Canopy Axis")

#Cook's D
SD_lm_focal_LCA <- lm(Canopy_long ~ sum_LCA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_LCA_cooks <- cooks.distance(SD_lm_focal_LCA) #calculating the cook.s D for each point
plot(SD_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_LCA_cooks[(SD_lm_focal_LCA_cooks > (3 * mean(SD_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_LCA_outliers <- SD_fixed_field_data_all_focal_trees[-c(14, 23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_LCA <- gls(Canopy_long ~ sum_LCA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_LCA_exp <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_LCA_gaus <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_LCA_spher <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_LCA_lin <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_LCA_ratio <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_LCA <- model.sel(SD_gls_focal_LCA, SD_gls_focal_LCA_exp, SD_gls_focal_LCA_gaus, SD_gls_focal_LCA_spher, SD_gls_focal_LCA_ratio) #without the linear correction model
SD_AIC_test_LCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_LCA$residuals))+
  geom_qq()

#shapiro-wilk test, not sign so normal residuals
shapiro.test(SD_gls_focal_LCA$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees , aes(x = SD_gls_focal_LCA$fitted, y = SD_gls_focal_LCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and LCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_LCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_LCA)

#non parametric Kendall's Tau Test
SD_tau_result_LCA <- cor.test(SD_fixed_field_data_all_focal_trees_no_LCA_outliers$sum_LCA_over_distance, SD_fixed_field_data_all_focal_trees_no_LCA_outliers$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_LCA)

# Calculate the trend line
SD_trend_line_LCA <- predict(loess(SD_fixed_field_data_all_focal_trees$Canopy_long ~ SD_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Extract fitted values from the GLS model
fitted_LCA <- fitted(SD_gls_focal_LCA)

# Create the data frame for plotting
SD_line_df_LCA <- data.frame(
  sum_LCA_over_distance = SD_fixed_field_data_all_focal_trees$sum_LCA_over_distance,
  fitted_LCA = fitted_LCA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (SD_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(data = SD_line_df_LCA, aes(x = sum_LCA_over_distance, y = fitted_LCA), color = "red") +
  labs(x = "Sum of LCA over Distance", y = "Long Canopy Axis ", title = "Trend Line Plot") +
  theme_minimal()

#CA

#plotting the linear model in ggplot 
ggplot(data = SD_fixed_field_data_all_focal_trees, (aes(x=sum_CA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
SD_lm_focal_CA <- lm(Canopy_area ~ sum_CA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_CA_cooks <- cooks.distance(SD_lm_focal_CA) #calculating the cook.s D for each point
plot(SD_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_CA_cooks[(SD_lm_focal_CA_cooks > (3 * mean(SD_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_CA_outliers <- SD_fixed_field_data_all_focal_trees[-c(22),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CA <- gls(Canopy_area ~ sum_CA_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CA_exp <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CA_gaus <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CA_spher <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CA_lin <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CA_ratio <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CA <- model.sel(SD_gls_focal_CA, SD_gls_focal_CA_exp, SD_gls_focal_CA_gaus, SD_gls_focal_CA_lin, SD_gls_focal_CA_ratio) #SD_gls_focal_CA_spher
SD_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, significant so residuals non-normal
shapiro.test(SD_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees, aes(x = SD_gls_focal_CA$fitted, y = SD_gls_focal_CA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_CA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_AIC_test_CA)

#non parametric Kendall's Tau Test
SD_tau_result_CA <- cor.test(SD_fixed_field_data_all_focal_trees$sum_CA_over_distance, 
                             SD_fixed_field_data_all_focal_trees$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CA)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_all_focal_trees$Canopy_area ~ SD_fixed_field_data_all_focal_trees$sum_CA_over_distance))

# Extract fitted values from the GLS model
fitted_CA <- fitted(SD_gls_focal_CA)

# Create the data frame for plotting
SD_line_df_CA <- data.frame(
  sum_CA_over_distance = SD_fixed_field_data_all_focal_trees$sum_CA_over_distance,
  fitted_CA = fitted_CA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = (SD_fixed_field_data_all_focal_trees$Canopy_area), color = "red")) +
  geom_line(data = SD_line_df_CA, aes(x = sum_CA_over_distance, y = fitted_CA), color = "red") +
  labs(x = "Sum of CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot 
ggplot(data = SD_fixed_field_data_all_focal_trees, (aes(x=sum_CS_over_distance, y=Crown_spread)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("Crown Spread Competition Metric")+ #CS over Distance
  ylab("Crown Spread (m)")+
  theme_classic()+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))


#Cook's D
SD_lm_focal_CS <- lm(Crown_spread ~ sum_CS_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_CS_cooks <- cooks.distance(SD_lm_focal_CS) #calculating the cook.s D for each point
plot(SD_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_CS_cooks[(SD_lm_focal_CS_cooks > (3 * mean(SD_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_CS_outliers <- SD_fixed_field_data_all_focal_trees[-c(22, 28),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CS <- gls(Crown_spread ~ sum_CS_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_exp <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_gaus <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_spher <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_lin <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_ratio <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CS <- model.sel(SD_gls_focal_CS, SD_gls_focal_CS_exp, SD_gls_focal_CS_gaus, SD_gls_focal_CS_spher, SD_gls_focal_CS_ratio) #without linear correlation
SD_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_CS$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_CS$residuals))+
  geom_qq()

# shapiro-wilk, not signficant, meaning not significantly different from normal
shapiro.test(SD_gls_focal_CS$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees , aes(x = SD_gls_focal_CS$fitted, y = SD_gls_focal_CS$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and CS over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_CS, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_CS)

#non parametric Kendall's Tau Test
SD_tau_result_CS <- cor.test(SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, SD_fixed_field_data_all_focal_trees$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CS)

# Calculate the trend line
SD_trend_line_CS <- predict(loess(SD_fixed_field_data_all_focal_trees$Crown_spread ~ SD_fixed_field_data_all_focal_trees$sum_CS_over_distance))

# Extract fitted values from the GLS model
fitted_CS <- fitted(SD_gls_focal_CS)

# Create the data frame for plotting
SD_line_df_CS <- data.frame(
  sum_CS_over_distance = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance,
  fitted_CS = fitted_CS
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (SD_fixed_field_data_all_focal_trees$Crown_spread), color = "red")) +
  geom_line(data = SD_line_df_CS, aes(x = sum_CS_over_distance, y = fitted_CS), color = "red") +
  labs(x = "Sum of CS over Distance", y = "Crown Spread", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#plotting the linear model in ggplot 
ggplot(data = SD_fixed_field_data_all_focal_trees, (aes(x=sum_DBH_over_distance, y=DBH_ag)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("DBH Competition Metric")+
  ylab("DBH (cm)")+
  ylim(c(0,1))+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))


#Cook's D
SD_lm_focal_DBH <- lm(DBH_ag ~ sum_DBH_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_DBH_cooks <- cooks.distance(SD_lm_focal_DBH) #calculating the cook.s D for each point
plot(SD_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_DBH_cooks[(SD_lm_focal_DBH_cooks > (3 * mean(SD_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_DBH_outliers <- SD_fixed_field_data_all_focal_trees[-c(22),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_DBH <- gls(DBH_ag ~ sum_DBH_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_exp <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_gaus <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_spher <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_lin <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_ratio <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_DHB <- model.sel(SD_gls_focal_DBH, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_exp, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_gaus, SD_gls_focal_DBH_spher, SD_gls_focal_DBH_ratio) 
SD_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_DBH$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_DBH$residuals))+
  geom_qq()

# shapiro-wilk, not significant so normal
shapiro.test(SD_gls_focal_DBH$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees , aes(x = SD_gls_focal_DBH$fitted, y = SD_gls_focal_DBH$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and DBH over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_DBH, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_DBH)

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CDBH)

# Calculate the trend line
SD_trend_line_DBH <- predict(loess(SD_fixed_field_data_all_focal_trees$DBH_ag ~ SD_fixed_field_data_all_focal_trees$sum_DBH_over_distance))

# Extract fitted values from the GLS model
fitted_DBH <- fitted(SD_gls_focal_DBH)

# Create the data frame for plotting
SD_line_df_DBH <- data.frame(
  sum_DBH_over_distance = SD_fixed_field_data_all_focal_trees$sum_DBH_over_distance,
  fitted_DBH = fitted_DBH
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_DBH_over_distance, y = (SD_fixed_field_data_all_focal_trees$DBH_ag), color = "red")) +
  geom_line(data = SD_line_df_DBH, aes(x = sum_DBH_over_distance, y = fitted_DBH), color = "red") +
  labs(x = "Sum of DBH over Distance", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

