# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei compete or facilitate with one another%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%with permutations of the slope test%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by the distance to other individuals of the same species 
# either due to competition or facilitation. 
# If they are impacted by facilitation, we would expect closer trees would be bigger. 
# If they are impacted by competition, we would expect closer trees to be smaller. 
# To test this, we used also created generalized linear models to see if for focal trees, there was a 
# relationship between how much competition the trees face (based on the 
# size of the neighbors over their distance to the focal trees) and the size of the focal trees.

# We create a focal_function() to generate focal trees and calculate the competition metrics
# for each tree and for each population. We then created a slope_tests() function to run the focal_function() 500
# times with different seeds to randomly generate a new set of focal trees, find the best fitting GLS model, and 
# then run the slope/Kendall's Tau to support the reliability of our findings.  

# The script is broken into sections of 
# 1) loading and processing the packages and spatial/size/shape data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations and loading in the river outline shapefiles, 
# 2) Creating the focal_function() to generate focal trees and calculate the competition metrics
#for each tree and for each population
# 3) Creating the slope_tests() function to run the focal_function() 500
#times with different seeds to randomly generate a new set of focal trees, find the best fitting GLS model, and 
#then run the slope/Kendall's Tau to support the reliability of our findings. 
# 2) Using the slope_tests() function to generate the generalized linear regressions to see if tree size seem related to local competition.
# We also store slope test (parametric) and Kendall's Tau (non-parametric) results with and without Bonferroni Corrections and generate histograms 
# and descriptive summary statistics to view the distribution of p-values to check how robust the findings are. 

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

# loading in the processed tree data 
source("./analyses/Data_Processing_Script.R")
# 
# # loading in the tree data (size, elevation, lat/lon, ID, size/shape)
# 
# fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R
# 
# # creating the point shapefiles of the tree locations for each population in UTM 12 N
# 
# #creating a point shapefile of all points with lat lon coordinates and other attributes in WGS 1984
# #sf objects are dataframes with rows representing simple features with attributes and a simple feature geometry list-column (sfc)
# fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
#                                           coords = c("long", "lat"), crs = 4326)
# 
# #creating a transformed point shapefile with UTM 12 N an equal area projection
# fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) 
# 
# #storing point shapefiles for the trees by population
# 
# LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LM") %>%
#   st_as_sf()
# 
# LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "LC") %>%
#   st_as_sf()
# 
# SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
#   filter(Locality == "SD") %>%
#   st_as_sf()
# 
# #create dataframe with X and Y UTM coordinates
# 
# fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with separate x and y columns from the UTM 12N transformation
# fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
#   cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe
# 
# # creating a dataframe with the 5 average nearest neighbors (ANN) for each individual tree/row
# fixed_field_data_processed_NN_UTM <- fixed_field_data_processed_sf_trans_coordinates %>%  #creates a dataframe with the ANN of the closest 5 individual trees for each individual
#   mutate(dist1 = nndist(X = X.1, Y= Y, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
#   mutate(dist2 = nndist(X = X.1, Y= Y, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
#   mutate(dist3 = nndist(X = X.1, Y= Y, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
#   mutate(dist4 = nndist(X = X.1, Y= Y, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
#   mutate(dist5 = nndist(X = X.1, Y= Y, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
#   rowwise()%>% #so that in the next part we take the averages across rows
#   mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5))) # %>% #creates a column of the average distances (1-5) of each individual
# 
# # Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns
# 
# LM_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
#   filter(Locality == "LM")
# 
# LC_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
#   filter(Locality == "LC")
# 
# SD_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
#   filter(Locality == "SD")



#### Creating the Generalized Linear Model Functions ####

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


focal_function <- function(population, seed_input){
  
  
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
  tree_grid_cropped <- st_make_grid(dataframe_cropped, cellsize = (((40*mean(dataframe$DBH_ag))*2)*2)) #ASH NOTE: why is this *2*2 and not *4???
  
  #creating an x_sequential column that is 1 through the number of LC points
  dataframe_sf <- dataframe_sf %>%
    mutate(X_sequential = 1:nrow(dataframe_sf))
  
  set.seed(seed_input) #setting the seed
  
  #randomly selecting a focal point from each grid cell with trees within them
  list_grids_and_points <- st_contains(tree_grid_cropped, dataframe_sf, sparse =T) #find which points are within which grid cells, make sure row number in the data frame of grid cells corresponds to the order of the points dataframe within st_contains
  
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
  
  #filtering the dataframe with all trait data so it consists of only the focal trees --> to do this, we using the row number of overall tree point dataframe
  list_grids_and_focal_trees_df <- as.data.frame(unlist(list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe --> the row number = the cell number that focal tree was in, and the content within that observation is which row of the original dataframe the focal tree was in
  colnames(list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
  
  # making a final df of only focal trees 
  list_grids_and_focal_trees_fixed <- list_grids_and_focal_trees_df %>% 
    mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
    filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them
  
  #filtering original point data (with all the traits we care about still attached) to be just the focal trees (by the row number)
  focal_tree_dataframe <- dataframe_sf %>%
    filter(X_sequential %in% list_grids_and_focal_trees_fixed$tree_row_num) %>%
    mutate(focal_tree_row = row_number()) #adding the row number of each focal tree ...
  
  #creating the buffers around the focal points that represent the neighborhood within which we are calculating competition metrics
  focal_tree_buffers <-st_buffer(focal_tree_dataframe$geometry, 40*mean(focal_tree_dataframe$DBH_ag)) %>%
    st_as_sf() %>%
    mutate(focal_tree_row = row_number()) #create a column with row numbers
  
  #this loop iterates over each focal tree and calculates the sum for all of the neighbors of each size metric divided by the distance to the focal tree generating the competition values for each metric for each tree
  
  focal_tree_dataframe_with_competition <- tibble() #creating the empty tibble
  
  for (focal_tree in focal_tree_dataframe$focal_tree_row){ #for each focal tree 
    
    #isolate just the buffer of the given focal tree
    focal_tree_buffer <- focal_tree_buffers %>% 
      filter(focal_tree_row == focal_tree) 
    
    #return a list of trees that are within the focal tree buffer based on the row number of that tree in the dataframe_sf
    pts_in_buffer <- st_contains(focal_tree_buffer, dataframe_sf, sparse = T)[[1]] 
    
    #filtering to the data to only be the trees within the buffer.
    all_trees_in_buffer <- dataframe_sf %>%
      filter(X_sequential %in% pts_in_buffer) 
    
    #obtain the row of the focal tree in dataframe_sf 
    focal_tree_X_sequential <- focal_tree_dataframe[focal_tree,]$X_sequential
    
    #create a dataframe with only the focal tree 
    focal_tree_in_buffer <- focal_tree_dataframe %>%
      filter(X_sequential %in% focal_tree_X_sequential) 
    
    #isolating the neighbor tree data (without the focal tree)
    neighbor_trees_in_buffer <- all_trees_in_buffer %>%
      filter(X_sequential %notin% focal_tree_X_sequential) #filtering out tree data for the neighbor trees 
    
    #if there are no neighbors, sets the competition values to 0
    if(nrow(neighbor_trees_in_buffer) == 0){
      neighbor_trees_in_buffer_competiton_calc <- data.frame(
        SCA_over_distance = 0, #create a new variable for short canopy axis over distance to focal tree set to 0
        LCA_over_distance = 0, #create a new variable for long canopy axis over distance to focal tree set to 0
        CA_over_distance = 0, #create a new variable for canopy area over distance to focal tree set to 0
        CS_over_distance = 0, #create a new variable for crown spread over distance to focal tree set to 0
        DBH_over_distance = 0) #create a new variable for DBH over distance to focal tree set to 0
    } else{ #if there are indeed neighbors actually calculate the competition metric for all variables
      
      # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
      neighbor_trees_in_buffer_competiton_calc <-  neighbor_trees_in_buffer %>% 
        st_as_sf() %>%
        mutate(focal_distance = as.numeric(st_distance(geometry, focal_tree_in_buffer$geometry, by_element = F))) %>% #calculate the distance between the focal tree and each tree that neighbors 
        mutate(focal_distance = case_when(focal_distance == 0 ~ 0.0000016, 
                                          focal_distance != 0 ~ focal_distance)) %>% #replace values of 0 (if the coords are the same for multiple trees) with a value an order of magnitude smaller than the smallest distance in our dataset to avoid undefined values from division
        mutate(SCA_over_distance = Canopy_short/focal_distance) %>% #creating a column with the short canopy axis size value divided by the tree's distance from the focal tree
        mutate(LCA_over_distance = Canopy_long/focal_distance) %>%
        mutate(CA_over_distance = Canopy_area/focal_distance) %>%
        mutate(CS_over_distance = Crown_spread/focal_distance) %>%
        mutate(DBH_over_distance = DBH_ag/focal_distance) %>%
        as.data.frame() %>% #make a data frame so we don't have to deal with geometry anymore
        dplyr::select(c(SCA_over_distance, LCA_over_distance,CA_over_distance, CS_over_distance, DBH_over_distance)) %>% #keep only the columns relevant to competition metrics
        summarise(across(everything(), sum)) #get the sum of all the relevant columns
    } #closing the else
    
    focal_tree_row_with_competition <- cbind(focal_tree_in_buffer, neighbor_trees_in_buffer_competiton_calc) #bind all the info we have about the focal tree with the new competition metrics
    
    focal_tree_dataframe_with_competition <- rbind(focal_tree_dataframe_with_competition, focal_tree_row_with_competition) #bind the info about each focal tree to one another to make a final dataset with the competition metrics for all focal trees
    
  }
  
  return(list(box_sf, box_sf_cropped, dataframe_cropped,
              tree_grid_cropped, focal_tree_buffers, focal_tree_dataframe_with_competition))
}

# To test if the results of the slope tests for each population and soil metric are robust, 
  # 1) we re-run the focal_function() 500 times with a different seed input to generate a new dataframe of new randomly
        #selected focal trees and calculated competition metrics. 
  # 2) for each iteration, we find the best fitting GLS model and store the slope test results and Kendall's Tau results for later use

slope_tests <- function(population, variable) {
  
  #setting up empty lists to story slopes and p-values
  slope_permutations <- c()
  pvalue_permutations <- c()
  tau_perm <- c()
  tau_p_value_perm <- c()
  
  seed_input <- 1
  
  # for loop generating permutations of the slopes and p-values for different randomly generated p-values
  for (i in 1:500){
    
    if (population == "LM") {
      #generating the focal tree and neighbor data
      focal_results <- focal_function("LM", seed_input)
      
      #storing the focal tree dataframes
      focal_tree_dataframe_sf <- focal_results[[6]] #focal tree dataframe as a spatial object
      focal_tree_dataframe <- as.data.frame(focal_tree_dataframe_sf)  #focal tree dataframe as a spatial object
      
    } else if (population == "LC") {
      #generating the focal tree and neighbor data
      focal_results <- focal_function("LC", seed_input)
      
      #storing the focal tree dataframes
      focal_tree_dataframe_sf <- focal_results[[6]] #focal tree dataframe as a spatial object
      focal_tree_dataframe <- as.data.frame(focal_tree_dataframe_sf)  #focal tree dataframe as a spatial object
      
    } else if (population == "SD") {
      #generating the focal tree and neighbor data
      focal_results <- focal_function("SD", seed_input)
      
      #storing the focal tree dataframes
      focal_tree_dataframe_sf <- focal_results[[6]] #focal tree dataframe as a spatial object
      focal_tree_dataframe <- as.data.frame(focal_tree_dataframe_sf)  #focal tree dataframe as a spatial object
      
    }
    
    if (variable == "SCA"){
      metric = "SCA_over_distance"
      size_metric = "Canopy_short"
    } else if (variable == "LCA"){
      metric = "LCA_over_distance"
      size_metric = "Canopy_long"
    } else if (variable == "CA"){
      metric = "CA_over_distance"
      size_metric = "Canopy_area"
    } else if (variable == "CS"){
      metric = "CS_over_distance"
      size_metric = "Crown_spread"
    } else if (variable == "DBH"){
      metric = "DBH_over_distance"
      size_metric = "DBH_ag"
    } 
    
    #setting the seed for random focal tree generation 
    seed_input <- seed_input + 1
    
    #creating x and y columns of the UTM 12N 
    focal_tree_dataframe$X.1 <- st_coordinates(focal_tree_dataframe_sf)[,1]
    focal_tree_dataframe$Y <- st_coordinates(focal_tree_dataframe_sf)[,2]
    
    #Cook's D
    lm_focal <- lm(focal_tree_dataframe[[size_metric]] ~ focal_tree_dataframe[[metric]], data = focal_tree_dataframe)
    lm_focal_cooks <- cooks.distance(lm_focal) #calculating the cook.s D for each point
   # plot(lm_focal_cooks, type = 'h') #checking to see which cook's D are unsually high
    influential <- lm_focal_cooks[(lm_focal_cooks > 0.5)] #remove points with cooks D that are bigger than 3 times the mean cook's D
    
    #removing outliers based on which points were deemed influential, meaning they change the slope of the linear model too much
    if(length(influential) != 0){
      focal_tree_dataframe_no_outliers <- focal_tree_dataframe[-c(as.numeric(names(influential))),]
    } else {
      focal_tree_dataframe_no_outliers <- focal_tree_dataframe
    }
    
    
    #creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
    gls_focal_SCA <- try(gls(as.formula(paste(size_metric, "~", metric)), data = focal_tree_dataframe_no_outliers), silent = T)
    gls_focal_SCA_exp <- try(gls(as.formula(paste(size_metric, "~", metric)), correlation = corExp(form = ~X.1 + Y), data = focal_tree_dataframe_no_outliers), silent = T)
    gls_focal_SCA_gaus <- try(gls(as.formula(paste(size_metric, "~", metric)), correlation = corGaus(form = ~X.1 + Y), data = focal_tree_dataframe_no_outliers), silent = T)
    gls_focal_SCA_spher <- try(gls(as.formula(paste(size_metric, "~", metric)), correlation = corSpher(form = ~X.1 + Y), data = focal_tree_dataframe_no_outliers), silent = T)
    gls_focal_SCA_lin <- try(gls(as.formula(paste(size_metric, "~", metric)), correlation = corLin(form = ~X.1 + Y), data = focal_tree_dataframe_no_outliers), silent = T)
    gls_focal_SCA_ratio <- try(gls(as.formula(paste(size_metric, "~", metric)), correlation = corRatio(form = ~X.1 + Y), data = focal_tree_dataframe_no_outliers), silent = T)
    
    model_list <- list(
      gls_focal_SCA = gls_focal_SCA,
      gls_focal_SCA_exp = gls_focal_SCA_exp,
      gls_focal_SCA_gaus = gls_focal_SCA_gaus,
      gls_focal_SCA_spher = gls_focal_SCA_spher,
      gls_focal_SCA_lin = gls_focal_SCA_lin,
      gls_focal_SCA_ratio = gls_focal_SCA_ratio
    )
    
    # Keep only the successfully fitted models
    valid_models <- Filter(function(x) !inherits(x, "try-error"), model_list)
    
    #ordering models by which ones have the lowest Akaike information criterion, lowest AIC is the best predictive model
    AIC_test <- model.sel(valid_models) 
    
    #extracting the model that had the lowest AIC value
    summary_gls_focal <- summary(get(rownames(AIC_test[AIC_test$AICc == min(AIC_test$AICc)])))
    
    #storing slope test slope and p-value
    slope_permutations <- c(slope_permutations, summary_gls_focal$coefficients[[2]])
    pvalue_permutations <- c(pvalue_permutations, summary_gls_focal$tTable[[8]])
    
    #non parametric Kendall's Tau Test for if the the residuals are not normal or if scatterplot seems non-linear 
    tau_result <- cor.test(focal_tree_dataframe_no_outliers[[metric]], 
                           focal_tree_dataframe_no_outliers[[size_metric]],  method = "kendall")
    
    # Store Kendall's tau and its associated p-value
    tau_perm <- c(tau_perm, tau_result$statistic)
    tau_p_value_perm <- c(tau_p_value_perm, tau_result$p.value)
    

  }
  
  return(list(slope_permutations, pvalue_permutations, tau_perm, tau_p_value_perm))
}

#### Creating the Generalized Linear Effects Models ####

### LM ###

#SCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LM_SCA <- slope_tests("LM", "SCA") #focal_results <- focal_function("LM")

#save the results from the permutations
LM_SCA_slope_permutations_results <- slope_tests_LM_SCA[[1]] #slope test statistics
LM_SCA_pvalue_permutations_results <- slope_tests_LM_SCA[[2]] #p-values for the slope test
LM_SCA_tau_perm_results <- slope_tests_LM_SCA[[3]] #tau/slope statistics
LM_SCA_tau_p_value_perm_results <- slope_tests_LM_SCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LM_SCA_pvalue_permutations_results_bonf <- p.adjust(LM_SCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LM_SCA_tau_p_value_perm_results_bonf <- p.adjust(LM_SCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LM_SCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LM_SCA_slope_permutations_results) #mean
median(LM_SCA_slope_permutations_results) #median 
sd(LM_SCA_slope_permutations_results) #standard deviation
range(LM_SCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_SCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LM_SCA_pvalue_permutations_results_bonf) #mean
median(LM_SCA_pvalue_permutations_results_bonf) #median 
sd(LM_SCA_pvalue_permutations_results_bonf) #standard deviation
range(LM_SCA_pvalue_permutations_results_bonf) #range

#without Bonf correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_SCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LM_SCA_pvalue_permutations_results) #mean
median(LM_SCA_pvalue_permutations_results) #median 
sd(LM_SCA_pvalue_permutations_results) #standard deviation
range(LM_SCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LM_SCA_tau_perm_results)) +
  labs(x = "tau")

mean(LM_SCA_tau_perm_results) #mean
median(LM_SCA_tau_perm_results) #median 
sd(LM_SCA_tau_perm_results) #standard deviation
range(LM_SCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_SCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LM_SCA_tau_p_value_perm_results_bonf) #mean
median(LM_SCA_tau_p_value_perm_results_bonf) #median 
sd(LM_SCA_tau_p_value_perm_results_bonf) #standard deviation
range(LM_SCA_tau_p_value_perm_results_bonf) #range

#without correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_SCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LM_SCA_tau_p_value_perm_results) #mean
median(LM_SCA_tau_p_value_perm_results) #median 
sd(LM_SCA_tau_p_value_perm_results) #standard deviation
range(LM_SCA_tau_p_value_perm_results) #range

#LCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LM_LCA <- slope_tests("LM", "LCA") #focal_results <- focal_function("LM")

#save the results from the permutations
LM_LCA_slope_permutations_results <- slope_tests_LM_LCA[[1]] #slope test statistics
LM_LCA_pvalue_permutations_results <- slope_tests_LM_LCA[[2]] #p-values for the slope test
LM_LCA_tau_perm_results <- slope_tests_LM_LCA[[3]] #tau/slope statistics
LM_LCA_tau_p_value_perm_results <- slope_tests_LM_LCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LM_LCA_pvalue_permutations_results_bonf <- p.adjust(LM_LCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LM_LCA_tau_p_value_perm_results_bonf <- p.adjust(LM_LCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LM_LCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LM_LCA_slope_permutations_results) #mean
median(LM_LCA_slope_permutations_results) #median 
sd(LM_LCA_slope_permutations_results) #standard deviation
range(LM_LCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_LCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LM_LCA_pvalue_permutations_results_bonf) #mean
median(LM_LCA_pvalue_permutations_results_bonf) #median 
sd(LM_LCA_pvalue_permutations_results_bonf) #standard deviation
range(LM_LCA_pvalue_permutations_results_bonf) #range

#without Bonf correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_LCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LM_LCA_pvalue_permutations_results) #mean
median(LM_LCA_pvalue_permutations_results) #median 
sd(LM_LCA_pvalue_permutations_results) #standard deviation
range(LM_LCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LM_LCA_tau_perm_results)) +
  labs(x = "tau")

mean(LM_LCA_tau_perm_results) #mean
median(LM_LCA_tau_perm_results) #median 
sd(LM_LCA_tau_perm_results) #standard deviation
range(LM_LCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_LCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LM_LCA_tau_p_value_perm_results_bonf) #mean
median(LM_LCA_tau_p_value_perm_results_bonf) #median 
sd(LM_LCA_tau_p_value_perm_results_bonf) #standard deviation
range(LM_LCA_tau_p_value_perm_results_bonf) #range

#without Bonferonni Correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_LCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LM_LCA_tau_p_value_perm_results) #mean
median(LM_LCA_tau_p_value_perm_results) #median 
sd(LM_LCA_tau_p_value_perm_results) #standard deviation
range(LM_LCA_tau_p_value_perm_results) #range


#CA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LM_CA <- slope_tests("LM", "CA") #focal_results <- focal_function("LM")

#save the results from the permutations
LM_CA_slope_permutations_results <- slope_tests_LM_CA[[1]] #slope test statistics
LM_CA_pvalue_permutations_results <- slope_tests_LM_CA[[2]] #p-values for the slope test
LM_CA_tau_perm_results <- slope_tests_LM_CA[[3]] #tau/slope statistics
LM_CA_tau_p_value_perm_results <- slope_tests_LM_CA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LM_CA_pvalue_permutations_results_bonf <- p.adjust(LM_CA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LM_CA_tau_p_value_perm_results_bonf <- p.adjust(LM_CA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LM_CA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LM_CA_slope_permutations_results) #mean
median(LM_CA_slope_permutations_results) #median 
sd(LM_CA_slope_permutations_results) #standard deviation
range(LM_CA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_CA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LM_CA_pvalue_permutations_results_bonf) #mean
median(LM_CA_pvalue_permutations_results_bonf) #median 
sd(LM_CA_pvalue_permutations_results_bonf) #standard deviation
range(LM_CA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_CA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LM_CA_pvalue_permutations_results) #mean
median(LM_CA_pvalue_permutations_results) #median 
sd(LM_CA_pvalue_permutations_results) #standard deviation
range(LM_CA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LM_CA_tau_perm_results)) +
  labs(x = "tau")

mean(LM_CA_tau_perm_results) #mean
median(LM_CA_tau_perm_results) #median 
sd(LM_CA_tau_perm_results) #standard deviation
range(LM_CA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_CA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LM_CA_tau_p_value_perm_results_bonf) #mean
median(LM_CA_tau_p_value_perm_results_bonf) #median 
sd(LM_CA_tau_p_value_perm_results_bonf) #standard deviation
range(LM_CA_tau_p_value_perm_results_bonf) #range

#without Bonferonni Correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_CA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LM_CA_tau_p_value_perm_results) #mean
median(LM_CA_tau_p_value_perm_results) #median 
sd(LM_CA_tau_p_value_perm_results) #standard deviation
range(LM_CA_tau_p_value_perm_results) #range


#CS

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LM_CS <- slope_tests("LM", "CS") #focal_results <- focal_function("LM")

#save the results from the permutations
LM_CS_slope_permutations_results <- slope_tests_LM_CS[[1]] #slope test statistics
LM_CS_pvalue_permutations_results <- slope_tests_LM_CS[[2]] #p-values for the slope test
LM_CS_tau_perm_results <- slope_tests_LM_CS[[3]] #tau/slope statistics
LM_CS_tau_p_value_perm_results <- slope_tests_LM_CS[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LM_CS_pvalue_permutations_results_bonf <- p.adjust(LM_CS_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LM_CS_tau_p_value_perm_results_bonf <- p.adjust(LM_CS_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LM_CS_slope_permutations_results)) +
  labs(x = "Slope")

mean(LM_CS_slope_permutations_results) #mean
median(LM_CS_slope_permutations_results) #median 
sd(LM_CS_slope_permutations_results) #standard deviation
range(LM_CS_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_CS_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LM_CS_pvalue_permutations_results_bonf) #mean
median(LM_CS_pvalue_permutations_results_bonf) #median 
sd(LM_CS_pvalue_permutations_results_bonf) #standard deviation
range(LM_CS_pvalue_permutations_results_bonf) #range

#without Bonferonni Correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_CS_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LM_CS_pvalue_permutations_results) #mean
median(LM_CS_pvalue_permutations_results) #median 
sd(LM_CS_pvalue_permutations_results) #standard deviation
range(LM_CS_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LM_CS_tau_perm_results)) +
  labs(x = "tau")

mean(LM_CS_tau_perm_results) #mean
median(LM_CS_tau_perm_results) #median 
sd(LM_CS_tau_perm_results) #standard deviation
range(LM_CS_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_CS_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LM_CS_tau_p_value_perm_results_bonf) #mean
median(LM_CS_tau_p_value_perm_results_bonf) #median 
sd(LM_CS_tau_p_value_perm_results_bonf) #standard deviation
range(LM_CS_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_CS_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LM_CS_tau_p_value_perm_results) #mean
median(LM_CS_tau_p_value_perm_results) #median 
sd(LM_CS_tau_p_value_perm_results) #standard deviation
range(LM_CS_tau_p_value_perm_results) #range


#DBH

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LM_DBH <- slope_tests("LM", "DBH") #focal_results <- focal_function("LM")

#save the results from the permutations
LM_DBH_slope_permutations_results <- slope_tests_LM_DBH[[1]] #slope test statistics
LM_DBH_pvalue_permutations_results <- slope_tests_LM_DBH[[2]] #p-values for the slope test
LM_DBH_tau_perm_results <- slope_tests_LM_DBH[[3]] #tau/slope statistics
LM_DBH_tau_p_value_perm_results <- slope_tests_LM_DBH[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LM_DBH_pvalue_permutations_results_bonf <- p.adjust(LM_DBH_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LM_DBH_tau_p_value_perm_results_bonf <- p.adjust(LM_DBH_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LM_DBH_slope_permutations_results)) +
  labs(x = "Slope")

mean(LM_DBH_slope_permutations_results) #mean
median(LM_DBH_slope_permutations_results) #median 
sd(LM_DBH_slope_permutations_results) #standard deviation
range(LM_DBH_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_DBH_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LM_DBH_pvalue_permutations_results_bonf) #mean
median(LM_DBH_pvalue_permutations_results_bonf) #median 
sd(LM_DBH_pvalue_permutations_results_bonf) #standard deviation
range(LM_DBH_pvalue_permutations_results_bonf) #range

#without Bonferroni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LM_DBH_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LM_DBH_pvalue_permutations_results) #mean
median(LM_DBH_pvalue_permutations_results) #median 
sd(LM_DBH_pvalue_permutations_results) #standard deviation
range(LM_DBH_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LM_DBH_tau_perm_results)) +
  labs(x = "tau")

mean(LM_DBH_tau_perm_results) #mean
median(LM_DBH_tau_perm_results) #median 
sd(LM_DBH_tau_perm_results) #standard deviation
range(LM_DBH_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_DBH_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LM_DBH_tau_p_value_perm_results_bonf) #mean
median(LM_DBH_tau_p_value_perm_results_bonf) #median 
sd(LM_DBH_tau_p_value_perm_results_bonf) #standard deviation
range(LM_DBH_tau_p_value_perm_results_bonf) #range

#without Bonferroni corrections
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LM_DBH_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LM_DBH_tau_p_value_perm_results) #mean
median(LM_DBH_tau_p_value_perm_results) #median 
sd(LM_DBH_tau_p_value_perm_results) #standard deviation
range(LM_DBH_tau_p_value_perm_results) #range


### LC ###

#SCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LC_LCA <- slope_tests("LC", "SCA") #focal_results <- focal_function("LC")

#save the results from the permutations
LC_SCA_slope_permutations_results <- slope_tests_LC_LCA[[1]] #slope test statistics
LC_SCA_pvalue_permutations_results <- slope_tests_LC_LCA[[2]] #p-values for the slope test
LC_SCA_tau_perm_results <- slope_tests_LC_LCA[[3]] #tau/slope statistics
LC_SCA_tau_p_value_perm_results <- slope_tests_LC_LCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LC_SCA_pvalue_permutations_results_bonf <- p.adjust(LC_SCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LC_SCA_tau_p_value_perm_results_bonf <- p.adjust(LC_SCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LC_SCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LC_SCA_slope_permutations_results) #mean
median(LC_SCA_slope_permutations_results) #median 
sd(LC_SCA_slope_permutations_results) #standard deviation
range(LC_SCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_SCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LC_SCA_pvalue_permutations_results_bonf) #mean
median(LC_SCA_pvalue_permutations_results_bonf) #median 
sd(LC_SCA_pvalue_permutations_results_bonf) #standard deviation
range(LC_SCA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_SCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LC_SCA_pvalue_permutations_results) #mean
median(LC_SCA_pvalue_permutations_results) #median 
sd(LC_SCA_pvalue_permutations_results) #standard deviation
range(LC_SCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LC_SCA_tau_perm_results)) +
  labs(x = "tau")

mean(LC_SCA_tau_perm_results) #mean
median(LC_SCA_tau_perm_results) #median 
sd(LC_SCA_tau_perm_results) #standard deviation
range(LC_SCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_SCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LC_SCA_tau_p_value_perm_results_bonf) #mean
median(LC_SCA_tau_p_value_perm_results_bonf) #median 
sd(LC_SCA_tau_p_value_perm_results_bonf) #standard deviation
range(LC_SCA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_SCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LC_SCA_tau_p_value_perm_results) #mean
median(LC_SCA_tau_p_value_perm_results) #median 
sd(LC_SCA_tau_p_value_perm_results) #standard deviation
range(LC_SCA_tau_p_value_perm_results) #range

#LCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LC_LCA <- slope_tests("LC", "LCA") #focal_results <- focal_function("LC")

#save the results from the permutations
LC_LCA_slope_permutations_results <- slope_tests_LC_LCA[[1]] #slope test statistics
LC_LCA_pvalue_permutations_results <- slope_tests_LC_LCA[[2]] #p-values for the slope test
LC_LCA_tau_perm_results <- slope_tests_LC_LCA[[3]] #tau/slope statistics
LC_LCA_tau_p_value_perm_results <- slope_tests_LC_LCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LC_LCA_pvalue_permutations_results_bonf <- p.adjust(LC_LCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LC_LCA_tau_p_value_perm_results_bonf <- p.adjust(LC_LCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LC_LCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LC_LCA_slope_permutations_results) #mean
median(LC_LCA_slope_permutations_results) #median 
sd(LC_LCA_slope_permutations_results) #standard deviation
range(LC_LCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_LCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LC_LCA_pvalue_permutations_results_bonf) #mean
median(LC_LCA_pvalue_permutations_results_bonf) #median 
sd(LC_LCA_pvalue_permutations_results_bonf) #standard deviation
range(LC_LCA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_LCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LC_LCA_pvalue_permutations_results) #mean
median(LC_LCA_pvalue_permutations_results) #median 
sd(LC_LCA_pvalue_permutations_results) #standard deviation
range(LC_LCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LC_LCA_tau_perm_results)) +
  labs(x = "tau")

mean(LC_LCA_tau_perm_results) #mean
median(LC_LCA_tau_perm_results) #median 
sd(LC_LCA_tau_perm_results) #standard deviation
range(LC_LCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_LCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LC_LCA_tau_p_value_perm_results_bonf) #mean
median(LC_LCA_tau_p_value_perm_results_bonf) #median 
sd(LC_LCA_tau_p_value_perm_results_bonf) #standard deviation
range(LC_LCA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_LCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LC_LCA_tau_p_value_perm_results) #mean
median(LC_LCA_tau_p_value_perm_results) #median 
sd(LC_LCA_tau_p_value_perm_results) #standard deviation
range(LC_LCA_tau_p_value_perm_results) #range


#CA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LC_CA <- slope_tests("LC", "CA") #focal_results <- focal_function("LC")

#save the results from the permutations
LC_CA_slope_permutations_results <- slope_tests_LC_CA[[1]] #slope test statistics
LC_CA_pvalue_permutations_results <- slope_tests_LC_CA[[2]] #p-values for the slope test
LC_CA_tau_perm_results <- slope_tests_LC_CA[[3]] #tau/slope statistics
LC_CA_tau_p_value_perm_results <- slope_tests_LC_CA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LC_CA_pvalue_permutations_results_bonf <- p.adjust(LC_CA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LC_CA_tau_p_value_perm_results_bonf <- p.adjust(LC_CA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LC_CA_slope_permutations_results)) +
  labs(x = "Slope")

mean(LC_CA_slope_permutations_results) #mean
median(LC_CA_slope_permutations_results) #median 
sd(LC_CA_slope_permutations_results) #standard deviation
range(LC_CA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_CA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LC_CA_pvalue_permutations_results_bonf) #mean
median(LC_CA_pvalue_permutations_results_bonf) #median 
sd(LC_CA_pvalue_permutations_results_bonf) #standard deviation
range(LC_CA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_CA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LC_CA_pvalue_permutations_results) #mean
median(LC_CA_pvalue_permutations_results) #median 
sd(LC_CA_pvalue_permutations_results) #standard deviation
range(LC_CA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LC_CA_tau_perm_results)) +
  labs(x = "tau")

mean(LC_CA_tau_perm_results) #mean
median(LC_CA_tau_perm_results) #median 
sd(LC_CA_tau_perm_results) #standard deviation
range(LC_CA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_CA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LC_CA_tau_p_value_perm_results_bonf) #mean
median(LC_CA_tau_p_value_perm_results_bonf) #median 
sd(LC_CA_tau_p_value_perm_results_bonf) #standard deviation
range(LC_CA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_CA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LC_CA_tau_p_value_perm_results) #mean
median(LC_CA_tau_p_value_perm_results) #median 
sd(LC_CA_tau_p_value_perm_results) #standard deviation
range(LC_CA_tau_p_value_perm_results) #range


#CS

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LC_CS <- slope_tests("LC", "CS") #focal_results <- focal_function("LC")

#save the results from the permutations
LC_CS_slope_permutations_results <- slope_tests_LC_CS[[1]] #slope test statistics
LC_CS_pvalue_permutations_results <- slope_tests_LC_CS[[2]] #p-values for the slope test
LC_CS_tau_perm_results <- slope_tests_LC_CS[[3]] #tau/slope statistics
LC_CS_tau_p_value_perm_results <- slope_tests_LC_CS[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LC_CS_pvalue_permutations_results_bonf <- p.adjust(LC_CS_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LC_CS_tau_p_value_perm_results_bonf <- p.adjust(LC_CS_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LC_CS_slope_permutations_results)) +
  labs(x = "Slope")

mean(LC_CS_slope_permutations_results) #mean
median(LC_CS_slope_permutations_results) #median 
sd(LC_CS_slope_permutations_results) #standard deviation
range(LC_CS_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_CS_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LC_CS_pvalue_permutations_results_bonf) #mean
median(LC_CS_pvalue_permutations_results_bonf) #median 
sd(LC_CS_pvalue_permutations_results_bonf) #standard deviation
range(LC_CS_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_CS_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LC_CS_pvalue_permutations_results) #mean
median(LC_CS_pvalue_permutations_results) #median 
sd(LC_CS_pvalue_permutations_results) #standard deviation
range(LC_CS_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LC_CS_tau_perm_results)) +
  labs(x = "tau")

mean(LC_CS_tau_perm_results) #mean
median(LC_CS_tau_perm_results) #median 
sd(LC_CS_tau_perm_results) #standard deviation
range(LC_CS_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_CS_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LC_CS_tau_p_value_perm_results_bonf) #mean
median(LC_CS_tau_p_value_perm_results_bonf) #median 
sd(LC_CS_tau_p_value_perm_results_bonf) #standard deviation
range(LC_CS_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_CS_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LC_CS_tau_p_value_perm_results) #mean
median(LC_CS_tau_p_value_perm_results) #median 
sd(LC_CS_tau_p_value_perm_results) #standard deviation
range(LC_CS_tau_p_value_perm_results) #range


#DBH

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_LC_DBH <- slope_tests("LC", "DBH") #focal_results <- focal_function("LC")

#save the results from the permutations
LC_DBH_slope_permutations_results <- slope_tests_LC_DBH[[1]] #slope test statistics
LC_DBH_pvalue_permutations_results <- slope_tests_LC_DBH[[2]] #p-values for the slope test
LC_DBH_tau_perm_results <- slope_tests_LC_DBH[[3]] #tau/slope statistics
LC_DBH_tau_p_value_perm_results <- slope_tests_LC_DBH[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
LC_DBH_pvalue_permutations_results_bonf <- p.adjust(LC_DBH_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
LC_DBH_tau_p_value_perm_results_bonf <- p.adjust(LC_DBH_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = LC_DBH_slope_permutations_results)) +
  labs(x = "Slope")

mean(LC_DBH_slope_permutations_results) #mean
median(LC_DBH_slope_permutations_results) #median 
sd(LC_DBH_slope_permutations_results) #standard deviation
range(LC_DBH_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_DBH_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(LC_DBH_pvalue_permutations_results_bonf) #mean
median(LC_DBH_pvalue_permutations_results_bonf) #median 
sd(LC_DBH_pvalue_permutations_results_bonf) #standard deviation
range(LC_DBH_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = LC_DBH_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(LC_DBH_pvalue_permutations_results) #mean
median(LC_DBH_pvalue_permutations_results) #median 
sd(LC_DBH_pvalue_permutations_results) #standard deviation
range(LC_DBH_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = LC_DBH_tau_perm_results)) +
  labs(x = "tau")

mean(LC_DBH_tau_perm_results) #mean
median(LC_DBH_tau_perm_results) #median 
sd(LC_DBH_tau_perm_results) #standard deviation
range(LC_DBH_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_DBH_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(LC_DBH_tau_p_value_perm_results_bonf) #mean
median(LC_DBH_tau_p_value_perm_results_bonf) #median 
sd(LC_DBH_tau_p_value_perm_results_bonf) #standard deviation
range(LC_DBH_tau_p_value_perm_results_bonf) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = LC_DBH_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(LC_DBH_tau_p_value_perm_results) #mean
median(LC_DBH_tau_p_value_perm_results) #median 
sd(LC_DBH_tau_p_value_perm_results) #standard deviation
range(LC_DBH_tau_p_value_perm_results) #range

### SD ###

#SCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_SD_SCA <- slope_tests("SD", "SCA") #focal_results <- focal_function("SD")

#save the results from the permutations
SD_SCA_slope_permutations_results <- slope_tests_SD_SCA[[1]] #slope test statistics
SD_SCA_pvalue_permutations_results <- slope_tests_SD_SCA[[2]] #p-values for the slope test
SD_SCA_tau_perm_results <- slope_tests_SD_SCA[[3]] #tau/slope statistics
SD_SCA_tau_p_value_perm_results <- slope_tests_SD_SCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
SD_SCA_pvalue_permutations_results_bonf <- p.adjust(SD_SCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
SD_SCA_tau_p_value_perm_results_bonf <- p.adjust(SD_SCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = SD_SCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(SD_SCA_slope_permutations_results) #mean
median(SD_SCA_slope_permutations_results) #median 
sd(SD_SCA_slope_permutations_results) #standard deviation
range(SD_SCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_SCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(SD_SCA_pvalue_permutations_results_bonf) #mean
median(SD_SCA_pvalue_permutations_results_bonf) #median 
sd(SD_SCA_pvalue_permutations_results_bonf) #standard deviation
range(SD_SCA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_SCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(SD_SCA_pvalue_permutations_results) #mean
median(SD_SCA_pvalue_permutations_results) #median 
sd(SD_SCA_pvalue_permutations_results) #standard deviation
range(SD_SCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = SD_SCA_tau_perm_results)) +
  labs(x = "tau")

mean(SD_SCA_tau_perm_results) #mean
median(SD_SCA_tau_perm_results) #median 
sd(SD_SCA_tau_perm_results) #standard deviation
range(SD_SCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_SCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(SD_SCA_tau_p_value_perm_results_bonf) #mean
median(SD_SCA_tau_p_value_perm_results_bonf) #median 
sd(SD_SCA_tau_p_value_perm_results_bonf) #standard deviation
range(SD_SCA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_SCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(SD_SCA_tau_p_value_perm_results) #mean
median(SD_SCA_tau_p_value_perm_results) #median 
sd(SD_SCA_tau_p_value_perm_results) #standard deviation
range(SD_SCA_tau_p_value_perm_results) #range

#LCA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_SD_LCA <- slope_tests("SD", "LCA") #focal_results <- focal_function("LC")

#save the results from the permutations
SD_LCA_slope_permutations_results <- slope_tests_SD_LCA[[1]] #slope test statistics
SD_LCA_pvalue_permutations_results <- slope_tests_SD_LCA[[2]] #p-values for the slope test
SD_LCA_tau_perm_results <- slope_tests_SD_LCA[[3]] #tau/slope statistics
SD_LCA_tau_p_value_perm_results <- slope_tests_SD_LCA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
SD_LCA_pvalue_permutations_results_bonf <- p.adjust(SD_LCA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
SD_LCA_tau_p_value_perm_results_bonf <- p.adjust(SD_LCA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = SD_LCA_slope_permutations_results)) +
  labs(x = "Slope")

mean(SD_LCA_slope_permutations_results) #mean
median(SD_LCA_slope_permutations_results) #median 
sd(SD_LCA_slope_permutations_results) #standard deviation
range(SD_LCA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_LCA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(SD_LCA_pvalue_permutations_results_bonf) #mean
median(SD_LCA_pvalue_permutations_results_bonf) #median 
sd(SD_LCA_pvalue_permutations_results_bonf) #standard deviation
range(SD_LCA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_LCA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(SD_LCA_pvalue_permutations_results) #mean
median(SD_LCA_pvalue_permutations_results) #median 
sd(SD_LCA_pvalue_permutations_results) #standard deviation
range(SD_LCA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = SD_LCA_tau_perm_results)) +
  labs(x = "tau")

mean(SD_LCA_tau_perm_results) #mean
median(SD_LCA_tau_perm_results) #median 
sd(SD_LCA_tau_perm_results) #standard deviation
range(SD_LCA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_LCA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(SD_LCA_tau_p_value_perm_results_bonf) #mean
median(SD_LCA_tau_p_value_perm_results_bonf) #median 
sd(SD_LCA_tau_p_value_perm_results_bonf) #standard deviation
range(SD_LCA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_LCA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(SD_LCA_tau_p_value_perm_results) #mean
median(SD_LCA_tau_p_value_perm_results) #median 
sd(SD_LCA_tau_p_value_perm_results) #standard deviation
range(SD_LCA_tau_p_value_perm_results) #range


#CA

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_SD_CA <- slope_tests("SD", "CA") #focal_results <- focal_function("LC")

#save the results from the permutations
SD_CA_slope_permutations_results <- slope_tests_SD_CA[[1]] #slope test statistics
SD_CA_pvalue_permutations_results <- slope_tests_SD_CA[[2]] #p-values for the slope test
SD_CA_tau_perm_results <- slope_tests_SD_CA[[3]] #tau/slope statistics
SD_CA_tau_p_value_perm_results <- slope_tests_SD_CA[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
SD_CA_pvalue_permutations_results_bonf <- p.adjust(SD_CA_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
SD_CA_tau_p_value_perm_results_bonf <- p.adjust(SD_CA_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = SD_CA_slope_permutations_results)) +
  labs(x = "Slope")

mean(SD_CA_slope_permutations_results) #mean
median(SD_CA_slope_permutations_results) #median 
sd(SD_CA_slope_permutations_results) #standard deviation
range(SD_CA_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_CA_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(SD_CA_pvalue_permutations_results_bonf) #mean
median(SD_CA_pvalue_permutations_results_bonf) #median 
sd(SD_CA_pvalue_permutations_results_bonf) #standard deviation
range(SD_CA_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_CA_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(SD_CA_pvalue_permutations_results) #mean
median(SD_CA_pvalue_permutations_results) #median 
sd(SD_CA_pvalue_permutations_results) #standard deviation
range(SD_CA_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = SD_CA_tau_perm_results)) +
  labs(x = "tau")

mean(SD_CA_tau_perm_results) #mean
median(SD_CA_tau_perm_results) #median 
sd(SD_CA_tau_perm_results) #standard deviation
range(SD_CA_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_CA_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(SD_CA_tau_p_value_perm_results_bonf) #mean
median(SD_CA_tau_p_value_perm_results_bonf) #median 
sd(SD_CA_tau_p_value_perm_results_bonf) #standard deviation
range(SD_CA_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_CA_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(SD_CA_tau_p_value_perm_results) #mean
median(SD_CA_tau_p_value_perm_results) #median 
sd(SD_CA_tau_p_value_perm_results) #standard deviation
range(SD_CA_tau_p_value_perm_results) #range


#CS

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests_SD_CS <- slope_tests("SD", "CS") #focal_results <- focal_function("LC")

#save the results from the permutations
SD_CS_slope_permutations_results <- slope_tests_SD_CS[[1]] #slope test statistics
SD_CS_pvalue_permutations_results <- slope_tests_SD_CS[[2]] #p-values for the slope test
SD_CS_tau_perm_results <- slope_tests_SD_CS[[3]] #tau/slope statistics
SD_CS_tau_p_value_perm_results <- slope_tests_SD_CS[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
SD_CS_pvalue_permutations_results_bonf <- p.adjust(SD_CS_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
SD_CS_tau_p_value_perm_results_bonf <- p.adjust(SD_CS_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = SD_CS_slope_permutations_results)) +
  labs(x = "Slope")

mean(SD_CS_slope_permutations_results) #mean
median(SD_CS_slope_permutations_results) #median 
sd(SD_CS_slope_permutations_results) #standard deviation
range(SD_CS_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_CS_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(SD_CS_pvalue_permutations_results_bonf) #mean
median(SD_CS_pvalue_permutations_results_bonf) #median 
sd(SD_CS_pvalue_permutations_results_bonf) #standard deviation
range(SD_CS_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_CS_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(SD_CS_pvalue_permutations_results) #mean
median(SD_CS_pvalue_permutations_results) #median 
sd(SD_CS_pvalue_permutations_results) #standard deviation
range(SD_CS_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = SD_CS_tau_perm_results)) +
  labs(x = "tau")

mean(SD_CS_tau_perm_results) #mean
median(SD_CS_tau_perm_results) #median 
sd(SD_CS_tau_perm_results) #standard deviation
range(SD_CS_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_CS_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(SD_CS_tau_p_value_perm_results_bonf) #mean
median(SD_CS_tau_p_value_perm_results_bonf) #median 
sd(SD_CS_tau_p_value_perm_results_bonf) #standard deviation
range(SD_CS_tau_p_value_perm_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_CS_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(SD_CS_tau_p_value_perm_results) #mean
median(SD_CS_tau_p_value_perm_results) #median 
sd(SD_CS_tau_p_value_perm_results) #standard deviation
range(SD_CS_tau_p_value_perm_results) #range


#DBH

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
slope_tests <- slope_tests("SD", "DBH") #focal_results <- focal_function("LC")

#save the results from the permutations
SD_DBH_slope_permutations_results <- slope_tests[[1]] #slope test statistics
SD_DBH_pvalue_permutations_results <- slope_tests[[2]] #p-values for the slope test
SD_DBH_tau_perm_results <- slope_tests[[3]] #tau/slope statistics
SD_DBH_tau_p_value_perm_results <- slope_tests[[4]] #p-value statistic for the Kendall's Tau Test

#use a Bonferroni Correction to control for multiple testing error
SD_DBH_pvalue_permutations_results_bonf <- p.adjust(SD_DBH_pvalue_permutations_results, method = "bonferroni") #p-values for the slope test
SD_DBH_tau_p_value_perm_results_bonf <- p.adjust(SD_DBH_tau_p_value_perm_results, method = "bonferroni") #p-value statistic for the Kendall's Tau Test

#Looking at the distribution and descriptive summary of the slope test statistics
ggplot()+
  geom_histogram(aes(x = SD_DBH_slope_permutations_results)) +
  labs(x = "Slope")

mean(SD_DBH_slope_permutations_results) #mean
median(SD_DBH_slope_permutations_results) #median 
sd(SD_DBH_slope_permutations_results) #standard deviation
range(SD_DBH_slope_permutations_results) #range

#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_DBH_pvalue_permutations_results_bonf)) +
  labs(x = "Slope Test P-Values")

mean(SD_DBH_pvalue_permutations_results_bonf) #mean
median(SD_DBH_pvalue_permutations_results_bonf) #median 
sd(SD_DBH_pvalue_permutations_results_bonf) #standard deviation
range(SD_DBH_pvalue_permutations_results_bonf) #range

#without Bonferonni correction
#Looking at the distribution and descriptive summary of the p-values for the slope test 
ggplot()+
  geom_histogram(aes(x = SD_DBH_pvalue_permutations_results)) +
  labs(x = "Slope Test P-Values")

mean(SD_DBH_pvalue_permutations_results) #mean
median(SD_DBH_pvalue_permutations_results) #median 
sd(SD_DBH_pvalue_permutations_results) #standard deviation
range(SD_DBH_pvalue_permutations_results) #range

#Looking at the distribution and descriptive summary of the tau/slope statistic
ggplot()+
  geom_histogram(aes(x = SD_DBH_tau_perm_results)) +
  labs(x = "tau")

mean(SD_DBH_tau_perm_results) #mean
median(SD_DBH_tau_perm_results) #median 
sd(SD_DBH_tau_perm_results) #standard deviation
range(SD_DBH_tau_perm_results) #range

#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_DBH_tau_p_value_perm_results_bonf)) +
  labs(x = "tau p-value")

mean(SD_DBH_tau_p_value_perm_results_bonf) #mean
median(SD_DBH_tau_p_value_perm_results_bonf) #median 
sd(SD_DBH_tau_p_value_perm_results_bonf) #standard deviation
range(SD_DBH_tau_p_value_perm_results_bonf) #range

#without correction
#Looking at the distribution and descriptive summary of the p-value statistic for the Kendall's Tau Test
ggplot()+
  geom_histogram(aes(x = SD_DBH_tau_p_value_perm_results)) +
  labs(x = "tau p-value")

mean(SD_DBH_tau_p_value_perm_results) #mean
median(SD_DBH_tau_p_value_perm_results) #median 
sd(SD_DBH_tau_p_value_perm_results) #standard deviation
range(SD_DBH_tau_p_value_perm_results) #range



