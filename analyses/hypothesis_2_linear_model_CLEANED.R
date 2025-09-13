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
  tree_grid_cropped <- st_make_grid(dataframe_cropped, cellsize = (((40*mean(dataframe$DBH_ag))*4))) 
  
  #creating an x_sequential column that is 1 through the number of LC points
  dataframe_sf <- dataframe_sf %>%
    mutate(X_sequential = 1:nrow(dataframe_sf))
  
  #setting the seed
  set.seed(25) 
  
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

####LM ####

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
focal_results <- focal_function("LM")

#assigning necessary dataframes and objects 
LM_box_sf <- focal_results[[1]] #bounding box
LM_box_sf_cropped <- focal_results[[2]] #bounding box cropped by 20m
LM_fixed_field_data_processed_sf_cropped <- focal_results[[3]] # cropped tree data
LM_tree_grid_cropped <- focal_results[[4]] #grid with 40*mean population DBH as grid size
LM_focal_tree_buffers <- focal_results[[5]] #focal tree buffers
LM_focal_tree_dataframe_sf <- focal_results[[6]] #focal tree dataframe as a spatial object
LM_focal_tree_dataframe <- as.data.frame(LM_focal_tree_dataframe_sf)  #focal tree dataframe as a spatial object

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
  geom_sf(data= LM_focal_tree_dataframe_sf, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LM_tree_grid_cropped) +
  geom_sf(data=LM_focal_tree_buffers, aes(color = focal_tree_row))+
  geom_sf(data=LM_focal_tree_dataframe_sf, aes(color = focal_tree_row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = LM_focal_tree_buffers)+
  geom_sf(data = LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_focal_tree_dataframe_sf, color = 'blue')


#descriptive statistics for the focal tree sum of size/shape metrics over distance

#histograms
ggplot(LM_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LM_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LM_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(LM_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(LM_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
LM_field_data_focal_summarized_focal <- LM_focal_tree_dataframe %>%
  dplyr::select(SCA_over_distance, LCA_over_distance, CS_over_distance, CA_over_distance, DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(LM_field_data_focal_summarized_focal)









#8 categories for direction


#all points 

#removing NAs in SCA and aspect to allow us to run tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_aspect_raster_15_data_pts_8_categorical)

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
all_points_aov_SCA_aspect_8 <- aov(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_SCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_SCA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test
hist(all_points_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test, if significant, have to run a non-parametric test

#residuals are not normal

# checking equal variances with Fligner-Killeen Test, Levene's Rest, and Rule of Thumb Test

#Fligner-Killeen, more useful than the other equal variance tests when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_SCA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                    all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_SCA, na.rm = T) / min(all_points_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two, the test did not pass

#variances are equally distributed

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

















#### creating the generalized linear effects model ####

#conditions are lINES: linearity, independence, normal distribution of residuals, equal variance, simple random sample

#creating x and y columns of the UTM 12N 
LM_focal_tree_dataframe$X.1 <- st_coordinates(LM_focal_tree_dataframe_sf)[,1]
LM_focal_tree_dataframe$Y <- st_coordinates(LM_focal_tree_dataframe_sf)[,2]

#SCA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_focal_tree_dataframe, (aes(x=SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
LM_lm_focal_SCA <- lm(Canopy_short ~ SCA_over_distance, data = LM_focal_tree_dataframe)
LM_lm_focal_SCA_cooks <- cooks.distance(LM_lm_focal_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential, meaning they change the slope of the linear model too much
LM_focal_tree_dataframe_no_SCA_outliers <- LM_focal_tree_dataframe[-c(13),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_SCA <- gls(Canopy_short ~ SCA_over_distance, data = LM_focal_tree_dataframe)
LM_gls_focal_SCA_exp <- gls(Canopy_short ~ SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_SCA_gaus <- gls(Canopy_short ~ SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_SCA_spher <- gls(Canopy_short ~ SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_SCA_lin <- gls(Canopy_short ~ SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_SCA_ratio <- gls(Canopy_short ~ SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion, lowest AIC is the best predictive model
LM_AIC_test <- model.sel(LM_gls_focal_SCA, LM_gls_focal_SCA_exp, LM_gls_focal_SCA_gaus, LM_gls_focal_SCA_spher, LM_gls_focal_SCA_ratio) #LM_gls_focal_SCA_lin
LM_AIC_test

# While LM_gls_focal_SCA_lin has the lowest AIC, but we have had trouble with it, so we are using the second best option LM_gls_focal_SCA 

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_focal_tree_dataframe, aes(x= LM_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plot
ggplot(LM_focal_tree_dataframe, aes(sample = LM_gls_focal_SCA$residuals))+
  geom_qq()

#Shapiro test
shapiro.test(LM_gls_focal_SCA$residuals) #not significant so it is normal 

#checking equal variance
ggplot(data = LM_focal_tree_dataframe, aes(x = LM_gls_focal_SCA$fitted, y = LM_gls_focal_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#checking we have appropriately removed the spatial autocorrelation
semivario <- Variogram( LM_gls_focal_SCA_ratio, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LM_gls_focal_SCA)

#non parametric Kendall's Tau Test for if the the residuals are not normal or if scatterplot seems non-linear 
LM_tau_result_SCA <- cor.test(LM_focal_tree_dataframe$SCA_over_distance, 
                              LM_focal_tree_dataframe$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_SCA)

# Calculate the trend line
LM_trend_line_SCA <- predict(loess(LM_focal_tree_dataframe$Canopy_short ~ LM_focal_tree_dataframe$SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_focal_tree_dataframe$SCA_over_distance, y = (LM_focal_tree_dataframe$Canopy_short), color = "blue")) +
  geom_line(aes(x = LM_focal_tree_dataframe$SCA_over_distance, y = LM_trend_line_SCA), color = "red") +
  labs(x = "SCA over Distance", y = "Short Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#LCA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_focal_tree_dataframe, (aes(x=LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
LM_lm_focal_LCA <- lm(Canopy_long ~ LCA_over_distance, data = LM_focal_tree_dataframe)
LM_lm_focal_LCA_cooks <- cooks.distance(LM_lm_focal_LCA) #calculating the cook.s D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_LCA_cooks[(LM_lm_focal_LCA_cooks > (3 * mean(LM_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_focal_tree_dataframe_no_LCA_outliers <- LM_focal_tree_dataframe[-c(13),]

#Cook's D
LM_lm_focal_LCA <- lm(Canopy_long ~ LCA_over_distance, data = LM_focal_tree_dataframe_no_LCA_outliers)
LM_lm_focal_LCA_cooks <- cooks.distance(LM_lm_focal_LCA) #calculating the cook.s D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_LCA_cooks[(LM_lm_focal_LCA_cooks > (3 * mean(LM_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_focal_tree_dataframe_no_LCA_outliers <- LM_focal_tree_dataframe[-c(9, 13),]


#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_LCA <- gls(Canopy_long ~ LCA_over_distance, data = LM_focal_tree_dataframe)
LM_gls_focal_LCA_exp <- gls(Canopy_long ~ LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_LCA_gaus <- gls(Canopy_long ~ LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_LCA_spher <- gls(Canopy_long ~ LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_LCA_lin <- gls(Canopy_long ~ LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_LCA_ratio <- gls(Canopy_long ~ LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_LCA <- model.sel(LM_gls_focal_LCA, LM_gls_focal_LCA_exp, LM_gls_focal_LCA_gaus, LM_gls_focal_LCA_spher, LM_gls_focal_LCA_lin, LM_gls_focal_LCA_ratio)
LM_AIC_test_LCA

#LM_gls_focal_LCA has lowest AIC

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_focal_tree_dataframe, aes(x= LM_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plots
ggplot(LM_focal_tree_dataframe, aes(sample = LM_gls_focal_LCA$residuals))+
  geom_qq()

#shaprio wilk test, not significant so our residuals are normally distributed
shapiro.test(LM_gls_focal_LCA$residuals) 

#checking equal variance
ggplot(data = LM_focal_tree_dataframe , aes(x = LM_gls_focal_LCA$fitted, y = LM_gls_focal_LCA$residuals))+
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
LM_tau_result_LCA <- cor.test(LM_focal_tree_dataframe$LCA_over_distance, 
                              LM_focal_tree_dataframe$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_LCA)

# Calculate the trend line
LM_trend_line_LCA <- predict(loess(LM_focal_tree_dataframe$Canopy_long ~ LM_focal_tree_dataframe$SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_focal_tree_dataframe$LCA_over_distance, y = (LM_focal_tree_dataframe$Canopy_long), color = "blue")) +
  geom_line(aes(x = LM_focal_tree_dataframe$LCA_over_distance, y = LM_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#CA

#plotting the linear model in ggplot for SCA
ggplot(data = LM_focal_tree_dataframe, (aes(x=CA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
LM_lm_focal_CA <- lm(Canopy_area ~ CA_over_distance, data = LM_focal_tree_dataframe)
LM_lm_focal_CA_cooks <- cooks.distance(LM_lm_focal_CA) #calculating the cook.s D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_CA_cooks[(LM_lm_focal_CA_cooks > (3 * mean(LM_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_focal_tree_dataframe_no_CA_outliers <- LM_focal_tree_dataframe[-c(9, 13),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_CA <- gls(Canopy_long ~ CA_over_distance, data = LM_focal_tree_dataframe)
LM_gls_focal_CA_exp <- gls(Canopy_short ~ CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CA_gaus <- gls(Canopy_short ~ CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CA_spher <- gls(Canopy_short ~ CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CA_lin <- gls(Canopy_short ~ CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CA_ratio <- gls(Canopy_short ~ CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CA <- model.sel(LM_gls_focal_CA, LM_gls_focal_CA_exp, LM_gls_focal_CA_gaus, LM_gls_focal_CA_spher, LM_gls_focal_CA_ratio) #LM_gls_focal_CA_lin
LM_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_focal_tree_dataframe, aes(x= LM_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_focal_tree_dataframe, aes(sample = LM_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, not sign so the residuals are normally distributed
shapiro.test(LM_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LM_focal_tree_dataframe , aes(x = LM_gls_focal_CA$fitted, y = LM_gls_focal_CA$residuals))+
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
LM_tau_result_CA <- cor.test(LM_focal_tree_dataframe$CA_over_distance, LM_focal_tree_dataframe$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CA)

# Calculate the trend line
LM_trend_line_CA <- predict(loess(LM_focal_tree_dataframe$Canopy_area ~ LM_focal_tree_dataframe$CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_focal_tree_dataframe$CA_over_distance, y = (LM_focal_tree_dataframe$Canopy_area), color = "blue")) +
  geom_line(aes(x = LM_focal_tree_dataframe$CA_over_distance, y = LM_trend_line_CA), color = "red") +
  labs(x = "CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot for SCA
ggplot(data = LM_focal_tree_dataframe, (aes(x=CS_over_distance, y=Crown_spread)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CS over Distance")+
  ylab("Crown Spread")

#Cook's D
LM_lm_focal_CS <- lm(Crown_spread ~ CS_over_distance, data = LM_focal_tree_dataframe)
LM_lm_focal_CS_cooks <- cooks.distance(LM_lm_focal_CS) #calculating the cook.s D for each point
plot(LM_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_CS_cooks[(LM_lm_focal_CS_cooks > (3 * mean(LM_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_focal_tree_dataframe_no_CS_outliers <- LM_focal_tree_dataframe[-c(9, 13),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_CS <- gls(Canopy_long ~ CS_over_distance, data = LM_focal_tree_dataframe)
LM_gls_focal_CS_exp <- gls(Canopy_short ~ CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CS_gaus <- gls(Canopy_short ~ CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CS_spher <- gls(Canopy_short ~ CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CS_lin <- gls(Canopy_short ~ CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_CS_ratio <- gls(Canopy_short ~ CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CS <- model.sel(LM_gls_focal_CS, LM_gls_focal_CS_exp, LM_gls_focal_CS_gaus, LM_gls_focal_CS_spher, LM_gls_focal_CS_ratio) #without linear correlation LM_gls_focal_CS_lin
LM_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_focal_tree_dataframe, aes(x= LM_gls_focal_CS$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_focal_tree_dataframe, aes(sample = LM_gls_focal_CS$residuals))+
  geom_qq()

# shapiro-wilk, not sign so residuals are normally distributed
shapiro.test(LM_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LM_focal_tree_dataframe , aes(x = LM_gls_focal_CS$fitted, y = LM_gls_focal_CS$residuals))+
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
LM_tau_result_CS <- cor.test(LM_focal_tree_dataframe$Cs_over_distance, LM_focal_tree_dataframe$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_CS)

# Calculate the trend line
LM_trend_line_CS <- predict(loess(LM_focal_tree_dataframe$Crown_spread ~ LM_focal_tree_dataframe$CS_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_focal_tree_dataframe$CS_over_distance, y = (LM_focal_tree_dataframe$Crown_spread), color = "blue")) +
  geom_line(aes(x = LM_focal_tree_dataframe$CS_over_distance, y = LM_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#plotting the linear model in ggplot
ggplot(data = LM_focal_tree_dataframe, (aes(x=DBH_over_distance, y=DBH_ag)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("DBH over Distance")+
  ylab("DBH")

#Cook's D
LM_lm_focal_DBH <- lm(DBH_ag ~ DBH_over_distance, data = LM_focal_tree_dataframe)
LM_lm_focal_DBH_cooks <- cooks.distance(LM_lm_focal_DBH) #calculating the cook.s D for each point
plot(LM_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_DBH_cooks[(LM_lm_focal_DBH_cooks > (3 * mean(LM_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LM_focal_tree_dataframe_no_CS_outliers <- LM_focal_tree_dataframe[-c(9,13, 20),]
plot(LM_focal_tree_dataframe_no_CS_outliers$DBH_over_distance)

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_DBH <- gls(Canopy_short ~ DBH_over_distance, data = LM_focal_tree_dataframe)
LM_gls_focal_DBH_exp <- gls(Canopy_short ~ DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_DBH_gaus <- gls(Canopy_short ~ DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_DBH_spher <- gls(Canopy_short ~ DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_DBH_lin <- gls(Canopy_short ~ DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_focal_tree_dataframe)
LM_gls_focal_DBH_ratio <- gls(Canopy_short ~ DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_DHB <- model.sel(LM_gls_focal_DBH, LM_gls_focal_DBH_exp, LM_gls_focal_DBH_gaus, LM_gls_focal_DBH_spher, LM_gls_focal_DBH_ratio) #without linear control LM_gls_focal_DBH_lin
LM_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_focal_tree_dataframe, aes(x= LM_gls_focal_DBH$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_focal_tree_dataframe, aes(sample = LM_gls_focal_DBH$residuals))+
  geom_qq()

# shapiro-wilk, not sign so the residuals are normally distributed
shapiro.test(LM_gls_focal_DBH$residuals) 

#checking equal variance
ggplot(data = LM_focal_tree_dataframe , aes(x = LM_gls_focal_DBH$fitted, y = LM_gls_focal_DBH$residuals))+
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
LM_tau_result_DBH <- cor.test(LM_focal_tree_dataframe$DBH_over_distance, LM_focal_tree_dataframe$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LM_tau_result_DBH)

# Calculate the trend line
LM_trend_line_DBH <- predict(loess(LM_focal_tree_dataframe$DBH_ag ~ LM_focal_tree_dataframe$DBH_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LM_focal_tree_dataframe$DBH_over_distance, y = (LM_focal_tree_dataframe$DBH_ag), color = "blue")) +
  geom_line(aes(x = LM_focal_tree_dataframe$DBH_over_distance, y = LM_trend_line_DBH), color = "red") +
  labs(x = "DBH over Distance", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()


#### LC ####

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
LC_focal_results <- focal_function("LC")

#assigning necessary dataframes and objects 
LC_box_sf <- LC_focal_results[[1]] #bounding box
LC_box_sf_cropped <- LC_focal_results[[2]] #bounding box cropped by 20m
LC_fixed_field_data_processed_sf_cropped <- LC_focal_results[[3]] # cropped tree data
LC_tree_grid_cropped <- LC_focal_results[[4]] #grid with 40*mean population DBH as grid size
LC_focal_tree_buffers <- LC_focal_results[[5]] #focal tree buffers
LC_focal_tree_dataframe_sf <- LC_focal_results[[6]] #focal tree data
LC_focal_tree_dataframe <- as.data.frame(LC_focal_tree_dataframe_sf) #focal tree data

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
  geom_sf(data= LC_focal_tree_dataframe_sf, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LC_tree_grid_cropped) +
  geom_sf(data=LC_focal_tree_buffers, aes(color = focal_tree_row))+
  geom_sf(data=LC_focal_tree_dataframe_sf, aes(color = focal_tree_row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = LC_focal_tree_buffers)+
  geom_sf(data = LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_focal_tree_dataframe_sf, color = 'blue')

#plotting the tree points and their competition metrics
ggplot()+
  geom_sf(data=LC_focal_tree_dataframe_sf, aes(color = SCA_over_distance))

#descriptive statistics

#histograms
ggplot(LC_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LC_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(LC_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(LC_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(LC_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
LC_field_data_focal_summarized_focal <- LC_focal_tree_dataframe %>%
  dplyr::select(SCA_over_distance, LCA_over_distance, CS_over_distance, CA_over_distance, DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(LC_field_data_focal_summarized_focal)


#creating the generalized linear effects model

#creating x and y columns of the UTM 12N 
LC_focal_tree_dataframe$X.1 <- st_coordinates(LC_focal_tree_dataframe_sf)[,1]
LC_focal_tree_dataframe$Y <- st_coordinates(LC_focal_tree_dataframe_sf)[,2]

#SCA

#plotting the linear model in ggplot for SCA
ggplot(data = LC_focal_tree_dataframe, (aes(x=SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")


#Cook's D
LC_lm_focal_SCA <- lm(Canopy_short ~ LCA_over_distance, data = LC_focal_tree_dataframe)
LC_lm_focal_SCA_cooks <- cooks.distance(LC_lm_focal_SCA) #calculating the cook.s D for each point
plot(LC_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_SCA_cooks[(LC_lm_focal_SCA_cooks > (3 * mean(LC_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_focal_tree_dataframe_no_SCA_outliers <- LC_focal_tree_dataframe[-c(19),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_SCA <- gls(Canopy_short ~ SCA_over_distance, data = LC_focal_tree_dataframe)
LC_gls_focal_SCA_exp <- gls(Canopy_short ~ SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_SCA_gaus <- gls(Canopy_short ~ SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_SCA_spher <- gls(Canopy_short ~ SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_SCA_lin <- gls(Canopy_short ~ SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_SCA_ratio <- gls(Canopy_short ~ SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test <- model.sel(LC_gls_focal_SCA, LC_gls_focal_SCA_exp, LC_gls_focal_SCA_gaus, LC_gls_focal_SCA_spher, LC_gls_focal_SCA_lin, LC_gls_focal_SCA_ratio)
LC_AIC_test

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_focal_tree_dataframe, aes(x= LC_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_focal_tree_dataframe, aes(sample = LC_gls_focal_SCA$residuals))+
  geom_qq()

#shapiro-wilk test, not sign so normal residuals
shapiro.test(LC_gls_focal_SCA$residuals) 

#checking equal variance
ggplot(data = LC_focal_tree_dataframe , aes(x = LC_gls_focal_SCA$fitted, y = LC_gls_focal_SCA$residuals))+
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
LC_tau_result_SCA <- cor.test(LC_focal_tree_dataframe$SCA_over_distance, 
                              LC_focal_tree_dataframe$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_SCA)

# Calculate the trend line
LC_trend_line_SCA <- predict(loess(LC_focal_tree_dataframe$Canopy_short ~ LC_focal_tree_dataframe$SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_focal_tree_dataframe$SCA_over_distance, y = (LC_focal_tree_dataframe$Canopy_short), color = "blue")) +
  geom_line(aes(x = LC_focal_tree_dataframe$SCA_over_distance, y = LC_trend_line_SCA), color = "red") +
  labs(x = "SCA over Distance", y = "Short Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#LCA

#plotting the linear model in ggplot for SCA
ggplot(data = LC_focal_tree_dataframe, (aes(x=LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("LCA over Distance")+
  ylab("Long Canopy Axis")

#Cook's D
LC_lm_focal_LCA <- lm(Canopy_long ~ LCA_over_distance, data = LC_focal_tree_dataframe)
LC_lm_focal_LCA_cooks <- cooks.distance(LC_lm_focal_LCA) #calculating the cook.s D for each point
plot(LC_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_LCA_cooks[(LC_lm_focal_LCA_cooks > (3 * mean(LC_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_focal_tree_dataframe_no_LCA_outliers <- LC_focal_tree_dataframe[-c(19),]
plot(LC_focal_tree_dataframe_no_LCA_outliers$LCA_over_distance)

#unlogged version of generalized linear model
LC_gls_focal_LCA <- gls(Canopy_long ~ LCA_over_distance, data = LC_focal_tree_dataframe)
LC_gls_focal_LCA_exp <- gls(Canopy_long ~ LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_LCA_gaus <- gls(Canopy_long ~ LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_LCA_spher <- gls(Canopy_long ~ LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_LCA_lin <- gls(Canopy_long ~ LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_LCA_ratio <- gls(Canopy_long ~ LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_LCA <- model.sel(LC_gls_focal_LCA, LC_gls_focal_LCA_exp, LC_gls_focal_LCA_gaus, LC_gls_focal_LCA_spher) #LC_gls_focal_LCA_ratio
LC_AIC_test_LCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_focal_tree_dataframe, aes(x= LC_gls_focal_LCA_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance",
       subtitle = "Using Gaussian Control for Spatial Autocorrelation")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_focal_tree_dataframe, aes(sample = LC_gls_focal_LCA_gaus$residuals))+
  geom_qq()

#shapiro-wilk test, significant so non-normal residuals
shapiro.test(LC_gls_focal_LCA_gaus$residuals) 

#checking equal variance
ggplot(data = LC_focal_tree_dataframe, aes(x = LC_gls_focal_LCA_gaus$fitted, y = LC_gls_focal_LCA_gaus$residuals))+
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
LC_tau_result_LCA <- cor.test(LC_focal_tree_dataframe$LCA_over_distance, 
                              LC_focal_tree_dataframe$Canopy_long,  
                              method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_LCA)

# Calculate the trend line
LC_trend_line_LCA <- predict(loess(LC_focal_tree_dataframe$Canopy_long ~ LC_focal_tree_dataframe$LCA_over_distance))

# Extract fitted values from the GLS model
fitted_canopy <- fitted(LC_gls_focal_LCA_gaus)

# Create the data frame for plotting
line_df <- data.frame(
  LCA_over_distance = LC_focal_tree_dataframe$LCA_over_distance,
  fitted_canopy = fitted_canopy
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_focal_tree_dataframe$LCA_over_distance, y = (LC_focal_tree_dataframe$Canopy_long), color = "blue")) +
  geom_line(data = line_df, aes(x = LCA_over_distance, y = fitted_canopy), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#CA

#plotting the linear model in ggplot
ggplot(data = LC_focal_tree_dataframe, (aes(x=CA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
LC_lm_focal_CA <- lm(Canopy_area ~ CA_over_distance, data = LC_focal_tree_dataframe)
LC_lm_focal_CA_cooks <- cooks.distance(LC_lm_focal_LCA) #calculating the cook.s D for each point
plot(LC_lm_focal_CA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_CA_cooks[(LC_lm_focal_CA_cooks > (3 * mean(LC_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_focal_tree_dataframe_no_CA_outliers <- LC_focal_tree_dataframe[-c(19),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_CA <- gls(Canopy_area ~ CA_over_distance, data = LC_focal_tree_dataframe)
LC_gls_focal_CA_exp <- gls(Canopy_area ~ CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CA_gaus <- gls(Canopy_area ~ CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CA_spher <- gls(Canopy_area ~ CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CA_lin <- gls(Canopy_area ~ CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CA_ratio <- gls(Canopy_area ~ CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CA <- model.sel(LC_gls_focal_CA, LC_gls_focal_CA_exp, LC_gls_focal_CA_gaus, LC_gls_focal_CA_spher, LC_gls_focal_CA_lin, LC_gls_focal_CA_ratio)
LC_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_focal_tree_dataframe, aes(x= LC_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_focal_tree_dataframe, aes(sample = LC_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, sign for when residuals so we are using a Kendall's Tau Correlation Test non-parametric data
shapiro.test(LC_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LC_focal_tree_dataframe , aes(x = LC_gls_focal_CA$fitted, y = LC_gls_focal_CA$residuals))+
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
LC_tau_result_CA <- cor.test(LC_focal_tree_dataframe$CA_over_distance, 
                             LC_focal_tree_dataframe$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CA)

# Calculate the trend line
LC_trend_line_CA <- predict(loess(LC_focal_tree_dataframe$Canopy_area ~ LC_focal_tree_dataframe$CA_over_distance))

# Create a trend line plot

# Extract fitted values from the GLS model
fitted_canopy <- fitted(LC_gls_focal_CA)

# Create the data frame for plotting
line_df <- data.frame(
  CA_over_distance = LC_focal_tree_dataframe$CA_over_distance,
  fitted_canopy = fitted_canopy
)

#plotting
ggplot() +
  geom_point(aes(x = LC_focal_tree_dataframe$CA_over_distance, y = (LC_focal_tree_dataframe$Canopy_area), color = "red")) +
  geom_line(aes(x = LC_focal_tree_dataframe$CA_over_distance, y = fitted_canopy), color = "red") +
  labs(x = "CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot for SCA
ggplot(data = LC_focal_tree_dataframe, (aes(x=CS_over_distance, y=Crown_spread)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CS over Distance")+
  ylab("Crown Spread")

#Cook's D
LC_lm_focal_CS <- lm(Crown_spread ~ CS_over_distance, data = LC_focal_tree_dataframe)
LC_lm_focal_CS_cooks <- cooks.distance(LC_lm_focal_CS) #calculating the cook.s D for each point
plot(LC_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_CS_cooks[(LC_lm_focal_CS_cooks > (3 * mean(LC_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_focal_tree_dataframe_no_CS_outliers <- LC_focal_tree_dataframe[-c(19),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_CS <- gls(Crown_spread ~ CS_over_distance, data = LC_focal_tree_dataframe)
LC_gls_focal_CS_exp <- gls(Crown_spread ~ CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CS_gaus <- gls(Crown_spread ~ CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CS_spher <- gls(Crown_spread ~ CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CS_lin <- gls(Crown_spread ~ CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_CS_ratio <- gls(Crown_spread ~ CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CS <- model.sel(LC_gls_focal_CS, LC_gls_focal_CS_exp, LC_gls_focal_CS_gaus, LC_gls_focal_CS_ratio) #without linear correlation and without spherical
LC_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_focal_tree_dataframe, aes(x= LC_gls_focal_CS_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq norm
ggplot(LC_focal_tree_dataframe, aes(sample = LC_gls_focal_CS_gaus$residuals))+
  geom_qq()

# shapiro-wilk, n sign for both versions with and without outliers so used Kendall's Tau non-parametric test
shapiro.test(LC_gls_focal_CA$residuals) # shapiro-wilk, n sign for both versions with and without outliers so used Kendall's Tau Test non-parametric test

#checking equal variance
ggplot(data = LC_focal_tree_dataframe , aes(x = LC_gls_focal_CS_gaus$fitted, y = LC_gls_focal_CS_gaus$residuals))+
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
LC_tau_result_CS <- cor.test(LC_focal_tree_dataframe$CS_over_distance, 
                             LC_focal_tree_dataframe$Crown_spread,  
                             method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CS)

# Calculate the trend line
LC_trend_line_CS <- predict(loess(LC_focal_tree_dataframe$Crown_spread ~ LC_focal_tree_dataframe$CS_over_distance))

# Extract fitted values from the GLS model
fitted_crown <- fitted(LC_gls_focal_CS_gaus)

# Create the data frame for plotting
line_df <- data.frame(
  CA_over_distance = LC_focal_tree_dataframe$CS_over_distance,
  fitted_crown = fitted_crown
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_focal_tree_dataframe$CS_over_distance, y = (LC_focal_tree_dataframe$Crown_spread), color = "red")) +
  geom_line(data = line_df, aes(x = CS_over_distance, y = fitted_crown), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#plotting the linear model in ggplot 
ggplot(data = LC_focal_tree_dataframe, (aes(x=DBH_over_distance, y=DBH_ag)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("DBH over Distance")+
  ylab("DBH")

#Cook's D
LC_lm_focal_DBH <- lm(DBH_ag ~ DBH_over_distance, data = LC_focal_tree_dataframe)
LC_lm_focal_DBH_cooks <- cooks.distance(LC_lm_focal_DBH) #calculating the cook.s D for each point
plot(LC_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LC_lm_focal_DBH_cooks[(LC_lm_focal_DBH_cooks > (3 * mean(LC_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
LC_focal_tree_dataframe_no_DBH_outliers <- LC_focal_tree_dataframe[-c(16, 19),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_DBH <- gls(DBH_ag ~ DBH_over_distance, data = LC_focal_tree_dataframe)
LC_gls_focal_DBH_exp <- gls(DBH_ag ~ DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_DBH_gaus <- gls(DBH_ag ~ DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_DBH_spher <- gls(DBH_ag ~ DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_DBH_lin <- gls(DBH_ag ~ DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_focal_tree_dataframe)
LC_gls_focal_DBH_ratio <- gls(DBH_ag ~ DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_DHB <- model.sel(LC_gls_focal_DBH, LC_gls_focal_DBH_exp, LC_gls_focal_DBH_gaus, LC_gls_focal_DBH_spher, LC_gls_focal_DBH_ratio) #without linear correlation
LC_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_focal_tree_dataframe, aes(x= LC_gls_focal_DBH_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plot
ggplot(LC_focal_tree_dataframe, aes(sample = LC_gls_focal_DBH_gaus$residuals))+
  geom_qq()

# shapiro-wilk, not significant so normal residuals
shapiro.test(LC_gls_focal_DBH_gaus$residuals)

#checking equal variance
ggplot(data = LC_focal_tree_dataframe, aes(x = LC_gls_focal_DBH$fitted, y = LC_gls_focal_DBH$residuals))+
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
LC_tau_result_DBH <- cor.test(LC_focal_tree_dataframe$CS_over_distance, 
                              LC_focal_tree_dataframe$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_DBH)

# Calculate the trend line
LC_trend_line_DBH <- predict(loess(LC_focal_tree_dataframe$DBH_ag ~ LC_focal_tree_dataframe$DBH_over_distance))

# Extract fitted values from the GLS model
fitted_DBH <- fitted(LC_gls_focal_DBH_gaus)

# Create the data frame for plotting
LC_line_df_DBH <- data.frame(
  DBH_over_distance = LC_focal_tree_dataframe$DBH_over_distance,
  fitted_DBH = fitted_DBH
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_focal_tree_dataframe$CS_over_distance, y = (LC_focal_tree_dataframe$DBH_ag), color = "red")) +
  geom_line(data = LC_line_df_DBH, aes(x = DBH_over_distance, y = fitted_DBH), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


#### SD ####

#running the function to determine the focal trees, neighbors, and calculate the competition metrics for each focal tree
SD_focal_results <- focal_function("SD")

#assigning necessary dataframes and objects 
SD_box_sf <- SD_focal_results[[1]] #bounding box
SD_box_sf_cropped <- SD_focal_results[[2]] #bounding box cropped by 20m
SD_fixed_field_data_processed_sf_cropped <- SD_focal_results[[3]] # cropped tree data
SD_tree_grid_cropped <- SD_focal_results[[4]] #grid with 40*mean population DBH as grid size
SD_focal_tree_buffers <- SD_focal_results[[5]] #focal tree buffers
SD_focal_tree_dataframe_sf <- SD_focal_results[[6]] #focal tree data as a spatial object
SD_focal_tree_dataframe <- as.data.frame(SD_focal_tree_dataframe_sf) #focal tree data

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
  geom_sf(data= SD_focal_tree_dataframe_sf, aes(color = X))

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data=SD_focal_tree_buffers, aes(color = focal_tree_row))+
  geom_sf(data=SD_focal_tree_dataframe_sf, aes(color = focal_tree_row))

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = SD_focal_tree_buffers)+
  geom_sf(data = SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_focal_tree_dataframe_sf, color = 'blue')

#plotting the tree points and their competition metrics
ggplot()+
  geom_sf(data=SD_focal_tree_dataframe_sf, aes(color = SCA_over_distance))


#descriptive statistics

#histograms
ggplot(SD_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = SCA_over_distance))+
  xlab("Sum of Short Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(SD_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = LCA_over_distance))+
  xlab("Sum of Long Canopy Axis over Distance")+
  ylab("Frequency")

ggplot(SD_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CS_over_distance))+
  xlab("Sum of Canopy Spread over Distance")+
  ylab("Frequency")

ggplot(SD_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = CA_over_distance))+
  xlab("Sum of Canopy Area over Distance")+
  ylab("Frequency")

ggplot(SD_focal_tree_dataframe) + # Generate the base plot
  geom_histogram(aes(x = DBH_over_distance))+
  xlab("Sum of Aggregated DBH over Distance")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
SD_field_data_focal_summarized_focal <- SD_focal_tree_dataframe %>%
  dplyr::select(SCA_over_distance, LCA_over_distance, CS_over_distance, CA_over_distance, DBH_over_distance) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(SD_field_data_focal_summarized_focal)


#creating x and y columns of the UTM 12N 
SD_focal_tree_dataframe$X.1 <- st_coordinates(SD_focal_tree_dataframe_sf)[,1]
SD_focal_tree_dataframe$Y <- st_coordinates(SD_focal_tree_dataframe_sf)[,2]

#SCA

#plotting the linear model in ggplot 
ggplot(data = SD_focal_tree_dataframe, (aes(x=SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")

#Cook's D
SD_lm_focal_SCA <- lm(Canopy_short ~ SCA_over_distance, data = SD_focal_tree_dataframe)
SD_lm_focal_SCA_cooks <- cooks.distance(SD_lm_focal_SCA) #calculating the cook.s D for each point
plot(SD_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_SCA_cooks[(SD_lm_focal_SCA_cooks > (3 * mean(SD_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_focal_tree_dataframe_no_SCA_outliers <- SD_focal_tree_dataframe[-c(23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_SCA <- gls(Canopy_short ~ SCA_over_distance, data = SD_focal_tree_dataframe)
SD_gls_focal_SCA_exp <- gls(Canopy_short ~ SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_SCA_gaus <- gls(Canopy_short ~ SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_SCA_spher <- gls(Canopy_short ~ SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_SCA_lin <- gls(Canopy_short ~ SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_SCA_ratio <- gls(Canopy_short ~ SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_SCA <- model.sel(SD_gls_focal_SCA, SD_gls_focal_SCA_lin, SD_gls_focal_SCA_exp, SD_gls_focal_SCA_gaus, SD_gls_focal_SCA_spher, SD_gls_focal_SCA_ratio) 
SD_AIC_test_SCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_focal_tree_dataframe, aes(x= SD_gls_focal_SCA_exp$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq nrom
ggplot(SD_focal_tree_dataframe, aes(sample = SD_gls_focal_SCA_exp$residuals))+
  geom_qq()

#shapiro-welk test, not significant so normal residuals
shapiro.test(SD_gls_focal_SCA_exp$residuals)

#checking equal variance
ggplot(data = SD_focal_tree_dataframe , aes(x = SD_gls_focal_SCA_exp$fitted, y = SD_gls_focal_SCA_exp$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_SCA_exp, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_SCA_exp)

#non parametric Kendall's Tau Test
SD_tau_result_SCA <- cor.test(SD_focal_tree_dataframe$SCA_over_distance, 
                              SD_focal_tree_dataframe_no_SCA_outliers$Canopy_short,  
                              method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_SCA)

# Calculate the trend line
SD_trend_line_SCA <- predict(loess(SD_focal_tree_dataframe$Canopy_short ~ SD_focal_tree_dataframe$SCA_over_distance))

# Extract fitted values from the GLS model
fitted_SCA <- fitted(SD_gls_focal_SCA)

# Create the data frame for plotting
SD_line_df_SCA <- data.frame(
  SCA_over_distance = SD_focal_tree_dataframe$SCA_over_distance,
  fitted_SCA = fitted_SCA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_focal_tree_dataframe$SCA_over_distance, y = (SD_focal_tree_dataframe$Canopy_short), color = "blue")) +
  geom_line(data = SD_line_df_SCA, aes(x = SCA_over_distance, y = fitted_SCA), color = "red") +
  labs(x = "Sum of SCA over Distance", y = "Short Canopy Axis ", title = "Trend Line Plot") +
  theme_minimal()


#LCA

#plotting the linear model in ggplot for LCA
ggplot(data = SD_focal_tree_dataframe, (aes(x=LCA_over_distance, y=Canopy_long)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Long Canopy Axis")

#Cook's D
SD_lm_focal_LCA <- lm(Canopy_long ~ LCA_over_distance, data = SD_focal_tree_dataframe)
SD_lm_focal_LCA_cooks <- cooks.distance(SD_lm_focal_LCA) #calculating the cook.s D for each point
plot(SD_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_LCA_cooks[(SD_lm_focal_LCA_cooks > (3 * mean(SD_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_focal_tree_dataframe_no_LCA_outliers <- SD_focal_tree_dataframe[-c(23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_LCA <- gls(Canopy_long ~ LCA_over_distance, data = SD_focal_tree_dataframe)
SD_gls_focal_LCA_exp <- gls(Canopy_long ~ LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_LCA_gaus <- gls(Canopy_long ~ LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_LCA_spher <- gls(Canopy_long ~ LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_LCA_lin <- gls(Canopy_long ~ LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_LCA_ratio <- gls(Canopy_long ~ LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_LCA <- model.sel(SD_gls_focal_LCA, SD_gls_focal_LCA_exp, SD_gls_focal_LCA_gaus, SD_gls_focal_LCA_spher, SD_gls_focal_LCA_ratio) #without the linear correction model
SD_AIC_test_LCA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_focal_tree_dataframe, aes(x= SD_gls_focal_LCA_spher$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_focal_tree_dataframe, aes(sample = SD_gls_focal_LCA_spher$residuals))+
  geom_qq()

#shapiro-wilk test, not sign so normal residuals
shapiro.test(SD_gls_focal_LCA_spher$residuals) 

#checking equal variance
ggplot(data = SD_focal_tree_dataframe , aes(x = SD_gls_focal_LCA_spher$fitted, y = SD_gls_focal_LCA_spher$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and LCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_LCA_spher, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_LCA_spher)

#non parametric Kendall's Tau Test
SD_tau_result_LCA <- cor.test(SD_focal_tree_dataframe_no_LCA_outliers$LCA_over_distance, SD_focal_tree_dataframe_no_LCA_outliers$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_LCA)

# Calculate the trend line
SD_trend_line_LCA <- predict(loess(SD_focal_tree_dataframe$Canopy_long ~ SD_focal_tree_dataframe$LCA_over_distance))

# Extract fitted values from the GLS model
fitted_LCA <- fitted(SD_gls_focal_LCA)

# Create the data frame for plotting
SD_line_df_LCA <- data.frame(
  LCA_over_distance = SD_focal_tree_dataframe$LCA_over_distance,
  fitted_LCA = fitted_LCA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_focal_tree_dataframe$LCA_over_distance, y = (SD_focal_tree_dataframe$Canopy_long), color = "blue")) +
  geom_line(data = SD_line_df_LCA, aes(x = LCA_over_distance, y = fitted_LCA), color = "red") +
  labs(x = "Sum of LCA over Distance", y = "Long Canopy Axis ", title = "Trend Line Plot") +
  theme_minimal()

#CA

#plotting the linear model in ggplot 
ggplot(data = SD_focal_tree_dataframe, (aes(x=CA_over_distance, y=Canopy_area)))+ 
  geom_smooth(method='glm')+
  geom_point()+
  xlab("CA over Distance")+
  ylab("Canopy Area")

#Cook's D
SD_lm_focal_CA <- lm(Canopy_area ~ CA_over_distance, data = SD_focal_tree_dataframe)
SD_lm_focal_CA_cooks <- cooks.distance(SD_lm_focal_CA) #calculating the cook.s D for each point
plot(SD_lm_focal_LCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_CA_cooks[(SD_lm_focal_CA_cooks > (3 * mean(SD_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_focal_tree_dataframe_no_CA_outliers <- SD_focal_tree_dataframe[-c(23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CA <- gls(Canopy_area ~ CA_over_distance, data = SD_focal_tree_dataframe)
SD_gls_focal_CA_exp <- gls(Canopy_area ~ CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CA_gaus <- gls(Canopy_area ~ CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CA_spher <- gls(Canopy_area ~ CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CA_lin <- gls(Canopy_area ~ CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CA_ratio <- gls(Canopy_area ~ CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CA <- model.sel(SD_gls_focal_CA, SD_gls_focal_CA_exp, SD_gls_focal_CA_gaus, SD_gls_focal_CA_lin, SD_gls_focal_CA_ratio) #SD_gls_focal_CA_spher
SD_AIC_test_CA

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_focal_tree_dataframe, aes(x= SD_gls_focal_CA_lin$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_focal_tree_dataframe, aes(sample = SD_gls_focal_CA_lin$residuals))+
  geom_qq()

# shapiro-wilk, significant so residuals non-normal
shapiro.test(SD_gls_focal_CA_lin$residuals) 

#checking equal variance
ggplot(data = SD_focal_tree_dataframe, aes(x = SD_gls_focal_CA_lin$fitted, y = SD_gls_focal_CA_lin$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_CA_lin, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_CA_lin)

#non parametric Kendall's Tau Test
SD_tau_result_CA <- cor.test(SD_focal_tree_dataframe$CA_over_distance, 
                             SD_focal_tree_dataframe$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CA)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_focal_tree_dataframe$Canopy_area ~ SD_focal_tree_dataframe$CA_over_distance))

# Extract fitted values from the GLS model
fitted_CA <- fitted(SD_gls_focal_CA)

# Create the data frame for plotting
SD_line_df_CA <- data.frame(
  CA_over_distance = SD_focal_tree_dataframe$CA_over_distance,
  fitted_CA = fitted_CA
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_focal_tree_dataframe$CA_over_distance, y = (SD_focal_tree_dataframe$Canopy_area), color = "red")) +
  geom_line(data = SD_line_df_CA, aes(x = CA_over_distance, y = fitted_CA), color = "red") +
  labs(x = "Sum of CA over Distance", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#plotting the linear model in ggplot 
ggplot(data = SD_focal_tree_dataframe, (aes(x=CS_over_distance, y=Crown_spread)))+ 
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
SD_lm_focal_CS <- lm(Crown_spread ~ CS_over_distance, data = SD_focal_tree_dataframe)
SD_lm_focal_CS_cooks <- cooks.distance(SD_lm_focal_CS) #calculating the cook.s D for each point
plot(SD_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_CS_cooks[(SD_lm_focal_CS_cooks > (3 * mean(SD_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_focal_tree_dataframe_no_CS_outliers <- SD_focal_tree_dataframe[-c(23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CS <- gls(Crown_spread ~ CS_over_distance, data = SD_focal_tree_dataframe)
SD_gls_focal_CS_exp <- gls(Crown_spread ~ CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CS_gaus <- gls(Crown_spread ~ CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CS_spher <- gls(Crown_spread ~ CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CS_lin <- gls(Crown_spread ~ CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_CS_ratio <- gls(Crown_spread ~ CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CS <- model.sel(SD_gls_focal_CS, SD_gls_focal_CS_exp, SD_gls_focal_CS_gaus, SD_gls_focal_CS_spher, SD_gls_focal_CS_ratio) #without linear correlation
SD_AIC_test_CS

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_focal_tree_dataframe, aes(x= SD_gls_focal_CS_exp$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_focal_tree_dataframe, aes(sample = SD_gls_focal_CS_exp$residuals))+
  geom_qq()

# shapiro-wilk, not signficant, meaning not significantly different from normal
shapiro.test(SD_gls_focal_CS_exp$residuals) 

#checking equal variance
ggplot(data = SD_focal_tree_dataframe , aes(x = SD_gls_focal_CS_exp$fitted, y = SD_gls_focal_CS_exp$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and CS over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_CS_exp, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_CS_exp)

#non parametric Kendall's Tau Test
SD_tau_result_CS <- cor.test(SD_focal_tree_dataframe$CS_over_distance, SD_focal_tree_dataframe$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CS)

# Calculate the trend line
SD_trend_line_CS <- predict(loess(SD_focal_tree_dataframe$Crown_spread ~ SD_focal_tree_dataframe$CS_over_distance))

# Extract fitted values from the GLS model
fitted_CS <- fitted(SD_gls_focal_CS)

# Create the data frame for plotting
SD_line_df_CS <- data.frame(
  CS_over_distance = SD_focal_tree_dataframe$CS_over_distance,
  fitted_CS = fitted_CS
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_focal_tree_dataframe$CS_over_distance, y = (SD_focal_tree_dataframe$Crown_spread), color = "red")) +
  geom_line(data = SD_line_df_CS, aes(x = CS_over_distance, y = fitted_CS), color = "red") +
  labs(x = "Sum of CS over Distance", y = "Crown Spread", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#plotting the linear model in ggplot 
ggplot(data = SD_focal_tree_dataframe, (aes(x=DBH_over_distance, y=DBH_ag)))+ 
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
SD_lm_focal_DBH <- lm(DBH_ag ~ DBH_over_distance, data = SD_focal_tree_dataframe)
SD_lm_focal_DBH_cooks <- cooks.distance(SD_lm_focal_DBH) #calculating the cook.s D for each point
plot(SD_lm_focal_DBH_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_DBH_cooks[(SD_lm_focal_DBH_cooks > (3 * mean(SD_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_focal_tree_dataframe_no_DBH_outliers <- SD_focal_tree_dataframe[-c(23),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_DBH <- gls(DBH_ag ~ DBH_over_distance, data = SD_focal_tree_dataframe)
SD_gls_focal_DBH_exp <- gls(DBH_ag ~ DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_DBH_gaus <- gls(DBH_ag ~ DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_DBH_spher <- gls(DBH_ag ~ DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_DBH_lin <- gls(DBH_ag ~ DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_focal_tree_dataframe)
SD_gls_focal_DBH_ratio <- gls(DBH_ag ~ DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_focal_tree_dataframe)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_DHB <- model.sel(SD_gls_focal_DBH, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_exp, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_gaus, SD_gls_focal_DBH_spher, SD_gls_focal_DBH_ratio) 
SD_AIC_test_DHB

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_focal_tree_dataframe, aes(x= SD_gls_focal_DBH_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_focal_tree_dataframe, aes(sample = SD_gls_focal_DBH_gaus$residuals))+
  geom_qq()

# shapiro-wilk, not significant so normal
shapiro.test(SD_gls_focal_DBH_gaus$residuals) 

#checking equal variance
ggplot(data = SD_focal_tree_dataframe , aes(x = SD_gls_focal_DBH_gaus$fitted, y = SD_gls_focal_DBH_gaus$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and DBH over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_DBH_gaus, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_DBH_gaus)

#non parametric Kendall's Tau Test
SD_tau_result_DBH <- cor.test(SD_focal_tree_dataframe$DBH_over_distance, SD_focal_tree_dataframe$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_DBH)

# Calculate the trend line
SD_trend_line_DBH <- predict(loess(SD_focal_tree_dataframe$DBH_ag ~ SD_focal_tree_dataframe$DBH_over_distance))

# Extract fitted values from the GLS model
fitted_DBH <- fitted(SD_gls_focal_DBH)

# Create the data frame for plotting
SD_line_df_DBH <- data.frame(
  DBH_over_distance = SD_focal_tree_dataframe$DBH_over_distance,
  fitted_DBH = fitted_DBH
)

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_focal_tree_dataframe$DBH_over_distance, y = (SD_focal_tree_dataframe$DBH_ag), color = "red")) +
  geom_line(data = SD_line_df_DBH, aes(x = DBH_over_distance, y = fitted_DBH), color = "red") +
  labs(x = "Sum of DBH over Distance", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

