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
library(raster) #to use point distance
library(lme4) #to use the linear mixed effects model function


fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#transforming the shapefile of trees from WGS84 into equal area projection UTM 12N
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

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

#### Computing Average Nearest Neighbors for each tree ####

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
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5))) # %>% #creates a column of the average distances (1-5) of each individual
#dplyr::select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances



#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
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

###MAKING THE LINEAR MODEL###

#Creating a grid over the tree points 
LM_tree_grid <- st_make_grid(LM_fixed_field_data_processed_sf, cellsize = (((40*mean(LM_fixed_field_data_processed$DBH_ag))*2)*2))

#ASH NOTE: BELOW SEEMS TO BE RETURNING NON UNIQUE GEOMETRIES, PERHAPS BECAUSE EDGE CASES ARE COUNTED AS UNIQUE OBSERVATIONS --> NEED TO USE SOMETHING LIKE DISTINCT OR UNIQUE ON THE GEOMETRIES TO MAKE SURE THIS GETS FIXED
# FOUND THIS OUT BY PLOTTING FOCAL POINTS AND BUFFERS AGAINST THE GRID ITSELF

#create a tibble with the the number of trees within the grids that contain trees
LM_tree_grid_points_within <- st_contains(LM_tree_grid, LM_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 0) #filter out any grids without trees in them

#filter out the grids to only have the grid cells that contain trees
LM_tree_grid_inside <- LM_tree_grid %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% LM_tree_grid_points_within$row) #only keep polygons that match the row number of the grid cells with trees within them 

View(LM_tree_grid_inside)

#creating the buffer around all of the tree points
LM_tree_buffers <-st_buffer(LM_fixed_field_data_processed$geometry, 40*mean(LM_fixed_field_data_processed$DBH_ag))


#selecting the random points from each grid cell to be a focal point
set.seed(25) #set a seed for the sample random generator
focal_pts <- c() #create an empty vector of focal trees
grid_number <- c()
for (i in 1:nrow(LM_tree_grid_inside)){ #for the length of the grids with trees inside of them
  LM_tree_grid_inside_df <- as.data.frame(LM_tree_grid_inside)
  LM_tree_grid_inside_df_i <- LM_tree_grid_inside_df[i,] #isolate a row of the grid dataframe
  LM_tree_grid_inside_sf_i <- st_as_sf(LM_tree_grid_inside_df_i) #set the row as a simple feature
  all_pts <- st_covered_by(LM_tree_buffers, LM_tree_grid_inside_sf_i)
  #all_pts <- st_contains(LM_tree_grid_inside_sf_i, LM_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts <- all_pts %>%
    as.data.frame()
  if (nrow(possible_pts) == 0){
    focal_pts <- focal_pts
  } else {
    focal_pt <- sample(possible_pts$row.id, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
    focal_pts <- c(focal_pts, focal_pt)
    grid_number <- c(grid_number, i)
  }
}

View(all_pts)

#filtering out point data to be just the focal points
LM_fixed_field_data_processed_focal <- LM_fixed_field_data_processed %>%
  filter(X %in% focal_pts) 

#creating the buffer around the focal points
LM_focal_tree_buffers <-st_buffer(LM_fixed_field_data_processed_focal$geometry, 40*mean(LM_fixed_field_data_processed_focal$DBH_ag))


#EXPERIMENT
LM_tree_grid_inside_sf_i_contains_experiment <- st_contains(LM_tree_grid_inside_sf_i, LM_fixed_field_data_processed_focal, sparse = F)

#graphing the selected focal trees, the buffers, the grid
ggplot()+
  geom_sf(data = LM_tree_grid)+
  geom_sf(data = LM_tree_grid_inside, fill = "DodgerBlue") +
  # geom_sf(data = LM_tree_buffers)
  # geom_sf(data=possible_pt_buffer, color = 'red')+
  geom_sf(data=LM_focal_tree_buffers, color = "blue")+
  geom_sf(data= LM_fixed_field_data_processed_focal, aes(color = X))
# geom_sf(data = possible_pts_buffers.df$geometry, color = "red")
View(LM_tree_grid_inside)



ggplot()+
  geom_sf(data=river_LM_trans)+
  geom_sf(data=LM_focal_tree_buffers)+
  geom_sf(data=LM_fixed_field_data_processed_sf)+
  geom_sf(data=LM_fixed_field_data_processed_focal$geometry, fill ="red")

#calculating the size/distance for focal trees and neighbors within buffers for buffers with only the focal tree and with more 

#create a tibble with the the number of trees within the buffers that contain trees
LM_tree_buffers_points_within_0 <- st_contains(LM_focal_tree_buffers, LM_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 0) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
LM_tree_buffer_inside_0 <- LM_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% LM_tree_buffers_points_within_0$row) #only keep polygons that match the row number of the grid cells with trees within them 

#calculating the size/distance for focal trees and neighbors within buffers for buffers with more than just the focal tree

#create a tibble with the the number of trees within the buffers that contain trees
LM_tree_buffers_points_within <- st_contains(LM_focal_tree_buffers, LM_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 1) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
LM_tree_buffer_inside <- LM_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% LM_tree_buffers_points_within$row) #only keep buffers that match the row number of the buffers cells with trees within them 
View(LM_tree_buffer_inside)
#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"

#Checking that row number in focal dataset is the same as the buffer dataset
LM_fixed_field_data_processed_focal_row <- LM_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()

View(LM_fixed_field_data_processed_focal_row)

LM_tree_buffer_inside_0 <- mutate(LM_tree_buffer_inside_0, row = as.factor(row))

ggplot()+
  geom_sf(data = LM_tree_grid) +
  #geom_sf(data=river_LM_trans)+
  geom_sf(data=LM_tree_buffer_inside_0, aes(color = row))+
  #geom_sf(data=LM_tree_buffer_inside, fill = "red")+
  #geom_sf(data=LM_fixed_field_data_processed_sf)
  geom_sf(data=LM_fixed_field_data_processed_focal_row, aes(color = row))


#create a dataframe with the column for all of the focal distances
LM_fixed_field_data_processed_focal_dist <- LM_fixed_field_data_processed %>%
  add_column(focal_distance = NA) #add a column for distances of neighbors to focal tree

LM_fixed_field_data_all_focal_trees <- tibble()#creating the empty tibble 

#calculating the distances of each tree within the buffer to the focal tree and the competition metric values
for (i in 1:nrow(LM_fixed_field_data_processed_focal)){ #for the length of the buffers with trees inside of them
  row_num = i
  LM_tree_buffer_inside_df <- as.data.frame(LM_tree_buffer_inside)
  LM_tree_buffer_inside_df_i <- LM_tree_buffer_inside_df %>% 
    filter(row == row_num) #isolate a row of the buffer dataframe
  LM_tree_buffer_inside_sf_i <- st_as_sf(LM_tree_buffer_inside_df_i) #set the row as a simple feature
  all_pts_buffer <- st_contains(LM_tree_buffer_inside_sf_i, LM_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
  LM_fixed_field_data_processed_trees <- LM_fixed_field_data_processed %>%
    filter(X %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer
  LM_fixed_field_data_focal_tree <- LM_fixed_field_data_processed_focal %>%
    filter(X %in% LM_fixed_field_data_processed_trees$X) #create a dataframe with only the focal tree
  #START NEW FROM ASH
  LM_fixed_field_data_nieghbor_trees <- LM_fixed_field_data_processed_trees %>%
    filter(X %notin% LM_fixed_field_data_focal_tree$X) %>% #create a dataframe with only the neighbors of the focal tree
    mutate(focal_distance = as.numeric(st_distance(geometry,  LM_fixed_field_data_focal_tree$geometry, by_element = T))) %>% #calculate the distance between the focal tree and each tree that nieghbors it
    mutate(focal_distance = case_when(focal_distance == 0 ~ 0.0000016, 
                                      focal_distance != 0 ~ focal_distance)) %>% #replace values of 0 (if the coords are the same for multiple trees) with a value an order of magnitude smaller than the smallest distance in our dataset
    mutate(SCA_over_distance = Canopy_short/focal_distance) %>% #creating a column with the short canopy axis size value divided by the tree's distance from the focal tree
    mutate(LCA_over_distance = Canopy_long/focal_distance) %>%
    mutate(CA_over_distance = Canopy_area/focal_distance) %>%
    mutate(CS_over_distance = Crown_spread/focal_distance) %>%
    mutate(DBH_over_distance = DBH_ag/focal_distance)
  
  sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
  sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
  sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
  sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
  sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  
  
  for (y in 1:nrow(LM_fixed_field_data_nieghbor_trees)){ #adding the size values of each neighbor to a sum total of the neighbors size values
    sum_SCA_over_distance = sum_SCA_over_distance + LM_fixed_field_data_processed_neighbors$SCA_over_distance[y] #summing the SCA of each neighbor
    sum_LCA_over_distance = sum_LCA_over_distance + LM_fixed_field_data_processed_neighbors$LCA_over_distance[y] #summing the LCA of each neighbor
    sum_CA_over_distance = sum_CA_over_distance + LM_fixed_field_data_processed_neighbors$CA_over_distance[y] #summing the CA of each neighbor
    sum_CS_over_distance = sum_CS_over_distance + LM_fixed_field_data_processed_neighbors$CS_over_distance[y] #summing the CS of each neighbor
    sum_DBH_over_distance = sum_DBH_over_distance + LM_fixed_field_data_processed_neighbors$DBH_over_distance[y] #summing the DBH of each neighbor
  }
  
  #creating a tibble with all of the calculated sizes over distances 
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  LM_fixed_field_data_focal_tree <- cbind(LM_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  LM_fixed_field_data_all_focal_trees <- rbind(LM_fixed_field_data_all_focal_trees, LM_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
  
  
  
  #END NEW FROM ASH
  # for (x in 1:nrow(LM_fixed_field_data_processed_trees)){ ###BIG NOTE FROM ASH: YOU CAN'T USE THE SAME VALUE TO ITERATE ACROSS TWO NESTED FOR LOOPS #for each tree in each buffer, calculate the distance from the tree to the focal tree
  #   LM_fixed_field_data_processed_focal_dist$focal_distance[LM_fixed_field_data_processed_trees[x,]$X] <- crossdist(LM_fixed_field_data_processed_trees[x,]$X.1, LM_fixed_field_data_processed_trees[x,]$Y, 
  #                                                   LM_fixed_field_data_focal_tree$X.1, LM_fixed_field_data_focal_tree$Y)
}

#NOTE FROM ASH: YOU DON'T WANT A SINGLE DATASET WITH DISTS TO FOCAL TREES AS A SINGLE COLUMN BECAUSE SOME OF THESE ARE NEIGHBORS TO MULTIPLE TREES --> THIS ALLOWS THOSE TREES TO ONLY HAVE ONE VALUE

}

View(LM_fixed_field_data_processed_focal_dist)

#create columns with the size values divided by the distance to the focal tree values and turning infinite values into their regular values
LM_fixed_field_data_processed_focal_dist <- LM_fixed_field_data_processed_focal_dist %>%
  mutate(SCA_over_distance = Canopy_short/focal_distance) %>% #creating a column with the short canopy axis size value divided by the tree's distance from the focal tree
  mutate(SCA_over_distance = case_when(is.infinite(SCA_over_distance) ~ 0.0000016,
                                       !is.infinite(SCA_over_distance) ~ SCA_over_distance)) %>% #if the SCA_over_distance value is infinite because it is dividing by a zero, we set it to NA and if not, it remains the value it was before
  mutate(LCA_over_distance = Canopy_long/focal_distance) %>%
  mutate(LCA_over_distance = case_when(is.infinite(LCA_over_distance) ~ 0.0000016,
                                       !is.infinite(LCA_over_distance) ~ LCA_over_distance)) %>%
  mutate(CA_over_distance = Canopy_area/focal_distance) %>%
  mutate(CA_over_distance = case_when(is.infinite(CA_over_distance) ~ 0.0000016,
                                      !is.infinite(CA_over_distance) ~ CA_over_distance)) %>%
  mutate(CS_over_distance = Crown_spread/focal_distance) %>%
  mutate(CS_over_distance = case_when(is.infinite(CS_over_distance) ~ 0.0000016,
                                      !is.infinite(CS_over_distance) ~ CS_over_distance)) %>%
  mutate(DBH_over_distance = DBH_ag/focal_distance) %>%
  mutate(DBH_over_distance = case_when(is.infinite(DBH_over_distance) ~ 0.0000016,
                                       !is.infinite(DBH_over_distance) ~ DBH_over_distance))

View(LM_fixed_field_data_processed_focal_dist)

#creating a dataframe with  the sum of the neighbors sizes/their distances for each buffer/focal tree 

LM_fixed_field_data_all_focal_trees <- tibble()#creating the empty tibble 

for (i in 1:nrow(LM_tree_buffer_inside)){  #for the length of the buffers with trees inside of them
  LM_tree_buffer_inside_df <- as.data.frame(LM_tree_buffer_inside)
  LM_tree_buffer_inside_df_i <- LM_tree_buffer_inside_df[i,] #isolate a row of the buffer dataframe
  LM_tree_buffer_inside_sf_i <- st_as_sf(LM_tree_buffer_inside_df_i) #set the row as a simple feature
  all_pts_buffer <- st_contains(LM_tree_buffer_inside_sf_i, LM_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
  LM_fixed_field_data_processed_trees <- LM_fixed_field_data_processed_focal_dist %>%
    filter(X %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer
  LM_fixed_field_data_focal_tree <- LM_fixed_field_data_processed_focal %>%
    filter(X %in% LM_fixed_field_data_processed_trees$X) #create a dataframe with only the focal tree
  LM_fixed_field_data_processed_neighbors <- LM_fixed_field_data_processed_trees %>%
    filter(X %notin%  LM_fixed_field_data_focal_tree$X) #create a dataframe with only the neighbors within each buffer
  sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
  sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
  sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
  sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
  sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  
  for (i in 1:nrow(LM_fixed_field_data_processed_neighbors)){ #adding the size values of each neighbor to a sum total of the neighbors size values
    sum_SCA_over_distance = sum_SCA_over_distance + LM_fixed_field_data_processed_neighbors$SCA_over_distance[i] #summing the SCA of each neighbor
    sum_LCA_over_distance = sum_LCA_over_distance + LM_fixed_field_data_processed_neighbors$LCA_over_distance[i] #summing the LCA of each neighbor
    sum_CA_over_distance = sum_CA_over_distance + LM_fixed_field_data_processed_neighbors$CA_over_distance[i] #summing the CA of each neighbor
    sum_CS_over_distance = sum_CS_over_distance + LM_fixed_field_data_processed_neighbors$CS_over_distance[i] #summing the CS of each neighbor
    sum_DBH_over_distance = sum_DBH_over_distance + LM_fixed_field_data_processed_neighbors$DBH_over_distance[i] #summing the DBH of each neighbor
  }
  
  #creating a tibble with all of the calculated sizes over distances 
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  LM_fixed_field_data_focal_tree <- cbind(LM_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  LM_fixed_field_data_all_focal_trees <- rbind(LM_fixed_field_data_all_focal_trees, LM_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
}
View(LM_fixed_field_data_all_focal_trees)
dist
LM_fixed_field_data_all_focal_trees %>%
  filter(is.na(sum_SCA_over_distance))
ggplot()+
  geom_sf(data=LM_fixed_field_data_all_focal_trees, color = LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance)


