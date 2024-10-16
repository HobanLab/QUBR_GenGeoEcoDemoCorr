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
library(rstatix) #to run the Games-Howell Test


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



#### Load in topographic rasters ####


#Importing the cropped rasters for LM, LC, and SD and setting crs
CEM_15_utm_LM <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LM.tif"))
terra::crs(CEM_15_utm_LM) <- CRS("+init=epsg:26912")
names(CEM_15_utm_LM) <- "elev"

CEM_15_utm_LC <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_LC.tif"))
terra::crs(CEM_15_utm_LC) <- CRS("+init=epsg:26912")
names(CEM_15_utm_LC) <- "elev"

CEM_15_utm_SD <- raster(paste0("./data/15 m Elevation Raster/CEM_15_utm_SD.tif"))
terra::crs(CEM_15_utm_SD) <- CRS("+init=epsg:26912")
names(CEM_15_utm_SD) <- "elev"

#extracting the slope in degrees, using the queens method (neighbor = 8)
LM_slope_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'slope')
LC_slope_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'slope')
SD_slope_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'slope')

#extracting the aspect in degrees, using the queens method (neighbor = 8)
LC_aspect_raster_15 <- terra::terrain(CEM_15_utm_LM, unit = 'degrees', neighbors = 8, 'aspect')
LM_aspect_raster_15 <- terra::terrain(CEM_15_utm_LC, unit = 'degrees', neighbors = 8, 'aspect')
SD_aspect_raster_15 <- terra::terrain(CEM_15_utm_SD, unit = 'degrees', neighbors = 8, 'aspect')

#stacking the rasters so I can extract all the info I want at once
LM_topo_stack <- stack(LM_slope_raster_15, LC_aspect_raster_15, CEM_15_utm_LM)

LC_topo_stack <- stack(LC_slope_raster_15, LM_aspect_raster_15, CEM_15_utm_LC)

SD_topo_stack <- stack(SD_slope_raster_15, SD_aspect_raster_15, CEM_15_utm_SD)

#extracting values from topo rasters
LM_data_processed_with_topo <- cbind(LM_fixed_field_data_processed, extract(LM_topo_stack, LM_fixed_field_data_processed)) 

LC_data_processed_with_topo <- cbind(LC_fixed_field_data_processed, extract(LC_topo_stack, LC_fixed_field_data_processed)) 

SD_data_processed_with_topo <- cbind(SD_fixed_field_data_processed, extract(SD_topo_stack, SD_fixed_field_data_processed)) 

#### Selecting single tree per 15m grid cell so that there is no spatial autocorrelation here ####

##LM##

#creating a grid over the soil cells
LM_tree_grid_cropped <- st_make_grid(LM_topo_stack, cellsize = c(15, 15))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(CEM_15_utm_LM, xy = T), aes(x=x, y=y, fill = elev))#+
  geom_sf(data = LM_tree_grid_cropped, fill = NA)


#selecting a point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, fixed_field_data_processed_sf_trans_coordinates, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
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
LM_list_grids_and_point_trees_df <- as.data.frame(unlist(LM_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LM_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
LM_list_grids_and_trees_fixed <- LM_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
LM_data_processed_with_topo_thinned <- LM_data_processed_with_topo %>%
  filter(X %in% LM_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data= LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_data_processed_with_topo_thinned, color = "red")


##LC##

#creating a grid over the soil cells
LC_tree_grid_cropped <- st_make_grid(LC_topo_stack, cellsize = c(15, 15))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(CEM_15_utm_LC, xy = T), aes(x=x, y=y, fill = elev))+
  geom_sf(data = LC_tree_grid_cropped, fill = NA)


#selecting a point from each grid cell with trees within them
LC_list_grids_and_points <- st_contains(LC_tree_grid_cropped, fixed_field_data_processed_sf_trans_coordinates, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
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
LC_list_grids_and_point_trees_df <- as.data.frame(unlist(LC_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LC_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
LC_list_grids_and_trees_fixed <- LC_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
LC_data_processed_with_topo_thinned <- LC_data_processed_with_topo %>%
  filter(X %in% LC_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data= LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_data_processed_with_topo_thinned, color = "red")


##SD##

#creating a grid over the soil cells
SD_tree_grid_cropped <- st_make_grid(SD_topo_stack, cellsize = c(15, 15))

#plotting the grid over an example soil raster
ggplot()+
  geom_raster(data= as.data.frame(CEM_15_utm_SD, xy = T), aes(x=x, y=y, fill = elev))+
geom_sf(data = SD_tree_grid_cropped, fill = NA)


#selecting a point from each grid cell with trees within them
SD_list_grids_and_points <- st_contains(SD_tree_grid_cropped, fixed_field_data_processed_sf_trans_coordinates, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
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
SD_list_grids_and_point_trees_df <- as.data.frame(unlist(SD_list_grids_and_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(SD_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name 
SD_list_grids_and_trees_fixed <- SD_list_grids_and_point_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
SD_data_processed_with_topo_thinned <- SD_data_processed_with_topo %>%
  filter(X %in% SD_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = SD_tree_grid_cropped)+
  geom_sf(data= SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_data_processed_with_topo_thinned, color = "red")


all_data_processed_with_topo_thinned <- rbind(LM_data_processed_with_topo_thinned, LC_data_processed_with_topo_thinned, SD_data_processed_with_topo_thinned)


### Comparing the soil metrics between populations ###
