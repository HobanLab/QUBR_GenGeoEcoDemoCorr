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


fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 or equal area projection UTM 12N
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

#add average nearest neighbor for each individual column
fixed_field_data_processed_NN <- fixed_field_data_processed %>%
  mutate(dist1 = nndist(X = long, Y= lat, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
  mutate(dist2 = nndist(X = long, Y= lat, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
  mutate(dist3 = nndist(X = long, Y= lat, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
  mutate(dist4 = nndist(X = long, Y= lat, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
  mutate(dist5 = nndist(X = long, Y= lat, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
  rowwise()%>% #so that in the next part we take the averages across rows
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5)))%>% #creates a column of the average distances (1-5) of each individual
  select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances

mean(c(1.374369e-05,3.064129e-05,8.661858e-05,2.288619e-04, 2.312526e-04)) #Average nearest neighbors for first row
View(fixed_field_data_processed_NN)



#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "SD")