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
library(MuMIn) #to be able to use model.sel for fitting linear models with spatial autocorrelation

`%notin%` <- Negate(`%in%`) # Make a function that is the opposite of the %in% function

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#transforming the shapefile of trees from WGS84 into equal area projection UTM 12N
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

#creating shapefiles for each population, turning sf of all points into sfc

LM_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LM") %>%
  st_as_sf()

LC_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "LC") %>%
  st_as_sf()

SD_fixed_field_data_processed_sf <- fixed_field_data_processed_sf_transformed %>%
  filter(Locality == "SD") %>%
  st_as_sf()

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

mean(c(1.405577,3.354128,8.840866,25.245919,25.470333))
View(fixed_field_data_processed_NN_UTM)


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

#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)

#checking for duplicates
duplicates <- fixed_field_data_processed_NN_UTM %>% #creates a dataframe called duplicates that filters out X.1 and Y if they have duplicates
  filter(duplicated(X.1) == TRUE) %>%
  filter(duplicated(Y) == TRUE)
View(duplicates)
which(duplicated(fixed_field_data_processed_NN_UTM$X.1)) #finds which rows have duplicates and returns the row number
which(duplicated(fixed_field_data_processed_NN_UTM$Y))  #finds which rows have duplicates and returns the row number

  
#### Moran's I ####

#Test for all points

#regular Moran's I

#creating a matrix of distances between trees where the higher values are at the top
tree.dist <- as.matrix(dist(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))) #making a matrix of the distances between trees
tree.dist.inv <- 1/tree.dist #makes it so closer trees are higher in the matrix
diag(tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
tree.dist.inv[is.infinite(tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another


#computing the Global Moran's I statistic
Moran.I(fixed_field_data_processed_NN_UTM$Canopy_short, tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

tree.coord.matrix <- as.matrix(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))

max(fixed_field_data_processed_NN_UTM$X.1)

#creates nearest neighbor knn using a matrix of the tree coordinates within a specific radius of each tree
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist <- nb2listwdist(knn.dist, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
View(lw.dist)
#checks the neighbor weights for the first tree
lw.dist$weights[1]
#creating lags, which computes the average neighboring short canopy axis for each tree
fixed_field_data_processed_NN_UTM$lag.canopy.short <- lag.listw(lw.dist, fixed_field_data_processed_NN_UTM$Canopy_short)

# Create a regression model
M <- lm(lag.canopy.short ~ Canopy_short, fixed_field_data_processed_NN_UTM)

# Plot the lagged variable vs. the variable 
ggplot(data=fixed_field_data_processed_NN_UTM, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(fixed_field_data_processed_NN_UTM$Canopy_short, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))

#assessing statistical significance with a Monte-Carlo simulation
MC<- moran.mc(fixed_field_data_processed_NN_UTM$Canopy_short, lw.dist, nsim = 999)
MC

#plot of simulated Moran's I values against our value
plot(MC, main="", las=1)
MC$p.value

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local <- localmoran_perm(fixed_field_data_processed_NN_UTM$Canopy_short, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")

#calculating the p-values for each individual tree Moran's I, observed vs. expected
fixed_field_data_processed_NN_UTM$p  <- MC_local.df$`Pr(folded) Sim`

###Test for LM###

#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
LM.tree.dist <- as.matrix(dist(cbind(LM_fixed_field_data_processed$X.1, 
                                  LM_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
LM.tree.dist.inv <- 1/LM.tree.dist #makes it so closer trees are higher in the matrix
diag(LM.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
LM.tree.dist.inv[is.infinite(LM.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_short, LM.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

LM.tree.coord.matrix <- as.matrix(cbind(LM_fixed_field_data_processed$X.1, 
                                     LM_fixed_field_data_processed$Y))

#creates nearest neighbor knn using a matrix of the tree coordinates within a specific radius of each tree
knn.dist.LM <- dnearneigh(LM.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))


#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.LM <- nb2listwdist(knn.dist.LM, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
View(lw.dist.LM)

#checks the neighbor weights for the first tree
lw.dist.LM$weights[4]

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.LM <- lm(lag.canopy.short ~ Canopy_short, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_short, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))
View(sum(lw.dist.LM$weights))
by(lw.dist.LM$weights[4])

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.short <- moran.mc(LM_fixed_field_data_processed$Canopy_short, lw.dist.LM, nsim = 999)
MC.LM.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.LM.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.short <- localmoran_perm(LM_fixed_field_data_processed$Canopy_short, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.short.df <- as.data.frame(MC_local.LM.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.short  <- MC_local.LM.canopy.short.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

LM_fixed_field_data_processed_sign <- LM_fixed_field_data_processed %>%
  filter(pval_sig == T)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for SCA")

#attempting to zoom on the sizes of the significant point
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(size = Canopy_short)) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.short.adjusted))+
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  xlim(c(585700.6, 585903.6))+
  ylim(c(2654803,2654983))
max(LM_fixed_field_data_processed$X.1)
max(LM_fixed_field_data_processed$Y)

#looking at whether similarity is larger or smaller values 


###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_long, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.LM.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_long, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.long <- moran.mc(LM_fixed_field_data_processed$Canopy_long, lw.dist.LM, nsim = 999)
MC.LM.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.LM.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.long <- localmoran_perm(LM_fixed_field_data_processed$Canopy_long, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.long.df <- as.data.frame(MC_local.LM.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.long  <- MC_local.LM.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.long, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Crown_spread, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.LM.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Crown_spread, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.crown.spread <- moran.mc(LM_fixed_field_data_processed$Crown_spread, lw.dist.LM, nsim = 999)
MC.LM.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.LM.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.LM.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.crown.spread <- localmoran_perm(LM_fixed_field_data_processed$Crown_spread, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.crown.spread.df <- as.data.frame(MC_local.LM.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.crown.spread  <- MC_local.LM.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(LM_fixed_field_data_processed$p.crown.spread, 
                                                                 method = "fdr", n=length(LM_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_area, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.LM.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_area, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.area <- moran.mc(LM_fixed_field_data_processed$Canopy_area, lw.dist.LM, nsim = 999)
MC.LM.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.LM.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.area <- localmoran_perm(LM_fixed_field_data_processed$Canopy_area, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.area.df <- as.data.frame(MC_local.LM.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.area <- MC_local.LM.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.area, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$DBH_ag, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.LM, LM_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.LM.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$DBH_ag, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.dbh.ag <- moran.mc(LM_fixed_field_data_processed$DBH_ag, lw.dist.LM, nsim = 999)
MC.LM.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.LM.dbh.ag, main="", las=1, xlab = "DBH")
MC.LM.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.dbh.ag <- localmoran_perm(LM_fixed_field_data_processed$DBH_ag, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.dbh.ag.df <- as.data.frame(MC_local.LM.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.dbh.ag <- MC_local.LM.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(LM_fixed_field_data_processed$p.dbh.ag, 
                                                                 method = "fdr", n=length(LM_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for DBH")

###Test for LC###

#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
LC.tree.dist <- as.matrix(dist(cbind(LC_fixed_field_data_processed$X.1, 
                                     LC_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
LC.tree.dist.inv <- 1/LC.tree.dist #makes it so closer trees are higher in the matrix
diag(LC.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
LC.tree.dist.inv[is.infinite(LC.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_short, LC.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

LC.tree.coord.matrix <- as.matrix(cbind(LC_fixed_field_data_processed$X.1, 
                                        LC_fixed_field_data_processed$Y))


#creates nearest neighbor knn using a matrix of the tree coordinates within a specific radius of each tree
knn.dist.LC <- dnearneigh(LC.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))


#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.LC <- nb2listwdist(knn.dist.LC, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
View(lw.dist.LC)

#checks the neighbor weights for the first tree
lw.dist.LC$weights[1]

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.LC <- lm(lag.canopy.short ~ Canopy_short, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_short, listw = lw.dist.LC, n = length(lw.dist.LC$neighbours), S0 = Szero(lw.dist.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.short <- moran.mc(LC_fixed_field_data_processed$Canopy_short, lw.dist.LC, nsim = 999)
MC.LC.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.LC.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.short <- localmoran_perm(LC_fixed_field_data_processed$Canopy_short, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.short.df <- as.data.frame(MC_local.LC.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.short  <- MC_local.LC.canopy.short.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_long, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.LC.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_long, listw = lw.dist.LC, n = length(lw.dist.LC$neighbours), S0 = Szero(lw.dist.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.long <- moran.mc(LC_fixed_field_data_processed$Canopy_long, lw.dist.LC, nsim = 999)
MC.LC.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.LC.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.long <- localmoran_perm(LC_fixed_field_data_processed$Canopy_long, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.long.df <- as.data.frame(MC_local.LC.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.long  <- MC_local.LC.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.long, 
                                                                 method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Crown_spread, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.LC.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Crown_spread, listw = lw.dist.LC, n = length(lw.dist.LC$neighbours), S0 = Szero(lw.dist.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.crown.spread <- moran.mc(LC_fixed_field_data_processed$Crown_spread, lw.dist.LC, nsim = 999)
MC.LC.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.LC.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.LC.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.crown.spread <- localmoran_perm(LC_fixed_field_data_processed$Crown_spread, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.crown.spread.df <- as.data.frame(MC_local.LC.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.crown.spread  <- MC_local.LC.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(LC_fixed_field_data_processed$p.crown.spread, 
                                                                  method = "fdr", n=length(LC_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_area, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.LC.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_area, listw = lw.dist.LC, n = length(lw.dist.LC$neighbours), S0 = Szero(lw.dist.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.area <- moran.mc(LC_fixed_field_data_processed$Canopy_area, lw.dist.LC, nsim = 999)
MC.LC.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.LC.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.area <- localmoran_perm(LC_fixed_field_data_processed$Canopy_area, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.area.df <- as.data.frame(MC_local.LC.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.area <- MC_local.LC.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.area, 
                                                                 method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$DBH_ag, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.LC.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$DBH_ag, listw = lw.dist.LC, n = length(lw.dist.LC$neighbours), S0 = Szero(lw.dist.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.dbh.ag <- moran.mc(LC_fixed_field_data_processed$DBH_ag, lw.dist.LC, nsim = 999)
MC.LC.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.LC.dbh.ag, main="", las=1, xlab = "DBH")
MC.LC.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.dbh.ag <- localmoran_perm(LC_fixed_field_data_processed$DBH_ag, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.dbh.ag.df <- as.data.frame(MC_local.LC.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.dbh.ag <- MC_local.LC.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(LC_fixed_field_data_processed$p.dbh.ag, 
                                                            method = "fdr", n=length(LC_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

#Test for SD


#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
SD.tree.dist <- as.matrix(dist(cbind(SD_fixed_field_data_processed$X.1, 
                                     SD_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
SD.tree.dist.inv <- 1/SD.tree.dist #makes it so closer trees are higher in the matrix
diag(SD.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
SD.tree.dist.inv[is.infinite(SD.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_short, SD.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

SD.tree.coord.matrix <- as.matrix(cbind(SD_fixed_field_data_processed$X.1, 
                                        SD_fixed_field_data_processed$Y))

#creates nearest neighbor knn using a matrix of the tree coordinates within a specific radius of each tree
knn.dist.SD <- dnearneigh(SD.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))


#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.SD <- nb2listwdist(knn.dist.SD, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
View(lw.dist.SD)

#checks the neighbor weights for the first tree
lw.dist.SD$weights[1]

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.SD <- lm(lag.canopy.short ~ Canopy_short, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_short, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.short <- moran.mc(SD_fixed_field_data_processed$Canopy_short, lw.dist.SD, nsim = 999)
MC.SD.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.SD.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.short <- localmoran_perm(SD_fixed_field_data_processed$Canopy_short, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.short.df <- as.data.frame(MC_local.SD.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.short  <- MC_local.SD.canopy.short.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_long, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.SD.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_long, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.long <- moran.mc(SD_fixed_field_data_processed$Canopy_long, lw.dist.SD, nsim = 999)
MC.SD.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.SD.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.long <- localmoran_perm(SD_fixed_field_data_processed$Canopy_long, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.long.df <- as.data.frame(MC_local.SD.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.long  <- MC_local.SD.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.long, 
                                                                 method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Crown_spread, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.SD.crown.spread <- lm(lag.crown.spread ~ Crown_spread, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Crown_spread, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.crown.spread <- moran.mc(SD_fixed_field_data_processed$Crown_spread, lw.dist.SD, nsim = 999)
MC.SD.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.SD.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.SD.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.crown.spread <- localmoran_perm(SD_fixed_field_data_processed$Crown_spread, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.crown.spread.df <- as.data.frame(MC_local.SD.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.crown.spread  <- MC_local.SD.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(SD_fixed_field_data_processed$p.crown.spread, 
                                                                  method = "fdr", n=length(SD_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_area, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.SD.canopy.area <- lm(lag.canopy.area ~ Canopy_area, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_area, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.area <- moran.mc(SD_fixed_field_data_processed$Canopy_area, lw.dist.SD, nsim = 999)
MC.SD.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.SD.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.area <- localmoran_perm(SD_fixed_field_data_processed$Canopy_area, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.area.df <- as.data.frame(MC_local.SD.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.area <- MC_local.SD.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.area, 
                                                                 method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$DBH_ag, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.SD.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$DBH_ag, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.dbh.ag <- moran.mc(SD_fixed_field_data_processed$DBH_ag, lw.dist.SD, nsim = 999)
MC.SD.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.SD.dbh.ag, main="", las=1, xlab = "DBH")
MC.SD.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.dbh.ag <- localmoran_perm(SD_fixed_field_data_processed$DBH_ag, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.dbh.ag.df <- as.data.frame(MC_local.SD.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.dbh.ag <- MC_local.SD.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(SD_fixed_field_data_processed$p.dbh.ag, 
                                                            method = "fdr", n=length(SD_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CA")


#### Linear Model ####

#creating columns with transformations: logged all of the variables
fixed_field_data_processed_NN_UTM_log <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))%>%
  mutate(ANN = log(ANN))
View(fixed_field_data_processed_NN_UTM_log)

#creating columns with transformations: square root all of the variables
fixed_field_data_processed_NN_UTM_sqrt <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = sqrt(Canopy_short))%>%
  mutate(Canopy_long_lg = sqrt(Canopy_long))%>%
  mutate(Canopy_area_lg = sqrt(Canopy_area))%>%
  mutate(Crown_spread_lg = sqrt(Crown_spread))%>%
  mutate(DBH_ag_lg = sqrt(DBH_ag))%>%
  mutate(ANN = sqrt(ANN))
View(fixed_field_data_processed_NN_UTM_sqrt)

#creating columns with transformations: square root all of the variables
fixed_field_data_processed_NN_UTM_inverse <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = 1/(Canopy_short))%>%
  mutate(Canopy_long_lg = 1/(Canopy_long))%>%
  mutate(Canopy_area_lg = 1/(Canopy_area))%>%
  mutate(Crown_spread_lg = 1/(Crown_spread))%>%
  mutate(DBH_ag_lg = 1/(DBH_ag))%>%
  mutate(ANN = 1/(ANN))
View(fixed_field_data_processed_NN_UTM_inverse)

  
#Linear Model for all points
  
#The conditions for the linear effects model are that:
#there is lineary, the constant variation in the error, independent errors, and normally distributed errors

#LM

#creating a grid over the points with a 10 m edge buffer
LM_box <- st_bbox(LM_fixed_field_data_processed_sf)

#cropping the tree points down by 20 m on all sides

#creating a cropped bbox 
LM_box_sf <- LM_box %>% #turning the bbox into polygon
  st_as_sfc()
LM_box_spatial <- as(LM_box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
LM_box_spatial_cropped <- raster::crop(LM_box_spatial, extent((LM_box[1]+20), (LM_box[3]-20), (LM_box[2]+20), (LM_box[4]-20))) #cropping the xmin, xmax, ymin, and ymax by 20 m inside
LM_box_sf_cropped <-  LM_box_spatial_cropped %>% #turning the spatial polygon into a polygon
  st_as_sfc()

#cropping the points by the cropped box
LM_fixed_field_data_processed_sf_cropped<- st_crop(LM_fixed_field_data_processed_sf, LM_box_sf_cropped)

#Creating a grid over the cropped tree points 
LM_tree_grid_cropped <- st_make_grid(LM_fixed_field_data_processed_sf_cropped, cellsize = (((40*mean(LM_fixed_field_data_processed$DBH_ag))*2)*2))

#plotting the original box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=LM_box_sf)+ #old box
  geom_sf(data=LM_box_sf_cropped)+ #cropped box
  geom_sf(data=LM_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=LM_fixed_field_data_processed_sf_cropped, color = "red") #old points


#selecting a focal point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, LM_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(25) #setting the seed
LM_list_grids_and_focal_trees <- lapply(LM_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
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
LM_list_grids_and_focal_trees_df <- as.data.frame(unlist(LM_list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LM_list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
LM_list_grids_and_focal_trees_fixed <- LM_list_grids_and_focal_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
LM_fixed_field_data_processed_focal <- LM_fixed_field_data_processed_sf %>%
  filter(X %in% LM_list_grids_and_focal_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

#creating the buffer around the focal points
LM_focal_tree_buffers <-st_buffer(LM_fixed_field_data_processed_focal$geometry, 40*mean(LM_fixed_field_data_processed_focal$DBH_ag))

#graphing the selected focal trees, the buffers, the grid
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data=LM_focal_tree_buffers, color = "blue") +
  geom_sf(data= LM_fixed_field_data_processed_focal, aes(color = X))

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

#Checking that row number in focal dataset is the same as the buffer dataset
LM_fixed_field_data_processed_focal_row <- LM_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
LM_tree_buffer_inside_0 <- mutate(LM_tree_buffer_inside_0, row = as.factor(row)) #making sure the buffers have the same row number as the focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LM_tree_grid_cropped) +
  geom_sf(data=LM_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=LM_fixed_field_data_processed_focal_row, aes(color = row))

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
ggplot()+
  geom_sf(data = LM_focal_tree_buffers)+
  geom_sf(data = LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_fixed_field_data_processed_focal, color = 'blue')

#Checking that row number in focal dataset is the same as the buffer dataset
LM_fixed_field_data_processed_focal_row <- LM_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
LM_tree_buffer_inside <- mutate(LM_tree_buffer_inside, row = as.factor(row)) #making sure the buffers have the same row number as the isoalted focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LM_tree_grid_cropped) +
  geom_sf(data=LM_tree_buffer_inside, aes(color = row))+
  geom_sf(data=LM_fixed_field_data_processed_focal_row, aes(color = row))


LM_fixed_field_data_all_focal_trees <- tibble()#creating the empty tibble 

#calculating the distances of each tree within the buffer to the focal tree and the competition metric values
for (i in 1:nrow(LM_fixed_field_data_processed_focal)){ #for the length of the buffers with trees inside of them #LM_fixed_field_data_processed_focal_row
  row_num = i
  LM_tree_buffer_inside_df <- as.data.frame(LM_tree_buffer_inside_0) #uses data of non-isolated and isolated focal trees 
  LM_tree_buffer_inside_df_i <- LM_tree_buffer_inside_df %>% 
    filter(row == row_num) #isolate a row of the buffer dataframe
  LM_tree_buffer_inside_sf_i <- st_as_sf(LM_tree_buffer_inside_df_i) #set the row as a simple feature
  all_pts_buffer <- st_contains(LM_tree_buffer_inside_sf_i, LM_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
  LM_fixed_field_data_processed_trees <- LM_fixed_field_data_processed %>%
    filter(X %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer
  correct_focal <- LM_list_grids_and_focal_trees_fixed[i,]
  LM_fixed_field_data_focal_tree <- LM_fixed_field_data_processed_focal %>%
    filter(X %in% correct_focal$tree_row_num) #create a dataframe with only the focal tree
  
  LM_fixed_field_data_neighbor_trees <- LM_fixed_field_data_processed_trees %>%
    filter(X %notin% LM_fixed_field_data_focal_tree$X) 
  
  if(nrow(LM_fixed_field_data_neighbor_trees) == 0){
    sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
    sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
    sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
    sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
    sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  } else{
  
  LM_fixed_field_data_neighbor_trees <-  LM_fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
    mutate(focal_distance = as.numeric(st_distance(geometry,  LM_fixed_field_data_focal_tree$geometry, by_element = T))) %>% #calculate the distance between the focal tree and each tree that neighbors it
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

  
  for (y in 1:nrow(LM_fixed_field_data_neighbor_trees)){ #adding the size values of each neighbor to a sum total of the neighbors size values
    sum_SCA_over_distance = sum_SCA_over_distance + LM_fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
    sum_LCA_over_distance = sum_LCA_over_distance + LM_fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the LCA of each neighbor
    sum_CA_over_distance = sum_CA_over_distance + LM_fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
    sum_CS_over_distance = sum_CS_over_distance + LM_fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
    sum_DBH_over_distance = sum_DBH_over_distance + LM_fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
  }
  }
  #creating a tibble with all of the calculated sizes over distances 
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  LM_fixed_field_data_focal_tree <- cbind(LM_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  LM_fixed_field_data_all_focal_trees <- rbind(LM_fixed_field_data_all_focal_trees, LM_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
    
}

View(LM_fixed_field_data_all_focal_trees)

#plotting the tree points and their 
ggplot()+
  geom_sf(data=LM_fixed_field_data_all_focal_trees, aes(color = sum_SCA_over_distance))

#descriptive statistics

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

#conditions are lINES: linearity, independence, normal distribution of residuals, equal variance, simple random sample

#calculate leverage for each observation in the model
leverage <- as.data.frame(hatvalues(LM_lmem_focal_SCA))
levarage <- leverage %>%
  mutate(row = row_number())
plot(hatvalues(LM_lmem_focal_SCA), type = 'h')
mean(hatvalues(LM_lmem_focal_SCA))
unusual_lev <- 4/nrow(LM_fixed_field_data_processed_focal_dist)
which(leverage > unusual)
which(leverage$`hatvalues(LM_lmem_focal_SCA)` > .025)

#Cook's D
cookd(LM_lmem_focal_SCA)
cooksD <- cooks.distance(LM_lmem_focal_SCA)
plot(cooksD, type = 'h')
unsual_cooksD <- 0.5
which(cooksD > unsual_cooksD)

#checking linearity 

#plotting the linear model in ggplot for SCA
ggplot(data = LM_fixed_field_data_all_focal_trees, (aes(x=sum_SCA_over_distance, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("SCA over Distance")+
  ylab("Short Canopy Axis")


#creating the generalized linear effects model

#creating x and y columns of the UTM 12N 
LM_fixed_field_data_all_focal_trees$X.1 <- st_coordinates(LM_fixed_field_data_all_focal_trees)[,1]
LM_fixed_field_data_all_focal_trees$Y <- st_coordinates(LM_fixed_field_data_all_focal_trees)[,2]


#logged version of generalized linear model
LM_gls_focal_SCA_lg <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lg_exp <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), 
                               correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lg_gaus <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), 
                                correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lg_spher <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), 
                                 correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lg_lin <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), 
                               correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lg_ratio <- gls(log(Canopy_short + 1) ~ log(sum_SCA_over_distance + 1), 
                                 correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_lg_SCA <- model.sel(LM_gls_focal_SCA_lg, LM_gls_focal_SCA_lg_exp, LM_gls_focal_SCA_lg_gaus, 
                         LM_gls_focal_SCA_lg_spher, LM_gls_focal_SCA_lg_ratio)
View(LM_AIC_test_lg_SCA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_SCA_lg$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_SCA_lg$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_SCA_lg$fitted, y = LM_gls_focal_SCA_lg$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for logged SCA and logged SCA over Distance")

#checking we have appropriately removed the spatial autocorrelation
semivario <- Variogram( LM_gls_focal_SCA_lg, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm
summary(LM_gls_focal_SCA_lg)

#unlogged version of generlized linear model
LM_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test <- model.sel(LM_gls_focal_SCA, LM_gls_focal_SCA_exp, LM_gls_focal_SCA_gaus, LM_gls_focal_SCA_spher, LM_gls_focal_SCA_lin, LM_gls_focal_SCA_ratio)
View(LM_gls_focal_SCA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_SCA$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees , aes(x = LM_gls_focal_SCA$fitted, y = LM_gls_focal_SCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and SCA over Distance")

#checking we have appropriately removed the spatial autocorrelation
semivario <- Variogram( LM_gls_focal_SCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm
summary(LM_gls_focal_SCA)

