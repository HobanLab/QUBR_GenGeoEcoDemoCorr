# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei compete or facilitate with one another%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by the distance to other individuals of the same species 
# either due to competition or facilitation. 
      # If they are impacted by facilitation, we would expect closer trees would be bigger. 
      # If they are impacted by competition, we would expect closer trees to be smaller. 
# To test this, we used Global and Local Moran's I to determine whether values of SCA, LCA, CS, CA, and DBH 
# that were closer together were more similar in value or not. 
    # The global Moran's I looked for general spatial autocorrelation
    # The local Moran's I looked for areas were values were more similar than other areas. 
# We used also performed a linear regression to see if for focal trees, there was a 
# relationship between how much competition the trees face (based on the 
# size of the neighbors over their distance to the focal trees) and the size of the focal trees.

# The script is broken into sections of 
  # 1) loading and processing the packages and spatial/size/shape data for the trees in the Las Matancitas,
        #San Dionisio, and La Cobriza populations and loading in the river outline shapefiles, 
  # 2) creating a dataframe with the coordinates and average distance between each tree and their 5 nearest neighbors
  # 3) graphing descriptive summary histograms and calculating summary statistics (mean, median, sd, etc.) for each response 
            #variable (SCA, LCA, CA, CS, DBH),
  # 4) Running the global and local Moran's I analyses, 
  # 5) using linear regression to see if tree size seem related to local competition  


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
library(Kendall)# to use the Mann-Kendall test to look for non-parametric correlations in the data

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

#Loading in ArcGIS river shapefile and storing out polygons for each population

#Las Matancitas (LM)
river_LM <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LM River/LM_Rivers_Final.shp")
river_LM  <- river_LM$geometry[1]
plot(river_LM)

#La Cobriza (LC)
river_LC  <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/LC River/LC_Rivers_Final.shp")
river_LC  <- river_LC$geometry[1]
plot(river_LC)

#San Dionisio (SD)
river_SD <- st_read("./data/Shapefiles/FINAL River Shapefiles ArcGIS/SD River/SD_Rivers_Final.shp")
river_SD <- river_SD$geometry[1]
plot(river_SD)

#changing the coordinate reference system of the river polygons to be equal area projection (UTM 12N), uses meters as distance measurement 

river_LM_trans <- st_transform(river_LM, crs = 26912) 
river_LC_trans <- st_transform(river_LC, crs = 26912)
river_SD_trans <- st_transform(river_SD, crs = 26912)


#### Computing Average Nearest Neighbors for each tree ####

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
  #dplyr::select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances


# Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns

LM_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "SD")



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

#checking for duplicates and filtering them out
duplicates <- fixed_field_data_processed_NN_UTM %>% #creates a dataframe called duplicates that filters out X.1 and Y if they have duplicates
  filter(duplicated(X.1) == TRUE) %>%
  filter(duplicated(Y) == TRUE)

  
#### Moran's I ####

# For each response variable (SCA, LCA, CA, CS, DBH), we ran these analyses:
    # 1) Global Moran's I
          # a) We first created a matrix of the inverse distances between each tree 
          # b) computed the average response variable of the neighboring trees for each tree (lagged response variable)
          # c) plot the regression line of the lagged response variable vs. the actual response values for each tree
                  # positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
                  # negative slope, negative spatial autocorrelation, variation in size of trees close together
          # d) Calculate the Moran's I statistic and then use a Monte Carlo test and plot to see if it is significant spatial autocorrelation 
    # 2) Local Moran's I
          # a) calculate and plot the expected local Moran's I vs. local Moran's I for each tree 
          # b) calculate the number of trees with significant local Moran's I (p-values < 0.05)

#Test for all points

#regular Moran's I

#creating a matrix of distances between trees where the higher values are at the top
tree.dist <- as.matrix(dist(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))) #making a matrix of the distances between trees
tree.dist.inv <- 1/tree.dist #makes it so closer trees are higher in the matrix
diag(tree.dist.inv) <- 0 #makes it so trees have a 0 distance with themselves
tree.dist.inv[is.infinite(tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another


#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(fixed_field_data_processed_NN_UTM$Canopy_short, tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))


#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist <- nb2listwdist(knn.dist, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
fixed_field_data_processed_NN_UTM$lag.canopy.short <- lag.listw(lw.dist, fixed_field_data_processed_NN_UTM$Canopy_short)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M <- lm(lag.canopy.short ~ Canopy_short, fixed_field_data_processed_NN_UTM)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=fixed_field_data_processed_NN_UTM, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the global Moran's I statistic
moran(fixed_field_data_processed_NN_UTM$Canopy_short, listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))

#assessing statistical significance of Moran's I statistic with a Monte-Carlo simulation
MC<- moran.mc(fixed_field_data_processed_NN_UTM$Canopy_short, lw.dist, nsim = 999)
MC

#plot of simulated Moran's I values against our value
#if the known value is to the right of the peak, similar sizes/shapes are more clustered together
#if the known value is to the left of the peak, similar sizes/shapes are more dispersed
plot(MC, main="", las=1,  xlab = "Short Canopy Axis")
MC$p.value

#Local Moran's I 

#using the weighted neighbors to simulate size values at random
MC_local <- localmoran_perm(fixed_field_data_processed_NN_UTM$Canopy_short, lw.dist, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local Moran's I values vs. the expected at random
ggplot(data=MC_local.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")

#calculating the p-values for each individual tree Moran's I, observed vs. expected
fixed_field_data_processed_NN_UTM$p  <- MC_local.df$`Pr(folded) Sim`

length(which(fixed_field_data_processed_NN_UTM$p < 0.05))

###Test for LM###

#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
LM.tree.dist <- as.matrix(dist(cbind(LM_fixed_field_data_processed$X.1, 
                                  LM_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
LM.tree.dist.inv <- 1/LM.tree.dist #makes it so closer trees are higher in the matrix
diag(LM.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
LM.tree.dist.inv[is.infinite(LM.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LM_fixed_field_data_processed$Canopy_short, LM.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#creating a matrix of the tree locations
LM.tree.coord.matrix <- as.matrix(cbind(LM_fixed_field_data_processed$X.1, 
                                     LM_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.LM <- dnearneigh(LM.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.LM <- nb2listwdist(knn.dist.LM, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_short)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LM <- lm(lag.canopy.short ~ Canopy_short, LM_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_short, listw = lw.dist.LM, n = length(lw.dist.LM$neighbours), S0 = Szero(lw.dist.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.short <- moran.mc(LM_fixed_field_data_processed$Canopy_short, lw.dist.LM, nsim = 999)
MC.LM.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.LM.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values at random
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

#Number of trees with significant adjusted Moran's I
length(which(LM_fixed_field_data_processed$p.canopy.short.adjusted < 0.05))

#filtering out significant p-values
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

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LM_fixed_field_data_processed$Canopy_long, LM.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_long)
# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LM.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LM_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
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

#Number of trees with significant adjusted Moran's I
length(which(LM_fixed_field_data_processed$p.canopy.short.adjusted < 0.05))


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

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LM_fixed_field_data_processed$Crown_spread, LM.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring crown spread for each tree
LM_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Crown_spread)
# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LM.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LM_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LM.crown.spread <- localmoran_perm(LM_fixed_field_data_processed$Crown_spread, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.crown.spread.df <- as.data.frame(MC_local.LM.crown.spread)

##Ii is local Moran's I statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
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

#Number of trees with significant adjusted Moran's I
length(which(LM_fixed_field_data_processed$p.crown.spread.adjusted < 0.05))

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

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LM_fixed_field_data_processed$Canopy_area, LM.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$Canopy_area)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LM.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LM_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LM.canopy.area <- localmoran_perm(LM_fixed_field_data_processed$Canopy_area, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.area.df <- as.data.frame(MC_local.LM.canopy.area)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local moran statistic  
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

#Number of trees with significant adjusted Moran's I
length(which(LM_fixed_field_data_processed$p.canopy.area.adjusted < 0.05))


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

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LM_fixed_field_data_processed$DBH_ag, LM.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.dist.LM, LM_fixed_field_data_processed$DBH_ag)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LM.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LM_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LM.dbh.ag <- localmoran_perm(LM_fixed_field_data_processed$DBH_ag, lw.dist.LM, nsim = 9999, alternative = "greater")
MC_local.LM.dbh.ag.df <- as.data.frame(MC_local.LM.dbh.ag)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
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
#Number of trees with significant adjusted Moran's I
length(which(LM_fixed_field_data_processed$p.dbh.ag.adjusted < 0.05))

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

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LC_fixed_field_data_processed$Canopy_short, LC.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#creating a matrix of the tree locations
LC.tree.coord.matrix <- as.matrix(cbind(LC_fixed_field_data_processed$X.1, 
                                        LC_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.LC <- dnearneigh(LC.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))


#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.LC <- nb2listwdist(knn.dist.LC, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#checks the neighbor weights for the first tree
lw.dist.LC$weights[1]

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_short)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LC <- lm(lag.canopy.short ~ Canopy_short, LC_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LC.canopy.short <- localmoran_perm(LC_fixed_field_data_processed$Canopy_short, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.short.df <- as.data.frame(MC_local.LC.canopy.short)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local moran statistic  
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

length(which(LC_fixed_field_data_processed$p.canopy.short.adjusted < 0.05))

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LC_fixed_field_data_processed$Canopy_long, LC.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring long canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_long)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LC.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LC_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LC.canopy.long <- localmoran_perm(LC_fixed_field_data_processed$Canopy_long, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.long.df <- as.data.frame(MC_local.LC.canopy.long)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
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

#number of trees with significant local clustering
length(which(LC_fixed_field_data_processed$p.canopy.long.adjusted < 0.05))


###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LC_fixed_field_data_processed$Crown_spread, LC.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring crown spread for each tree
LC_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Crown_spread)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LC.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LC_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LC.crown.spread <- localmoran_perm(LC_fixed_field_data_processed$Crown_spread, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.crown.spread.df <- as.data.frame(MC_local.LC.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local Moran's I statistic  
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

#number of trees with significant of clustering
length(which(LC_fixed_field_data_processed$p.crown.spread.adjusted < 0.05))


###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LC_fixed_field_data_processed$Canopy_area, LC.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring canopy area for each tree
LC_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$Canopy_area)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LC.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LC_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LC.canopy.area <- localmoran_perm(LC_fixed_field_data_processed$Canopy_area, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.area.df <- as.data.frame(MC_local.LC.canopy.area)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
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

length(which(LC_fixed_field_data_processed$p.canopy.area.adjusted < 0.05))


###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(LC_fixed_field_data_processed$DBH_ag, LC.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.dist.LC, LC_fixed_field_data_processed$DBH_ag)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.LC.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LC_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.LC.dbh.ag <- localmoran_perm(LC_fixed_field_data_processed$DBH_ag, lw.dist.LC, nsim = 9999, alternative = "greater")
MC_local.LC.dbh.ag.df <- as.data.frame(MC_local.LC.dbh.ag)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
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

#number of trees with significant clustering
length(which(LC_fixed_field_data_processed$p.dbh.ag.adjusted < 0.05))



#Test for SD


#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
SD.tree.dist <- as.matrix(dist(cbind(SD_fixed_field_data_processed$X.1, 
                                     SD_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
SD.tree.dist.inv <- 1/SD.tree.dist #makes it so closer trees are higher in the matrix
diag(SD.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
SD.tree.dist.inv[is.infinite(SD.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(SD_fixed_field_data_processed$Canopy_short, SD.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#creating a matrix of the tree locations
SD.tree.coord.matrix <- as.matrix(cbind(SD_fixed_field_data_processed$X.1, 
                                        SD_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.SD <- dnearneigh(SD.tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))


#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.dist.SD <- nb2listwdist(knn.dist.SD, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                           alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_short)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.SD <- lm(lag.canopy.short ~ Canopy_short, SD_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
MC_local.SD.canopy.short <- localmoran_perm(SD_fixed_field_data_processed$Canopy_short, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.short.df <- as.data.frame(MC_local.SD.canopy.short)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
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

#checking how many points are significant
length(which(SD_fixed_field_data_processed$p.canopy.short.adjusted <0.05))

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(SD_fixed_field_data_processed$Canopy_long, SD.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring long canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_long)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.SD.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, SD_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

#using the weighted neighbors to simulate size values at random
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


#checking how many points are significant
length(which(SD_fixed_field_data_processed$p.canopy.long.adjusted <0.05))


###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(SD_fixed_field_data_processed$Crown_spread, SD.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring crown spread for each tree
SD_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Crown_spread)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.SD.crown.spread <- lm(lag.crown.spread ~ Crown_spread, SD_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")+
  theme(text = element_text(size = 20))

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Crown_spread, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.crown.spread <- moran.mc(SD_fixed_field_data_processed$Crown_spread, lw.dist.SD, nsim = 999)
MC.SD.crown.spread
 
#plot of simulated Moran's I values against our value
plot(MC.SD.crown.spread, main="", las=1, xlab = "Crown Spread", cex.lab = 1.5, col = "red" )
lines(MC.SD.crown.spread$statistic, col = "red")
MC.SD.crown.spread$p.value #extracting the pvalue
MC.SD.crown.spread$statistic
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

#checking how many points are significant
length(which(SD_fixed_field_data_processed$p.crown.spread.adjusted <0.05))


###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(SD_fixed_field_data_processed$Canopy_area, SD.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring canopy area for each tree
SD_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$Canopy_area)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.SD.canopy.area <- lm(lag.canopy.area ~ Canopy_area, SD_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
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

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
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


#checking how many points are significant
length(which(SD_fixed_field_data_processed$p.canopy.area.adjusted <0.05))



###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic with one function Moran.I()
Moran.I(SD_fixed_field_data_processed$DBH_ag, SD.tree.dist.inv)

#creating lags for each tree, which computes the average neighboring DBH for each tree
SD_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.dist.SD, SD_fixed_field_data_processed$DBH_ag)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M.SD.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, SD_fixed_field_data_processed)

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#creating lags for each tree, which computes the average neighboring short canopy axis for each tree
fixed_field_data_processed_NN_UTM$lag.canopy.short <- lag.listw(lw.dist, fixed_field_data_processed_NN_UTM$Canopy_short)

# Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
M <- lm(lag.canopy.short ~ Canopy_short, fixed_field_data_processed_NN_UTM)

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$DBH_ag, listw = lw.dist.SD, n = length(lw.dist.SD$neighbours), S0 = Szero(lw.dist.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.dbh.ag <- moran.mc(SD_fixed_field_data_processed$DBH_ag, lw.dist.SD, nsim = 999)
MC.SD.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.SD.dbh.ag, main="", las=1, xlab = "DBH")
MC.SD.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values at random
MC_local.SD.dbh.ag <- localmoran_perm(SD_fixed_field_data_processed$DBH_ag, lw.dist.SD, nsim = 9999, alternative = "greater")
MC_local.SD.dbh.ag.df <- as.data.frame(MC_local.SD.dbh.ag)

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
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


#checking how many points are significant
length(which(SD_fixed_field_data_processed$p.dbh.ag.adjusted <0.05))


#### Linear Model ####

# To see if trees that face more competition (they face closer and larger trees) are smaller, for each populatio,
      # 1) we create a bounding box around the points and then cropped the bounding box of the points by 20 m on all sides to avoid edge effects
      # 2) randomly select a tree (focal tree) from each grid cell with a tree in it
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
          # f) use a non-parametric Mann-Kendall Test to see if there is a significant correlation between the competition metric
                #and the size of the focal trees for the data with no outliers


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


#LM

#creating a grid over the points with a 10 m edge buffer

# creating a bounding box based on the point locations
LM_box <- st_bbox(LM_fixed_field_data_processed_sf)

#cropping the tree points down by 20 m on all sides

#creating a cropped bbox 
LM_box_sf <- LM_box %>% #turning the bbox into polygon
  st_as_sfc()
LM_box_spatial <- as(LM_box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
LM_box_spatial_cropped <- raster::crop(LM_box_spatial, extent((LM_box[1]+20), (LM_box[3]-20), (LM_box[2]+20), (LM_box[4]-20))) #cropping the bounding box xmin, xmax, ymin, and ymax by 20 m inside
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


#randomly selecting a focal point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, LM_fixed_field_data_processed_sf, sparse =T) #find which points are within which grid cells, make sure row number in the data frame of grid cells corresponds to the order of the points dataframe within st_contains
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

# in this loop, it iterates over each focal tree with neighbors and calculates the sum of the size metric 
#divided by the distance for all of the neighbors for each focal tree (competition values for each tree)

for (i in 1:nrow(LM_fixed_field_data_processed_focal)){ #for the length focal trees 
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
  
  #isolating the neighbor tree data
  LM_fixed_field_data_neighbor_trees <- LM_fixed_field_data_processed_trees %>%
    filter(X %notin% LM_fixed_field_data_focal_tree$X) 
  
  #if there are no neighbors, its sets the sum of the response variable divided by the distance of the tree to the focal tree to 0
  if(nrow(LM_fixed_field_data_neighbor_trees) == 0){
    sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
    sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
    sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
    sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
    sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  } else{
    
    # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
    LM_fixed_field_data_neighbor_trees <-  LM_fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
      mutate(focal_distance = as.numeric(st_distance(geometry,  LM_fixed_field_data_focal_tree$geometry, by_element = T))) %>% #calculate the distance between the focal tree and each tree that neighbors it
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
    for (y in 1:nrow(LM_fixed_field_data_neighbor_trees)){ 
      sum_SCA_over_distance = sum_SCA_over_distance + LM_fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
      sum_LCA_over_distance = sum_LCA_over_distance + LM_fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the LCA of each neighbor
      sum_CA_over_distance = sum_CA_over_distance + LM_fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
      sum_CS_over_distance = sum_CS_over_distance + LM_fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
      sum_DBH_over_distance = sum_DBH_over_distance + LM_fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
    }
  }
  #creating a tibble with all of the calculated sizes over distances and other tree attributes for each focal tree
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  LM_fixed_field_data_focal_tree <- cbind(LM_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  LM_fixed_field_data_all_focal_trees <- rbind(LM_fixed_field_data_all_focal_trees, LM_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
  
}
View(LM_fixed_field_data_all_focal_trees)


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

#creating the generalized linear effects model

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
LM_fixed_field_data_all_focal_trees_no_SCA_outliers <- LM_fixed_field_data_all_focal_trees[-c(10, 16),]

#removing the outlier sum of size over distance values skewing the results
LM_fixed_field_data_all_focal_trees <- LM_fixed_field_data_all_focal_trees %>%
  filter(sum_SCA_over_distance < 2)

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LM_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test <- model.sel(LM_gls_focal_SCA, LM_gls_focal_SCA_exp, LM_gls_focal_SCA_gaus, LM_gls_focal_SCA_spher, LM_gls_focal_SCA_lin, LM_gls_focal_SCA_ratio)
LM_AIC_test

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
shapiro.test(LM_gls_focal_SCA$residuals) #shapiro wilk test, not significant so it is normal 

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

#non parametric Mann-Kendall Test for the version without outliers
LM_tau_result_SCA <- cor.test(LM_fixed_field_data_all_focal_trees$sum_SCA_over_distance, LM_fixed_field_data_all_focal_trees$Canopy_short,  method = "kendall")

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


#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)LM_gls_focal_LCA <- gls(Canopy_short ~ sum_LCA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_exp <- gls(Canopy_short ~ sum_LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_gaus <- gls(Canopy_short ~ sum_LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_spher <- gls(Canopy_short ~ sum_LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_lin <- gls(Canopy_short ~ sum_LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_LCA_ratio <- gls(Canopy_short ~ sum_LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_LCA <- model.sel(LM_gls_focal_LCA, LM_gls_focal_LCA_exp, LM_gls_focal_LCA_gaus, LM_gls_focal_LCA_spher, LM_gls_focal_LCA_lin, LM_gls_focal_LCA_ratio)
View(LM_AIC_test_LCA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_LCA$residuals))+
  geom_qq()

shapiro.test(LM_gls_focal_LCA$residuals) #shaprio wilk test, not sign so our residuals are normally distribtued

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

#non parametric Mann-Kendall Test for the version without outliers
LM_tau_result_LCA <- cor.test(LM_fixed_field_data_all_focal_trees$sum_LCA_over_distance, LM_fixed_field_data_all_focal_trees$Canopy_long,  method = "kendall")

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


#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)LM_gls_focal_CA <- gls(Canopy_short ~ sum_CA_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_exp <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_gaus <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_spher <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_lin <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CA_ratio <- gls(Canopy_short ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CA <- model.sel(LM_gls_focal_CA, LM_gls_focal_CA_exp, LM_gls_focal_CA_gaus, LM_gls_focal_CA_spher, LM_gls_focal_CA_ratio,LM_gls_focal_CA_lin) #LM_gls_focal_CA_lin
View(LM_AIC_test_CA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees, aes(x= LM_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LM_fixed_field_data_all_focal_trees, aes(sample = LM_gls_focal_CA$residuals))+
  geom_qq()

shapiro.test(LM_gls_focal_CA$residuals) # shapiro-wilk, not sign so the residuals are normally distributed

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
summary(LM_AIC_test_CA)


#non parametric Mann-Kendall Test for the version without outliers
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

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)LM_gls_focal_CS <- gls(Canopy_short ~ sum_CS_over_distance, data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_exp <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_gaus <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_spher <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_lin <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)
LM_gls_focal_CS_ratio <- gls(Canopy_short ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_CS <- model.sel(LM_gls_focal_CS, LM_gls_focal_CS_lin, LM_gls_focal_CS_exp, LM_gls_focal_CS_gaus, LM_gls_focal_CS_spher, LM_gls_focal_CS_ratio) #without linear correlation
View(LM_AIC_test_CS)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees_no_CS_outliers, aes(x= LM_gls_focal_CS$residuals))+
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
summary(LM_AIC_test_CA)

#non parametric Mann-Kendall Test for the version without outliers
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
LM_gls_focal_DBH <- gls(Canopy_short ~ sum_DBH_over_distance, data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)
LM_gls_focal_DBH_exp <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)
LM_gls_focal_DBH_gaus <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)
LM_gls_focal_DBH_spher <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)
LM_gls_focal_DBH_lin <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)
LM_gls_focal_DBH_ratio <- gls(Canopy_short ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LM_fixed_field_data_all_focal_trees_no_CS_outliers)

#ordering models by which ones have the lowest Akaike information criterion
LM_AIC_test_DHB <- model.sel(LM_gls_focal_DBH, LM_gls_focal_DBH_exp, LM_gls_focal_DBH_lin, LM_gls_focal_DBH_gaus, LM_gls_focal_DBH_spher, LM_gls_focal_DBH_ratio)
View(LM_AIC_test_CS)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LM_fixed_field_data_all_focal_trees_no_CS_outliers, aes(x= LM_gls_focal_DBH$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LM_fixed_field_data_all_focal_trees_no_CS_outliers, aes(sample = LM_gls_focal_DBH$residuals))+
  geom_qq()

# shapiro-wilk, not sign so the residuals are normally distributed
shapiro.test(LM_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = LM_fixed_field_data_all_focal_trees_no_CS_outliers , aes(x = LM_gls_focal_DBH$fitted, y = LM_gls_focal_DBH$residuals))+
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
summary(LM_gls_focal_DBH)

#non parametric Mann-Kendall Test for the version without outliers
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

#creating x and y columns of the UTM 12N 
LC_fixed_field_data_all_focal_trees$X.1 <- st_coordinates(LC_fixed_field_data_all_focal_trees)[,1]
LC_fixed_field_data_all_focal_trees$Y <- st_coordinates(LC_fixed_field_data_all_focal_trees)[,2]


#creating a grid over the points with a 10 m edge buffer

# creating a bounding box based on the point locations
LC_box <- st_bbox(LC_fixed_field_data_processed_sf)

#cropping the tree points down by 20 m on all sides

#creating a cropped bbox 
LC_box_sf <- LC_box %>% #turning the bbox into polygon
  st_as_sfc()
LC_box_spatial <- as(LC_box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
LC_box_spatial_cropped <- raster::crop(LC_box_spatial, extent((LC_box[1]+20), (LC_box[3]-20), (LC_box[2]+20), (LC_box[4]-20))) #cropping the xmin, xmax, ymin, and ymax by 20 m inside
LC_box_sf_cropped <-  LC_box_spatial_cropped %>% #turning the spatial polygon into a polygon
  st_as_sfc()

#cropping the points by the cropped box
LC_fixed_field_data_processed_sf_cropped<- st_crop(LC_fixed_field_data_processed_sf, LC_box_sf_cropped)

#Creating a grid over the cropped tree points 
LC_tree_grid_cropped <- st_make_grid(LC_fixed_field_data_processed_sf_cropped, cellsize = (((40*mean(LC_fixed_field_data_processed$DBH_ag))*2)*2))

#plotting the original box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=LC_box_sf)+ #old box
  geom_sf(data=LC_box_sf_cropped)+ #cropped box 
  geom_sf(data=LC_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=LC_fixed_field_data_processed_sf_cropped, color = "red") #old points

#creating an x_sequential column that is 1 through the number of LC points
LC_fixed_field_data_processed_sf <- LC_fixed_field_data_processed_sf %>%
  mutate(X_sequential = 1:nrow(LC_fixed_field_data_processed_sf))
View(LC_fixed_field_data_processed_sf)

#selecting a focal point from each grid cell with trees within them
LC_list_grids_and_points <- st_contains(LC_tree_grid_cropped, LC_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(25) #setting the seed
LC_list_grids_and_focal_trees <- lapply(LC_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
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

View(LC_fixed_field_data_processed_sf)

#creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
LC_list_grids_and_focal_trees_df <- as.data.frame(unlist(LC_list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(LC_list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
LC_list_grids_and_focal_trees_fixed <- LC_list_grids_and_focal_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
  mutate(data_row = LC_fixed_field_data_processed_sf$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
LC_fixed_field_data_processed_focal <- LC_fixed_field_data_processed_sf %>%
  filter(X_sequential %in% LC_list_grids_and_focal_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

View(LC_fixed_field_data_processed_focal)

#creating the buffer around the focal points
LC_focal_tree_buffers <-st_buffer(LC_fixed_field_data_processed_focal$geometry, 40*mean(LC_fixed_field_data_processed_focal$DBH_ag))

#graphing the selected focal trees, the buffers, the grid
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data=LC_focal_tree_buffers, color = "blue") +
  geom_sf(data = LC_fixed_field_data_processed) +
  geom_sf(data= LC_fixed_field_data_processed_focal, aes(color = X))

#calculating the size/distance for focal trees and neighbors within buffers for buffers with only the focal tree and with more 

#create a tibble with the the number of trees within the buffers that contain trees
LC_tree_buffers_points_within_0 <- st_contains(LC_focal_tree_buffers, LC_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 0) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
LC_tree_buffer_inside_0 <- LC_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% LC_tree_buffers_points_within_0$row) #only keep polygons that match the row number of the grid cells with trees within them 

#Checking that row number in focal dataset is the same as the buffer dataset
LC_fixed_field_data_processed_focal_row <- LC_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
LC_tree_buffer_inside_0 <- mutate(LC_tree_buffer_inside_0, row = as.factor(row)) #making sure the buffers have the same row number as the focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LC_tree_grid_cropped) +
  geom_sf(data=LC_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=LC_fixed_field_data_processed_focal_row, aes(color = row))

#calculating the size/distance for focal trees and neighbors within buffers for buffers with more than just the focal tree

#create a tibble with the the number of trees within the buffers that contain trees
LC_tree_buffers_points_within <- st_contains(LC_focal_tree_buffers, LC_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 1) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
LC_tree_buffer_inside <- LC_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% LC_tree_buffers_points_within$row) #only keep buffers that match the row number of the buffers cells with trees within them 
View(LC_tree_buffer_inside)

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = LC_focal_tree_buffers)+
  geom_sf(data = LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_fixed_field_data_processed_focal, color = 'blue')

#Checking that row number in focal dataset is the same as the buffer dataset
LC_fixed_field_data_processed_focal_row <- LC_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
LC_tree_buffer_inside <- mutate(LC_tree_buffer_inside, row = as.factor(row)) #making sure the buffers have the same row number as the isoalted focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = LC_tree_grid_cropped) +
  geom_sf(data=LC_tree_buffer_inside, aes(color = row))+
  geom_sf(data=LC_fixed_field_data_processed_focal_row, aes(color = row))

LC_fixed_field_data_all_focal_trees <- tibble()#creating the empty tibble 

#calculating the distances of each tree within the buffer to the focal tree and the competition metric values
for (i in 1:nrow(LC_fixed_field_data_processed_focal)){ #for the length of the buffers with trees inside of them #LC_fixed_field_data_processed_focal_row
  row_num = i
  LC_tree_buffer_inside_df <- as.data.frame(LC_tree_buffer_inside_0) #uses data of non-isolated and isolated focal trees 
  LC_tree_buffer_inside_df_i <- LC_tree_buffer_inside_df %>% 
    filter(row == row_num) #isolate a row of the buffer dataframe
  LC_tree_buffer_inside_sf_i <- st_as_sf(LC_tree_buffer_inside_df_i) #set the row as a simple feature
  all_pts_buffer <- st_contains(LC_tree_buffer_inside_sf_i, LC_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
  LC_fixed_field_data_processed_trees <- LC_fixed_field_data_processed_sf %>%
    filter(X_sequential %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer.
  
  #correct sequence of focal trees
  correct_focal <- LC_fixed_field_data_processed_focal[i,]$X_sequential
  
  #create a dataframe with only the focal tree 
  LC_fixed_field_data_focal_tree <- LC_fixed_field_data_processed_trees %>%
    filter(X_sequential %in% correct_focal) 
  
  #isolating the neighbor tree data
  LC_fixed_field_data_neighbor_trees <- LC_fixed_field_data_processed_trees %>%
    filter(X_sequential %notin% LC_fixed_field_data_focal_tree$X_sequential) #filtering out tree data for the neighbor trees 

  #if there are no neighbors, its sets the sum of the response variable divided by the distance of the tree to the focal tree to 0
  if(nrow(LC_fixed_field_data_neighbor_trees) == 0){
    sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
    sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
    sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
    sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
    sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  } else{
    # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
    LC_fixed_field_data_neighbor_trees <-  LC_fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
      mutate(focal_distance = as.numeric(st_distance(geometry,  LC_fixed_field_data_focal_tree$geometry))) %>% #calculate the distance between the focal tree and each tree that neighbors it
      mutate(focal_distance = case_when(focal_distance == 0 ~ 0.0000016, 
                                        focal_distance != 0 ~ focal_distance)) %>% #replace values of 0 (if the coords are the same for multiple trees) with a value an order of magnitude smaller than the smallest distance in our dataset
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
    for (y in 1:nrow(LC_fixed_field_data_neighbor_trees)){ #adding the size values of each neighbor to a sum total of the neighbors size values
      sum_SCA_over_distance = sum_SCA_over_distance + LC_fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
      sum_LCA_over_distance = sum_LCA_over_distance + LC_fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the LCA of each neighbor
      sum_CA_over_distance = sum_CA_over_distance + LC_fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
      sum_CS_over_distance = sum_CS_over_distance + LC_fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
      sum_DBH_over_distance = sum_DBH_over_distance + LC_fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
    }
  }
  #creating a tibble with all of the calculated sizes over distances and other tree attributes for each focal tree
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  LC_fixed_field_data_focal_tree <- cbind(LC_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  LC_fixed_field_data_all_focal_trees <- rbind(LC_fixed_field_data_all_focal_trees, LC_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
  
}
View(LC_fixed_field_data_neighbor_trees)
View(LC_fixed_field_data_all_focal_trees)

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
LC_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
LC_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)
LC_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test <- model.sel(LC_gls_focal_SCA, LC_gls_focal_SCA_exp, LC_gls_focal_SCA_gaus, LC_gls_focal_SCA_spher, LC_gls_focal_SCA_lin, LC_gls_focal_SCA_ratio)
View(LC_AIC_test)

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
summary(LC_gls_focal_SCA_lin)

#non parametric Mann-Kendall Test for the version without outliers
LC_tau_result_SCA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance, LC_fixed_field_data_all_focal_trees$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_SCA)

# Calculate the trend line
LC_trend_line_SCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_long ~ LC_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

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

View(LC_fixed_field_data_all_focal_trees)

#unlogged version of generlized linear model
LC_gls_focal_LCA <- gls(Canopy_long ~ sum_LCA_over_distance, data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)
LC_gls_focal_LCA_exp <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)
LC_gls_focal_LCA_gaus <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)
LC_gls_focal_LCA_spher <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)
LC_gls_focal_LCA_lin <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)
LC_gls_focal_LCA_ratio <- gls(Canopy_long ~ sum_LCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_LCA_outliers)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_LCA <- model.sel(LC_gls_focal_LCA, LC_gls_focal_LCA_exp, LC_gls_focal_LCA_gaus, LC_gls_focal_LCA_spher, LC_gls_focal_LCA_lin, LC_gls_focal_LCA_ratio)
View(LC_AIC_test_LCA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_LCA_outliers, aes(x= LC_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(LC_fixed_field_data_all_focal_trees_no_LCA_outliers, aes(sample = LC_gls_focal_LCA$residuals))+
  geom_qq()

shapiro.test(LC_gls_focal_LCA$residuals) #shapiro-wilk test, sign so non-normal residuals

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees , aes(x = LC_gls_focal_LCA$fitted, y = LC_gls_focal_LCA$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and LCA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_LCA, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_LCA)
summary(LC_gls_focal_LCA_lin)

#non parametric Mann-Kendall Test for the version without outliers
LC_tau_result_LCA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, LC_fixed_field_data_all_focal_trees$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_LCA)

# Calculate the trend line
LC_trend_line_LCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_long ~ LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
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
LC_gls_focal_CA <- gls(Canopy_area ~ sum_CA_over_distance, data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)
LC_gls_focal_CA_exp <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)
LC_gls_focal_CA_gaus <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)
LC_gls_focal_CA_spher <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)
LC_gls_focal_CA_lin <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)
LC_gls_focal_CA_ratio <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CA_outliers)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CA <- model.sel(LC_gls_focal_CA, LC_gls_focal_CA_exp, LC_gls_focal_CA_gaus, LC_gls_focal_CA_spher, LC_gls_focal_CA_lin, LC_gls_focal_CA_ratio)
View(LC_AIC_test_CA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_CA_outliers, aes(x= LC_gls_focal_CA_lin$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(LC_fixed_field_data_all_focal_trees_no_CA_outliers, aes(sample = LC_gls_focal_CA_lin$residuals))+
  geom_qq()

# shapiro-wilk, not sign so normal residuals for no outliers, and sign for when residuals so we are using a Mann-Kendall non-parametric
shapiro.test(LC_gls_focal_CA_lin$residuals) 

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees_no_CA_outliers , aes(x = LC_gls_focal_CA_lin$fitted, y = LC_gls_focal_CA_lin$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_CA_lin, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_CA_lin)
summary(LC_AIC_test_CA)

#non parametric Mann-Kendall Test for the version without outliers
LC_tau_result_CA <- cor.test(LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, LC_fixed_field_data_all_focal_trees$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CA)

# Calculate the trend line
LC_trend_line_CA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_area ~ LC_fixed_field_data_all_focal_trees$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_area), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CA_over_distance, y = LC_trend_line_CA), color = "red") +
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
LC_gls_focal_CS <- gls(Crown_spread ~ sum_CS_over_distance, data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)
LC_gls_focal_CS_exp <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)
LC_gls_focal_CS_gaus <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)
LC_gls_focal_CS_spher <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)
LC_gls_focal_CS_lin <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)
LC_gls_focal_CS_ratio <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_CS_outliers)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_CS <- model.sel(LC_gls_focal_CS, LC_gls_focal_CS_exp, LC_gls_focal_CS_gaus, LC_gls_focal_CS_ratio) #without linear correlation and without spherical

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_CS_outliers, aes(x= LC_gls_focal_CS$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq norm
ggplot(LC_fixed_field_data_all_focal_trees_no_CS_outliers, aes(sample = LC_gls_focal_CS$residuals))+
  geom_qq()

# shapiro-wilk, n sign for both versions with and without outliers so used mann-kendall non-parametric test
shapiro.test(LC_gls_focal_CA$residuals) # shapiro-wilk, n sign for both versions with and without outliers so used mann-kendall non-parametric test

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees_no_CS_outliers , aes(x = LC_gls_focal_CS$fitted, y = LC_gls_focal_CS$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(LC_gls_focal_CS, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(LC_gls_focal_CS)
summary(LC_AIC_test_CA)

#non parametric Mann-Kendall Test for the version without outliers
LC_tau_result_CS <- cor.test(LC_fixed_field_data_all_focal_trees_no_CS_outliers$sum_CS_over_distance, LC_fixed_field_data_all_focal_trees_no_CS_outliers$Crown_spread,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_CS)

# Calculate the trend line
LC_trend_line_CS <- predict(loess(LC_fixed_field_data_all_focal_trees$Crown_spread ~ LC_fixed_field_data_all_focal_trees$sum_CS_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
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
LC_fixed_field_data_all_focal_trees_no_DBH_outliers <- LC_fixed_field_data_all_focal_trees[-c(3),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
LC_gls_focal_DBH <- gls(DBH_ag ~ sum_DBH_over_distance, data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)
LC_gls_focal_DBH_exp <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)
LC_gls_focal_DBH_gaus <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)
LC_gls_focal_DBH_spher <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)
LC_gls_focal_DBH_lin <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)
LC_gls_focal_DBH_ratio <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers)

#ordering models by which ones have the lowest Akaike information criterion
LC_AIC_test_DHB <- model.sel(LC_gls_focal_DBH, LC_gls_focal_DBH_exp, LC_gls_focal_DBH_lin, LC_gls_focal_DBH_gaus, LC_gls_focal_DBH_spher, LC_gls_focal_DBH_ratio) #without linear correlation

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_DBH_outliers, aes(x= LC_gls_focal_DBH_gaus$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_DBH_outliers, aes(sample = LC_gls_focal_DBH_gaus$residuals))+
  geom_qq()

# shapiro-wilk, significant so non-normal residuals for both
shapiro.test(LC_gls_focal_DBH_gaus$residuals)

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees_no_DBH_outliers, aes(x = LC_gls_focal_DBH$fitted, y = LC_gls_focal_DBH$residuals))+
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
summary(LC_gls_focal_DBH)
summary(LC_gls_focal_DBH_gaus)

#non parametric Mann-Kendall Test for the version without outliers
LC_tau_result_DBH <- cor.test(LC_fixed_field_data_all_focal_trees_no_DBH_outliers$sum_CS_over_distance, LC_fixed_field_data_all_focal_trees_no_DBH_outliers$DBH_ag,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(LC_tau_result_DBH)

# Calculate the trend line
LC_trend_line_DBH <- predict(loess(LC_fixed_field_data_all_focal_trees_no_DBH_outliers$DBH_ag ~ LC_fixed_field_data_all_focal_trees_no_DBH_outliers$sum_DBH_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


# SD

#creating a grid over the points with a 10 m edge buffer

# creating a bounding box based on the point locations
SD_box <- st_bbox(SD_fixed_field_data_processed_sf)

#cropping the tree points down by 20 m on all sides

#creating a cropped bbox 
SD_box_sf <- SD_box %>% #turning the bbox into polygon
  st_as_sfc()
SD_box_spatial <- as(SD_box_sf, 'Spatial') #turning the polygon into a spatial polygon to be able to use raster::crop
SD_box_spatial_cropped <- raster::crop(SD_box_spatial, extent((SD_box[1]+20), (SD_box[3]-20), (SD_box[2]+20), (SD_box[4]-20))) #cropping the xmin, xmax, ymin, and ymax by 20 m inside
SD_box_sf_cropped <-  SD_box_spatial_cropped %>% #turning the spatial polygon into a polygon
  st_as_sfc()

#cropping the points by the cropped box
SD_fixed_field_data_processed_sf_cropped<- st_crop(SD_fixed_field_data_processed_sf, SD_box_sf_cropped)

#Creating a grid over the cropped tree points 
SD_tree_grid_cropped <- st_make_grid(SD_fixed_field_data_processed_sf_cropped, cellsize = (((40*mean(SD_fixed_field_data_processed$DBH_ag))*2)*2))

#plotting the original box, cropped box, original tree points, and cropped tree points
ggplot()+
  geom_sf(data=SD_box_sf)+ #old box
  geom_sf(data=SD_box_sf_cropped)+ #cropped box 
  geom_sf(data=SD_fixed_field_data_processed_sf)+ #original points
  geom_sf(data=SD_fixed_field_data_processed_sf_cropped, color = "red") #old points

#creating an x_sequential column that is 1 through the number of SD points
SD_fixed_field_data_processed_sf <- SD_fixed_field_data_processed_sf %>%
  mutate(X_sequential = 1:nrow(SD_fixed_field_data_processed_sf))
View(SD_fixed_field_data_processed_sf)

#selecting a focal point from each grid cell with trees within them
SD_list_grids_and_points <- st_contains(SD_tree_grid_cropped, SD_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(25) #setting the seed
SD_list_grids_and_focal_trees <- lapply(SD_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
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

View(SD_fixed_field_data_processed_sf)

#creating a dataframe of all of the focal trees with their row number in the overall tree point dataframe and in which grid cell they are in
SD_list_grids_and_focal_trees_df <- as.data.frame(unlist(SD_list_grids_and_focal_trees)) #unlists the list of grid cells and what focal trees were within them and turns it into a dataframe
colnames(SD_list_grids_and_focal_trees_df) <- c("tree_row_num") #changes the column name 
SD_list_grids_and_focal_trees_fixed <- SD_list_grids_and_focal_trees_df %>% #filters out grid cells that do not have trees within them
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree.    #cell_num = row_number()
  mutate(data_row = SD_fixed_field_data_processed_sf$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
SD_fixed_field_data_processed_focal <- SD_fixed_field_data_processed_sf %>%
  filter(X_sequential %in% SD_list_grids_and_focal_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe 

View(SD_fixed_field_data_processed_focal)

#creating the buffer around the focal points
SD_focal_tree_buffers <-st_buffer(SD_fixed_field_data_processed_focal$geometry, 40*mean(SD_fixed_field_data_processed_focal$DBH_ag))

#graphing the selected focal trees, the buffers, the grid
ggplot()+
  geom_sf(data = SD_tree_grid_cropped)+
  geom_sf(data=SD_focal_tree_buffers, color = "blue") +
  geom_sf(data = SD_fixed_field_data_processed) +
  geom_sf(data= SD_fixed_field_data_processed_focal, aes(color = X))

#caSDulating the size/distance for focal trees and neighbors within buffers for buffers with only the focal tree and with more 

#create a tibble with the the number of trees within the buffers that contain trees
SD_tree_buffers_points_within_0 <- st_contains(SD_focal_tree_buffers, SD_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 0) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
SD_tree_buffer_inside_0 <- SD_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% SD_tree_buffers_points_within_0$row) #only keep polygons that match the row number of the grid cells with trees within them 

#Checking that row number in focal dataset is the same as the buffer dataset
SD_fixed_field_data_processed_focal_row <- SD_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
SD_tree_buffer_inside_0 <- mutate(SD_tree_buffer_inside_0, row = as.factor(row)) #making sure the buffers have the same row number as the focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data=SD_tree_buffer_inside_0, aes(color = row))+
  geom_sf(data=SD_fixed_field_data_processed_focal_row, aes(color = row))

#caSDulating the size/distance for focal trees and neighbors within buffers for buffers with more than just the focal tree

#create a tibble with the the number of trees within the buffers that contain trees
SD_tree_buffers_points_within <- st_contains(SD_focal_tree_buffers, SD_fixed_field_data_processed_sf, sparse =F) %>%
  rowSums() %>% #find how many trees are within each grid
  as_tibble() %>% 
  mutate(row = row_number()) %>% #assign a new column with row numbers 
  filter(value > 1) #filter out any buffers with only the focal tree

#filter out the buffers to only have the buffers that contain trees
SD_tree_buffer_inside <- SD_focal_tree_buffers %>%
  st_as_sf() %>% 
  mutate(row = row_number()) %>% #create a column with row numbers
  filter(row %in% SD_tree_buffers_points_within$row) #only keep buffers that match the row number of the buffers cells with trees within them 
View(SD_tree_buffer_inside)

#plotting the points with buffers with neighbors in it and without neighbors, "isolated focal trees"
ggplot()+
  geom_sf(data = SD_focal_tree_buffers)+
  geom_sf(data = SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_fixed_field_data_processed_focal, color = 'blue')

#Checking that row number in focal dataset is the same as the buffer dataset
SD_fixed_field_data_processed_focal_row <- SD_fixed_field_data_processed_focal %>%
  as.data.frame() %>%
  mutate(row = as.factor(row_number())) %>%
  st_as_sf()
SD_tree_buffer_inside <- mutate(SD_tree_buffer_inside, row = as.factor(row)) #making sure the buffers have the same row number as the isoalted focal data

#plotting the grid, the buffers with and without neighbors, and the focal trees, to see if the row numbers for the buffers match the row numbers for the focal tree points
ggplot()+
  geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data=SD_tree_buffer_inside, aes(color = row))+
  geom_sf(data=SD_fixed_field_data_processed_focal_row, aes(color = row))

ggplot()+
 # geom_sf(data = SD_tree_grid_cropped) +
  geom_sf(data = river_SD_trans)+
  #geom_sf(data=SD_tree_buffer_inside)+
  geom_sf(data= SD_fixed_field_data_processed)+
  geom_sf(data=SD_fixed_field_data_processed_focal_row, color = 'red')+
  theme_light()



SD_fixed_field_data_all_focal_trees <- tibble()#creating the empty tibble 

#calculating the distances of each tree within the buffer to the focal tree and the competition metric values
for (i in 1:nrow(SD_fixed_field_data_processed_focal)){ #for the length of the buffers with trees inside of them #SD_fixed_field_data_processed_focal_row
  row_num = i
  SD_tree_buffer_inside_df <- as.data.frame(SD_tree_buffer_inside_0) #uses data of non-isolated and isolated focal trees 
  SD_tree_buffer_inside_df_i <- SD_tree_buffer_inside_df %>% 
    filter(row == row_num) #isolate a row of the buffer dataframe
  SD_tree_buffer_inside_sf_i <- st_as_sf(SD_tree_buffer_inside_df_i) #set the row as a simple feature
  all_pts_buffer <- st_contains(SD_tree_buffer_inside_sf_i, SD_fixed_field_data_processed_sf, sparse = F) #assign true or falses to the trees based on whether they are within that polygon
  possible_pts_buffer <- which(all_pts_buffer == T) #keep only the rows of trees that are within the polygon
  SD_fixed_field_data_processed_trees <- SD_fixed_field_data_processed_sf %>%
    filter(X_sequential %in% possible_pts_buffer) #filtering to the data to only be the trees within the buffer. #before it was X %in% possible_pts_buffer
  
  #create column with correct tree number
  correct_focal <- SD_fixed_field_data_processed_focal[i,]$X_sequential
  
  #create a dataframe with only the focal tree
  SD_fixed_field_data_focal_tree <- SD_fixed_field_data_processed_trees %>%
    filter(X_sequential %in% correct_focal) 
  
  #isolating the neighbor tree data
  SD_fixed_field_data_neighbor_trees <- SD_fixed_field_data_processed_trees %>%
    filter(X_sequential %notin% SD_fixed_field_data_focal_tree$X_sequential) #filtering out tree data for the neighbor trees 
  
  #if there are no neighbors, its sets the sum of the response variable divided by the distance of the tree to the focal tree to 0
  if(nrow(SD_fixed_field_data_neighbor_trees) == 0){
    sum_SCA_over_distance = 0 #create a new variable for short canopy axis over distance to focal tree set to 0
    sum_LCA_over_distance = 0 #create a new variable for long canopy axis over distance to focal tree set to 0
    sum_CA_over_distance = 0 #create a new variable for canopy area over distance to focal tree set to 0
    sum_CS_over_distance = 0 #create a new variable for crown spread over distance to focal tree set to 0
    sum_DBH_over_distance = 0 #create a new variable for DBH over distance to focal tree set to 0
  } else{
    
    # for each neighbor tree, calculates the distance of the tree to the focal tree and find the shape/size metric divided by the distance
    SD_fixed_field_data_neighbor_trees <-  SD_fixed_field_data_neighbor_trees %>% #create a dataframe with only the neighbors of the focal tree
      mutate(focal_distance = as.numeric(st_distance(geometry,  SD_fixed_field_data_focal_tree$geometry))) %>% #caSDulate the distance between the focal tree and each tree that neighbors it
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
    for (y in 1:nrow(SD_fixed_field_data_neighbor_trees)){ #adding the size values of each neighbor to a sum total of the neighbors size values
      sum_SCA_over_distance = sum_SCA_over_distance + SD_fixed_field_data_neighbor_trees$SCA_over_distance[y] #summing the SCA of each neighbor
      sum_LCA_over_distance = sum_LCA_over_distance + SD_fixed_field_data_neighbor_trees$LCA_over_distance[y] #summing the SDA of each neighbor
      sum_CA_over_distance = sum_CA_over_distance + SD_fixed_field_data_neighbor_trees$CA_over_distance[y] #summing the CA of each neighbor
      sum_CS_over_distance = sum_CS_over_distance + SD_fixed_field_data_neighbor_trees$CS_over_distance[y] #summing the CS of each neighbor
      sum_DBH_over_distance = sum_DBH_over_distance + SD_fixed_field_data_neighbor_trees$DBH_over_distance[y] #summing the DBH of each neighbor
    }
  }
  #creating a tibble with all of the calculated sizes over distances and other tree attributes for each focal tree
  all_vals_tibble <- tibble(sum_SCA_over_distance, sum_LCA_over_distance, sum_CS_over_distance, sum_CA_over_distance, sum_DBH_over_distance)
  SD_fixed_field_data_focal_tree <- cbind(SD_fixed_field_data_focal_tree, all_vals_tibble) #bind the sizes over distances values within each buffer to the focal trees
  SD_fixed_field_data_all_focal_trees <- rbind(SD_fixed_field_data_all_focal_trees, SD_fixed_field_data_focal_tree) #add the focal trees with sum of size over distance values to the originally empty tibble
  
  
}
View(SD_fixed_field_data_neighbor_trees)
View(SD_fixed_field_data_all_focal_trees)

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
LC_fixed_field_data_all_focal_trees_no_SCA_outliers <- LC_fixed_field_data_all_focal_trees[-c(3, 24),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_SCA <- gls(Canopy_short ~ sum_SCA_over_distance, data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
SD_gls_focal_SCA_exp <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corExp(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
SD_gls_focal_SCA_gaus <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
SD_gls_focal_SCA_spher <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
SD_gls_focal_SCA_lin <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corLin(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)
SD_gls_focal_SCA_ratio <- gls(Canopy_short ~ sum_SCA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_SCA <- model.sel(SD_gls_focal_SCA, SD_gls_focal_SCA_lin, SD_gls_focal_SCA_exp, SD_gls_focal_SCA_gaus, SD_gls_focal_SCA_spher, SD_gls_focal_SCA_ratio) #without linear
View(SD_AIC_test_SCA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(LC_fixed_field_data_all_focal_trees_no_SCA_outliers, aes(x= SD_gls_focal_SCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. SCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

# qq nrom
ggplot(LC_fixed_field_data_all_focal_trees_no_SCA_outliers, aes(sample = SD_gls_focal_SCA$residuals))+
  geom_qq()

#shapiro-welk test, not sign so normal residuals
shapiro.test(SD_gls_focal_SCA$residuals)

#checking equal variance
ggplot(data = LC_fixed_field_data_all_focal_trees_no_SCA_outliers , aes(x = SD_gls_focal_SCA$fitted, y = SD_gls_focal_SCA_gaus$residuals))+
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

#non parametric Mann-Kendall Test
SD_tau_result_SCA <- cor.test(SD_fixed_field_data_all_focal_trees_no_LCA_outliers$sum_SCA_over_distance, SD_fixed_field_data_all_focal_trees_no_SCA_outliers$Canopy_short,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_SCA)

# Calculate the trend line
SD_trend_line_SCA <- predict(loess(SD_fixed_field_data_all_focal_trees$Canopy_short ~ SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = (SD_fixed_field_data_all_focal_trees$Canopy_short), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_all_focal_trees$sum_SCA_over_distance, y = SD_trend_line_SCA), color = "red") +
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
SD_AIC_test_LCA <- model.sel(SD_gls_focal_LCA, SD_gls_focal_LCA_exp, SD_gls_focal_LCA_gaus, SD_gls_focal_LCA_spher, SD_gls_focal_LCA_lin, SD_gls_focal_LCA_ratio)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_LCA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. LCA over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_LCA$residuals))+
  geom_qq()

#shapiro-wilk test, not sign so normal residuals for with outliers, and sign without outliers
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
summary(SD_gls_focal_LCA_lin)

#non parametric Mann-Kendall Test
SD_tau_result_LCA <- cor.test(SD_fixed_field_data_all_focal_trees_no_LCA_outliers$sum_LCA_over_distance, SD_fixed_field_data_all_focal_trees_no_LCA_outliers$Canopy_long,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_LCA)

# Calculate the trend line
SD_trend_line_LCA <- predict(loess(SD_fixed_field_data_all_focal_trees$Canopy_long ~ SD_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (SD_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = SD_trend_line_CS), color = "red") +
  labs(x = "Sum of LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
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
SD_fixed_field_data_all_focal_trees_no_CA_outliers <- SD_fixed_field_data_all_focal_trees[-c(3,14,23,24),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CA <- gls(Canopy_area ~ sum_CA_over_distance, data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)
SD_gls_focal_CA_exp <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)
SD_gls_focal_CA_gaus <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)
SD_gls_focal_CA_spher <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)
SD_gls_focal_CA_lin <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)
SD_gls_focal_CA_ratio <- gls(Canopy_area ~ sum_CA_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees_no_CA_outliers)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CA <- model.sel(SD_gls_focal_CA, SD_gls_focal_CA_exp, SD_gls_focal_CA_gaus, SD_gls_focal_CA_lin, SD_gls_focal_CA_ratio) #SD_gls_focal_CA_spher
View(SD_AIC_test_CA)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees_no_CA_outliers, aes(x= SD_gls_focal_CA$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Canopy Area over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees_no_CA_outliers, aes(sample = SD_gls_focal_CA$residuals))+
  geom_qq()

# shapiro-wilk, for both versions with and without outliers sign so residuals non-normal
shapiro.test(SD_gls_focal_CA$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees, aes(x = SD_gls_focal_CA_ratio$fitted, y = SD_gls_focal_CA_ratio$residuals))+
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
summary(SD_gls_focal_CA_ratio)
summary(SD_AIC_test_CA)

#non parametric Mann-Kendall Test
SD_tau_result_CA <- cor.test(SD_fixed_field_data_all_focal_trees_no_CA_outliers$sum_CA_over_distance, SD_fixed_field_data_all_focal_trees_no_CA_outliers$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CA)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_all_focal_trees_no_CA_outliers$Canopy_area ~ SD_fixed_field_data_all_focal_trees_no_CA_outliers$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (SD_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = SD_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
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
  
View(SD_fixed_field_data_all_focal_trees)
#Cook's D
SD_lm_focal_CS <- lm(Crown_spread ~ sum_CS_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_lm_focal_CS_cooks <- cooks.distance(SD_lm_focal_CS) #calculating the cook.s D for each point
plot(SD_lm_focal_CS_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- SD_lm_focal_CS_cooks[(SD_lm_focal_CS_cooks > (3 * mean(SD_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential
SD_fixed_field_data_all_focal_trees_no_CS_outliers <- SD_fixed_field_data_all_focal_trees[-c(3),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_CS <- gls(Crown_spread ~ sum_CS_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_exp <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_gaus <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_spher <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_lin <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_CS_ratio <- gls(Crown_spread ~ sum_CS_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_CS <- model.sel(SD_gls_focal_CS, SD_gls_focal_CS_exp, SD_gls_focal_CS_gaus, SD_gls_focal_CS_spher, SD_gls_focal_CS_ratio) #without linear correlation
View(SD_AIC_test_CS)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees_no_CS_outliers, aes(x= SD_gls_focal_CS_lin$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Crown Spread over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees_no_CS_outliers, aes(sample = SD_gls_focal_CS_lin$residuals))+
  geom_qq()

# shapiro-wilk, not signficant, meaning not signficantly different from normal
shapiro.test(SD_gls_focal_CS$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees_no_CS_outliers , aes(x = SD_gls_focal_CS_lin$fitted, y = SD_gls_focal_CS$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and CA over Distance")

fligner.test(Canopy_short ~ sum_CS_over_distance, data = SD_fixed_field_data_all_focal_trees_no_CS_outliers) #checks if residuals are normal, even if residuals not normal and having difficulty removing outliers, null is group variances are equal

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_CS_lin, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_CS)
summary(SD_gls_focal_CS_gaus)

#non parametric Mann-Kendall Test
SD_tau_result_CA <- cor.test(SD_fixed_field_data_all_focal_trees_no_CA_outliers$sum_CA_over_distance, SD_fixed_field_data_all_focal_trees_no_CA_outliers$Canopy_area,  method = "kendall")

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CA)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_all_focal_trees_no_CA_outliers$Canopy_area ~ SD_fixed_field_data_all_focal_trees_no_CA_outliers$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (SD_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = SD_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
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
SD_fixed_field_data_all_focal_trees_no_DBH_outliers <- SD_fixed_field_data_all_focal_trees[-c(3,23,24),]

#creating generalized linear model with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratices)
SD_gls_focal_DBH <- gls(DBH_ag ~ sum_DBH_over_distance, data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_exp <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corExp(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_gaus <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corGaus(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_spher <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corSpher(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_lin <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corLin(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)
SD_gls_focal_DBH_ratio <- gls(DBH_ag ~ sum_DBH_over_distance, correlation = corRatio(form = ~X.1 + Y), data = SD_fixed_field_data_all_focal_trees)

#ordering models by which ones have the lowest Akaike information criterion
SD_AIC_test_DHB <- model.sel(SD_gls_focal_DBH, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_exp, SD_gls_focal_DBH_lin, SD_gls_focal_DBH_gaus, SD_gls_focal_DBH_spher, SD_gls_focal_DBH_ratio) 
View(SD_AIC_test_CS)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(SD_fixed_field_data_all_focal_trees, aes(x= SD_gls_focal_DBH_ratio$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. DBH over Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qq norm
ggplot(SD_fixed_field_data_all_focal_trees, aes(sample = SD_gls_focal_DBH_ratio$residuals))+
  geom_qq()

# shapiro-wilk, not significant so normal
shapiro.test(SD_gls_focal_DBH_ratio$residuals) 

#checking equal variance
ggplot(data = SD_fixed_field_data_all_focal_trees , aes(x = SD_gls_focal_DBH_ratio$fitted, y = SD_gls_focal_DBH_ratio$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and DBH over Distance")

#plotting semivariogram, checking we have appropriately removed the spatial autocorrelation 
#(hovering around 1 indicates model controlled for spatial autocorrelation)
semivario <- Variogram(SD_gls_focal_DBH_ratio, form = ~X.1 + Y, resType = "normalized")
plot(semivario, smooth = TRUE)

#Slope Test visible in summary of the lm, lack of significant of slope indicates lack of impact from competition
#positive slope hints at facilitation
#negative slope hints at competition
summary(SD_gls_focal_DBH_ratio)
summary(SD_gls_focal_DBH_gaus)

# Print Kendall's tau and its associated p-value
print(SD_tau_result_CA)

# Calculate the trend line
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_all_focal_trees_no_CA_outliers$Canopy_area ~ SD_fixed_field_data_all_focal_trees_no_CA_outliers$sum_CA_over_distance))

# Create a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (SD_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = SD_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

