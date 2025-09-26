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

# The script is broken into sections of 
# 1) loading and processing the packages and spatial/size/shape data for the trees in the Las Matancitas,
#San Dionisio, and La Cobriza populations and loading in the river outline shapefiles, 
# 2) creating a dataframe with the coordinates and average distance between each tree and their 5 nearest neighbors
# 3) graphing descriptive summary histograms and calculating summary statistics (mean, median, sd, etc.) for each response 
#variable (SCA, LCA, CA, CS, DBH),
# 4) Running the global and local Moran's I analyses, 


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

#adding a sequential column, "X," to number each tree

fixed_field_data_processed <- fixed_field_data_processed %>%
  mutate(X = row_number())

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


# Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns

LM_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_sf_trans_coordinates %>%
  filter(Locality == "SD")



#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_sf_trans_coordinates) + # Generate the base plot
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
duplicates <- fixed_field_data_processed_sf_trans_coordinates %>% #creates a dataframe called duplicates that filters out X.1 and Y if they have duplicates
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

#creating function to compute the Global and Local Moran's I

morans_I <- function(population, variable){
  #assigning the size/shape metric to a variable
  metric <- variable
  
  #assigning the dataframes based on the population
  if (population == "LM") {
    dataframe <- LM_fixed_field_data_processed
  }
  
  if (population == "LC") {
    dataframe <- LC_fixed_field_data_processed
  }
  
  if (population == "SD") {
    dataframe <- SD_fixed_field_data_processed
  }
  
  ## Global Moran's I
  
  #creating a matrix of the tree locations
  tree.coord.matrix <- as.matrix(cbind(dataframe$X.1, 
                                       dataframe$Y))
  
  #creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
  knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$DBH_ag)))
  
  #inverse distance weighting with raw distance-based weights without applying any normalization
  lw.dist <- nb2listwdist(knn.dist, fixed_field_data_processed_sf_trans_coordinates, type="idw", style="W", 
                          alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
  
  #creating lags for each tree, which computes the average neighboring short canopy axis for each tree
  dataframe$lag.metric <- lag.listw(lw.dist, dataframe[[metric]])
  
  # Create a regression model of the lagged response variable (average amongst closest trees) vs. the known response variable 
  lm <- lm(lag.metric ~ dataframe[[metric]], dataframe)
  
  #computing the Moran's I statistic
  global.moran.I <- moran(dataframe[[metric]], listw = lw.dist, n = length(lw.dist$neighbours), S0 = Szero(lw.dist))
  global.moran.I
  
  #assessing statistical significance with a Monte-Carlo simulation
  MC.LM.metric <- moran.mc(dataframe[[metric]], lw.dist, nsim = 999)
  MC.LM.metric
  
  #plot of simulated Moran's I values against our value
  plot(MC.LM.metric, main="", las=1, xlab = metric)
  MC.LM.metric$p.value #extracting the pvalue
  
  print(paste0("Global Moran's I: ", global.moran.I$I))
  print(paste0("Global Moran's I Monte Carlo P-value: ", MC.LM.metric$p.value))
  
  
  ## Local Moran's I
  
  #using the weighted neighbors to simulate size values at random
  MC_local <- localmoran_perm(dataframe[[metric]], lw.dist, nsim = 9999, alternative = "greater")
  MC_local.df <- as.data.frame(MC_local)
  
  #calculating the p-values for each individual tree Moran's I, observed vs. expected
  dataframe$p.metric  <- MC_local.df$`Pr(folded) Sim`
  #adjusting the p-vlaues to take into account multiple tests
  dataframe$p.metric.adjusted <- p.adjust(dataframe$p.metric, 
                                          method = "fdr", n=length(dataframe$p.metric))
  
  #Number of trees with significant adjusted Moran's I
  sig.tree <- length(which(dataframe$p.metric.adjusted < 0.05))
  
  #filtering out significant p-values
  dataframe_sign <- dataframe %>%
    mutate(pval_sig = p.metric.adjusted <= .05) %>%
    filter(pval_sig == T)
  
  return(list(dataframe$lag.metric, lm, MC.LM.metric,
              MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I))
  
}


###Test for LM###

#creating a shapefile for plotting later
LM_box <- st_bbox(river_LM_trans)

#Short Canopy Axis

#global Moran's I

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#conducting Moran's I analysis
LM_SCA_Morans_I <- morans_I("LM", "Canopy_short")

#regression for ANN size metric vs. tree metric 
LM_SCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LM_SCA_Morans_I[[3]]

#creating a column for the lagged size metric
LM_fixed_field_data_processed$lag.canopy.short <- LM_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LM.canopy.short.df <- LM_SCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LM_SCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LM_fixed_field_data_processed_sign <- LM_SCA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LM_fixed_field_data_processed$p.canopy.short.adjusted <- LM_SCA_Morans_I[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting the local Moran's I 
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


###Long Canopy Axis

#global Moran's I

#conducting Moran's I analysis
LM_LCA_Morans_I <- morans_I("LM", "Canopy_long")

#regression for ANN size metric vs. tree metric 
LM_LCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LM_LCA_Morans_I[[3]]

#creating a column for the lagged size metric
LM_fixed_field_data_processed$lag.canopy.long <- LM_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")


#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LM.canopy.long.df <- LM_LCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LM_LCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LM_fixed_field_data_processed_sign <- LM_LCA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LM_fixed_field_data_processed$p.canopy.long.adjusted <- LM_LCA_Morans_I[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#plotting the local Moran's I 
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#conducting Moran's I analysis
LM_CS_Morans_I <- morans_I("LM", "Crown_spread")

#regression for ANN size metric vs. tree metric 
LM_CS_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LM_CS_Morans_I[[3]]

#creating a column for the lagged size metric
LM_fixed_field_data_processed$lag.crown.spread <- LM_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LM.crown.spread.df <- LM_CS_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LM_CS_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LM_fixed_field_data_processed_sign <- LM_CS_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LM_fixed_field_data_processed$p.crown.spread.adjusted <- LM_CS_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
ggplot(data=MC_local.LM.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#conducting Moran's I analysis
LM_CA_Morans_I <- morans_I("LM", "Canopy_area")

#regression for ANN size metric vs. tree metric 
LM_CA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LM_CA_Morans_I[[3]]

#creating a column for the lagged size metric
LM_fixed_field_data_processed$lag.canopy.area <- LM_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LM.canopy.area.df <- LM_CA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LM_CA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LM_fixed_field_data_processed_sign <- LM_CA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LM_fixed_field_data_processed$p.canopy.area.adjusted <- LM_CA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#plotting local Moran's I
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#conducting Moran's I analysis
LM_DBH_Morans_I <- morans_I("LM", "DBH_ag")

#regression for ANN size metric vs. tree metric 
LM_DBH_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LM_DBH_Morans_I[[3]]

#creating a column for the lagged size metric
LM_fixed_field_data_processed$lag.dbh.ag <- LM_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LM.dbh.ag.df <- LM_DBH_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LM_DBH_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LM_fixed_field_data_processed_sign <- LM_DBH_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LM_fixed_field_data_processed$p.dbh.ag.adjusted <- LM_DBH_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#plotting local Moran's I
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for DBH")


###Test for LC###

#making a shapefile for later plotting
LC_box <- st_bbox(river_LC_trans)

#Short Canopy Axis

#global Moran's I

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#conducting Moran's I analysis
LC_SCA_Morans_I <- morans_I("LC", "Canopy_short")

#regression for ANN size metric vs. tree metric 
LC_SCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LC_SCA_Morans_I[[3]]

#creating a column for the lagged size metric
LC_fixed_field_data_processed$lag.canopy.short <- LC_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LC.canopy.short.df <- LC_SCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LC_SCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LC_fixed_field_data_processed_sign <- LC_SCA_Morans_I[[6]]
 
#assigning the p-values of the adjusted local Moran's I to a dataframe
LC_fixed_field_data_processed$p.canopy.short.adjusted <- LC_SCA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting local Moran's I
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#conducting Moran's I analysis
LC_LCA_Morans_I <- morans_I("LC", "Canopy_long")

#regression for ANN size metric vs. tree metric 
LC_LCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LC_LCA_Morans_I[[3]]

#creating a column for the lagged size metric
LC_fixed_field_data_processed$lag.canopy.long <- LC_LCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LC.canopy.long.df <- LC_LCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LC_LCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LC_fixed_field_data_processed_sign <- LC_LCA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LC_fixed_field_data_processed$p.canopy.long.adjusted <- LC_LCA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#plotting local Moran's I
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#conducting Moran's I analysis
LC_CS_Morans_I <- morans_I("LC", "Crown_spread")

#regression for ANN size metric vs. tree metric 
LC_CS_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LC_CS_Morans_I[[3]]

#creating a column for the lagged size metric
LC_fixed_field_data_processed$lag.crown.spread <- LC_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LC.crown.spread.df <- LC_CS_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LC_CS_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LC_fixed_field_data_processed_sign <- LC_CS_Morans_I[[6]]
 
#assigning the p-values of the adjusted local Moran's I to a dataframe
LC_fixed_field_data_processed$p.crown.spread.adjusted <- LC_CS_Morans_I[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#plotting the Local Moran's I
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#conducting Moran's I analysis
LC_CA_Morans_I <- morans_I("LC", "Canopy_area")

#regression for ANN size metric vs. tree metric 
LC_CA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LC_CA_Morans_I[[3]]

#creating a column for the lagged size metric
LC_fixed_field_data_processed$lag.canopy.area <- LC_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LC.canopy.area.df <- LC_CA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LC_CA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LC_fixed_field_data_processed_sign <- LC_CA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LC_fixed_field_data_processed$p.canopy.area.adjusted <- LC_CA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#plotting Local Moran's I
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#conducting Moran's I analysis
LC_DBH_Morans_I <- morans_I("LC", "DBH_ag")

#regression for ANN size metric vs. tree metric 
LC_DBH_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
LC_DBH_Morans_I[[3]]

#creating a column for the lagged size metric
LC_fixed_field_data_processed$lag.dbh.ag <- LC_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.LC.dbh.ag.df <- LC_DBH_Morans_I[[4]]

#number of trees with Significant Local Moran's I
LC_DBH_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
LC_fixed_field_data_processed_sign <- LC_DBH_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
LC_fixed_field_data_processed$p.dbh.ag.adjusted <- LC_DBH_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#plotting local Moran's I
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

#Test for SD

#making a shapefile for later plotting
SD_box <- st_bbox(river_SD_trans)

#Short Canopy Axis

#global Moran's I

#Moran's I and Monte Carlo, using Lags, requires package: spdep

#conducting Moran's I analysis
SD_SCA_Morans_I <- morans_I("SD", "Canopy_short")

#regression for ANN size metric vs. tree metric 
SD_SCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
SD_SCA_Morans_I[[3]]

#creating a column for the lagged size metric
SD_fixed_field_data_processed$lag.canopy.short <- SD_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.SD.canopy.short.df <- SD_SCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
SD_SCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
SD_fixed_field_data_processed_sign <- SD_SCA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
SD_fixed_field_data_processed$p.canopy.short.adjusted <- SD_SCA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#plotting Local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#conducting Moran's I analysis
SD_LCA_Morans_I <- morans_I("SD", "Canopy_long")

#regression for ANN size metric vs. tree metric 
SD_LCA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
SD_LCA_Morans_I[[3]]

#creating a column for the lagged size metric
SD_fixed_field_data_processed$lag.canopy.long <- SD_LCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.SD.canopy.long.df <- SD_LCA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
SD_LCA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
SD_fixed_field_data_processed_sign <- SD_LCA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
SD_fixed_field_data_processed$p.canopy.long.adjusted <- SD_LCA_Morans_I[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#plotting Local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#conducting Moran's I analysis
SD_CS_Morans_I <- morans_I("SD", "Crown_spread")

#regression for ANN size metric vs. tree metric 
SD_CS_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
SD_CS_Morans_I[[3]]

#creating a column for the lagged size metric
SD_fixed_field_data_processed$lag.crown.spread <- SD_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")+
  theme(text = element_text(size = 20))

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.SD.crown.spread.df <- SD_CS_Morans_I[[4]]

#number of trees with Significant Local Moran's I
SD_CS_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
SD_fixed_field_data_processed_sign <- SD_CS_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
SD_fixed_field_data_processed$p.crown.spread.adjusted <- SD_CS_Morans_I[[7]]

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#conducting Moran's I analysis
SD_CA_Morans_I <- morans_I("SD", "Canopy_area")

#regression for ANN size metric vs. tree metric 
SD_CA_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
SD_CA_Morans_I[[3]]

#creating a column for the lagged size metric
SD_fixed_field_data_processed$lag.canopy.area <- SD_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.SD.canopy.area.df <- SD_CA_Morans_I[[4]]

#number of trees with Significant Local Moran's I
SD_CA_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
SD_fixed_field_data_processed_sign <- SD_CA_Morans_I[[6]]

#assigning the p-values of the adjusted local Moran's I to a dataframe
SD_fixed_field_data_processed$p.canopy.area.adjusted <- SD_CA_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#conducting Moran's I analysis
SD_DBH_Morans_I <- morans_I("SD", "DBH_ag")

#regression for ANN size metric vs. tree metric 
SD_DBH_Morans_I[[2]]

#Monte Carlo Simulation for Global Moran's I
SD_DBH_Morans_I[[3]]

plot(SD_DBH_Morans_I[[3]], col = "red", xlab = "DBH")

#creating a column for the lagged size metric
SD_fixed_field_data_processed$lag.dbh.ag <- SD_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#Local Moran's I 

#assigning a Monta Carlo dataframe for plotting
MC_local.SD.dbh.ag.df <- SD_DBH_Morans_I[[4]]

#number of trees with Significant Local Moran's I
SD_DBH_Morans_I[[5]]

#assigning the trees with the significant local Moran's I to a dataframe
SD_fixed_field_data_processed_sign <- SD_DBH_Morans_I[[6]]
 
#assigning the p-values of the adjusted local Moran's I to a dataframe
SD_fixed_field_data_processed$p.dbh.ag.adjusted <- SD_DBH_Morans_I[[7]]

##Ii is local Moran's I statistic, E.Ii is expected local Moran's I statistic, Vari.Ii is variance of local Moran's I statistic, Z. Ii standard deviation of local Moran's I statistic  
#plotting the local Moran's I values vs. the expected
ggplot(data=MC_local.SD.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for DBH", fill = "Significant P-Value") 

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  labs(color = "Adjusted P Value for DBH", fill = "Significant P-Value") +
  coord_sf(xlim = c(SD_fixed_field_data_processed_sign$X.1-50, SD_fixed_field_data_processed_sign$X.1+50), 
           ylim = c(SD_fixed_field_data_processed_sign$Y-50, SD_fixed_field_data_processed_sign$Y+50)) 
  
