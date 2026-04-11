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

# NOTE: Uncomment and run line 41, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

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

# loading in the processed tree data 
# NOTE: Uncomment and run line 41, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")

# Make a function that is the opposite of the %in% function
`%notin%` <- Negate(`%in%`) 

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
  
  #removing NAs in DBH (needed to compute nearest neighbor and the metric of focus) ADDED AFTER ADDITION OF NEW DATA IN 2026
  dataframe <- dataframe %>%
    filter(!is.na(DBH_ag)) %>%
    filter(!is.na(.data[[metric]]))
  
  #creating a matrix of the tree locations
  tree.coord.matrix <- as.matrix(cbind(dataframe$X.1, 
                                       dataframe$Y))
  
  #creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
  knn.dist <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(dataframe$DBH_ag)))
  
  #inverse distance weighting with raw distance-based weights without applying any normalization
  lw.dist <- nb2listwdist(knn.dist, dataframe, type="idw", style="W", #fixed_field_data_processed_sf_trans_coordinates
                          alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T) # had to set zero.policy to true because of empty neighbor sets
  
  #creating lags for each tree, which computes the average neighboring size metric for each tree
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
  
  #setting a seed for all of the results from using the morans_I() function because localmoran_perm() uses random permutations
  set.seed(25)
  
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
              MC_local.df, sig.tree, dataframe_sign, dataframe$p.metric.adjusted, global.moran.I, dataframe))
  
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

#storing the dataframe used for the Moran's Is with the NAs removed
LM_SCA_dataframe <- LM_SCA_Morans_I[[9]]


#creating a column for the lagged size metric
LM_SCA_dataframe$lag.canopy.short <- LM_SCA_Morans_I[[1]] #LM_fixed_field_data_processed

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_SCA_dataframe, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = "lm", col="blue")+
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
LM_SCA_dataframe$p.canopy.short.adjusted <- LM_SCA_Morans_I[[7]]

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
  geom_sf(data =LM_SCA_dataframe, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for SCA")

#attempting to zoom on the sizes of the significant point
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_SCA_dataframe, aes(size = Canopy_short)) +
  geom_sf(data =LM_SCA_dataframe, aes(color = p.canopy.short.adjusted))+
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

#storing the dataframe used for the Moran's Is with the NAs removed
LM_LCA_dataframe <- LM_LCA_Morans_I[[9]]

#creating a column for the lagged size metric
LM_LCA_dataframe$lag.canopy.long <- LM_LCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_LCA_dataframe, aes(x=Canopy_long, y=lag.canopy.long))+
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
LM_LCA_dataframe$p.canopy.long.adjusted <- LM_LCA_Morans_I[[7]]

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
  geom_sf(data =LM_LCA_dataframe, aes(color = p.canopy.long.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LM_CS_dataframe <- LM_CS_Morans_I[[9]]

#creating a column for the lagged size metric
LM_CS_dataframe$lag.crown.spread <- LM_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_CS_dataframe, aes(x=Crown_spread, y=lag.crown.spread))+
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
LM_CS_dataframe$p.crown.spread.adjusted <- LM_CS_Morans_I[[7]]

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
  geom_sf(data =LM_CS_dataframe, aes(color = p.crown.spread.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LM_CA_dataframe <- LM_CA_Morans_I[[9]]

#creating a column for the lagged size metric
LM_CA_dataframe$lag.canopy.area <- LM_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_CA_dataframe, aes(x=Canopy_area, y=lag.canopy.area))+
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
LM_CA_dataframe$p.canopy.area.adjusted <- LM_CA_Morans_I[[7]]

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
  geom_sf(data =LM_CA_dataframe, aes(color = p.canopy.area.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LM_DBH_dataframe <- LM_DBH_Morans_I[[9]]

#creating a column for the lagged size metric
LM_DBH_dataframe$lag.dbh.ag <- LM_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LM_DBH_dataframe, aes(x=DBH_ag, y=lag.dbh.ag))+
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
LM_DBH_dataframe$p.dbh.ag.adjusted <- LM_DBH_Morans_I[[7]]

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
  geom_sf(data =LM_DBH_dataframe, aes(color = p.dbh.ag.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LC_SCA_dataframe <- LC_SCA_Morans_I[[9]]

#creating a column for the lagged size metric
LC_SCA_dataframe$lag.canopy.short <- LC_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_SCA_dataframe, aes(x=Canopy_short, y=lag.canopy.short))+
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
LC_SCA_dataframe$p.canopy.short.adjusted <- LC_SCA_Morans_I[[7]]

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
  geom_sf(data =LC_SCA_dataframe, aes(color = p.canopy.short.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LC_LCA_dataframe <- LC_LCA_Morans_I[[9]]

#creating a column for the lagged size metric
LC_LCA_dataframe$lag.canopy.long <- LC_LCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_LCA_dataframe, aes(x=Canopy_long, y=lag.canopy.long))+
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
LC_LCA_dataframe$p.canopy.long.adjusted <- LC_LCA_Morans_I[[7]]

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
  geom_sf(data =LC_LCA_dataframe, aes(color = p.canopy.long.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LC_CS_dataframe <- LC_CS_Morans_I[[9]]

#creating a column for the lagged size metric
LC_CS_dataframe$lag.crown.spread <- LC_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_CS_dataframe, aes(x=Crown_spread, y=lag.crown.spread))+
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
LC_CS_dataframe$p.crown.spread.adjusted <- LC_CS_Morans_I[[7]]

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
  geom_sf(data =LC_CS_dataframe, aes(color = p.crown.spread.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LC_CA_dataframe <- LC_CA_Morans_I[[9]]

#creating a column for the lagged size metric
LC_CA_dataframe$lag.canopy.area <- LC_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_CA_dataframe, aes(x=Canopy_area, y=lag.canopy.area))+
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
LC_CA_dataframe$p.canopy.area.adjusted <- LC_CA_Morans_I[[7]]

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
  geom_sf(data =LC_CA_dataframe, aes(color = p.canopy.area.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
LC_DBH_dataframe <- LC_DBH_Morans_I[[9]]

#creating a column for the lagged size metric
LC_DBH_dataframe$lag.dbh.ag <- LC_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=LC_DBH_dataframe, aes(x=DBH_ag, y=lag.dbh.ag))+
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
LC_DBH_dataframe$p.dbh.ag.adjusted <- LC_DBH_Morans_I[[7]]

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
  geom_sf(data =LC_DBH_dataframe, aes(color = p.dbh.ag.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
SD_SCA_dataframe <- SD_SCA_Morans_I[[9]]

#creating a column for the lagged size metric
SD_SCA_dataframe$lag.canopy.short <- SD_SCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_SCA_dataframe, aes(x=Canopy_short, y=lag.canopy.short))+
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
SD_SCA_dataframe$p.canopy.short.adjusted <- SD_SCA_Morans_I[[7]]

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
  geom_sf(data =SD_SCA_dataframe, aes(color = p.canopy.short.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
SD_LCA_dataframe <- SD_LCA_Morans_I[[9]]

#creating a column for the lagged size metric
SD_LCA_dataframe$lag.canopy.long <- SD_LCA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_LCA_dataframe, aes(x=Canopy_long, y=lag.canopy.long))+
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
SD_LCA_dataframe$p.canopy.long.adjusted <- SD_LCA_Morans_I[[7]]

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
  geom_sf(data =SD_LCA_dataframe, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
SD_CS_dataframe <- SD_CS_Morans_I[[9]]

#creating a column for the lagged size metric
SD_CS_dataframe$lag.crown.spread <- SD_CS_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_CS_dataframe, aes(x=Crown_spread, y=lag.crown.spread))+
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
SD_CS_dataframe$p.crown.spread.adjusted <- SD_CS_Morans_I[[7]]

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
  geom_sf(data =SD_CS_dataframe, aes(color = p.crown.spread.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
SD_CA_dataframe <- SD_CA_Morans_I[[9]]

#creating a column for the lagged size metric
SD_CA_dataframe$lag.canopy.area <- SD_CA_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_CA_dataframe, aes(x=Canopy_area, y=lag.canopy.area))+
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
SD_CA_dataframe$p.canopy.area.adjusted <- SD_CA_Morans_I[[7]]

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
  geom_sf(data =SD_CA_dataframe, aes(color = p.canopy.area.adjusted)) +
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

#storing the dataframe used for the Moran's Is with the NAs removed
SD_DBH_dataframe <- SD_DBH_Morans_I[[9]]

plot(SD_DBH_Morans_I[[3]], col = "red", xlab = "DBH")

#creating a column for the lagged size metric
SD_DBH_dataframe$lag.dbh.ag <- SD_DBH_Morans_I[[1]]

# Plot the lagged response variable (average amongst closest trees) vs. the variable 
# positive slope, positive spatial autocorrelation, bigger trees are closer together and smaller trees are closer together
# negative slope, negative spatial autocorrelation, variation in size of trees close together
ggplot(data=SD_DBH_dataframe, aes(x=DBH_ag, y=lag.dbh.ag))+
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
SD_DBH_dataframe$p.dbh.ag.adjusted <- SD_DBH_Morans_I[[7]]

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
  geom_sf(data =SD_DBH_dataframe, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for DBH", fill = "Significant P-Value") 

#plotting the local Moran's I
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_DBH_dataframe, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed_sign, color = "red", aes(fill = "red")) +
  labs(color = "Adjusted P Value for DBH", fill = "Significant P-Value") +
  coord_sf(xlim = c(min(SD_fixed_field_data_processed_sign$X.1)-50, max(SD_fixed_field_data_processed_sign$X.1)-100), 
           ylim = c(min(SD_fixed_field_data_processed_sign$Y)+300, max(SD_fixed_field_data_processed_sign$Y))) 
  

#### Looking at Fruiting ####

# Fruiting is a binary, categorical variable, so we cannot use Moran's I and instead to check for spatial autocorrelation
# of fruiting we will use join count statistics. 

# Join Count Statistic tests whether a binary map (categories A/B) shows spatial clustering

###Test for LM###

#removing the NAs and turning the fruiting variable into a factor
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  filter(!is.na(DBH_ag)) %>%
  filter(!is.na(fruiting)) %>%
  mutate(fruiting = as.factor(fruiting))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(LM_fixed_field_data_processed$X.1, 
                                     LM_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.LM <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(LM_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist.LM <- nb2listwdist(knn.dist.LM, LM_fixed_field_data_processed, type="idw", style="W", 
             alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T)

#running the joint count function
joincount.test(fx = LM_fixed_field_data_processed$fruiting, listw = lw.dist.LM)

#running the monte carlo, permutation test for join count statistics 
set.seed(42) #need to set a seed first before the permutation test
joincount.mc(fx = LM_fixed_field_data_processed$fruiting, listw = lw.dist.LM, nsim = 1000, zero.policy=T,
             alternative="greater", spChk=NULL)

#plotting the fruiting trees
ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = fruiting)) +
  labs(color = "Fruiting")

###Test for LC###

#removing the NAs and turning the fruiting variable into a factor
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  filter(!is.na(DBH_ag)) %>%
  filter(!is.na(fruiting)) %>%
  mutate(fruiting = as.factor(fruiting))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(LC_fixed_field_data_processed$X.1, 
                                     LC_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.LC <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(LC_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist.LC <- nb2listwdist(knn.dist.LC, LC_fixed_field_data_processed, type="idw", style="W", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T)

#running the joint count function
joincount.test(fx = LC_fixed_field_data_processed$fruiting, listw = lw.dist.LC)

#running the monte carlo, permutation test for join count statistics 
set.seed(42) #need to set a seed first before the permutation test
joincount.mc(fx = LC_fixed_field_data_processed$fruiting, listw = lw.dist.LC, nsim = 1000, zero.policy=T,
             alternative="two.sided", spChk=NULL)

#plotting the fruiting trees
ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = fruiting)) +
  labs(color = "Fruiting")

###Test for LM###

#removing the NAs and turning the fruiting variable into a factor
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  filter(!is.na(DBH_ag)) %>%
  filter(!is.na(fruiting)) %>%
  mutate(fruiting = as.factor(fruiting))

#creating a matrix of the tree locations
tree.coord.matrix <- as.matrix(cbind(SD_fixed_field_data_processed$X.1, 
                                     SD_fixed_field_data_processed$Y))

#creates nearest neighbor matrix of the tree coordinates within 40 meters of the mean DBH of the population
knn.dist.SD <- dnearneigh(tree.coord.matrix, d1 = 0, d2 = (40*mean(SD_fixed_field_data_processed$DBH_ag)))

#inverse distance weighting with raw distance-based weights without applying any normalization
lw.dist.SD <- nb2listwdist(knn.dist.SD, SD_fixed_field_data_processed, type="idw", style="W", 
                        alpha = 1, dmax = NULL, longlat = NULL, zero.policy=T)

#running the joint count function
joincount.test(fx = SD_fixed_field_data_processed$fruiting, listw = lw.dist.SD)

#running the monte carlo, permutation test for join count statistics 
set.seed(42) #need to set a seed first before the permutation test
joincount.mc(fx = SD_fixed_field_data_processed$fruiting, listw = lw.dist.SD, nsim = 1000, zero.policy=T,
             alternative="two.sided", spChk=NULL)

#plotting the fruiting trees
ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = fruiting)) +
  labs(color = "Fruiting")


#### Ch-square Test Comparing Fruiting Proportions ####


#creating a contingency table
cont_table <- fixed_field_data_processed_sf_transformed %>%
  st_drop_geometry() %>%
  group_by(locality) %>%
  summarise(is_fruiting = sum(fruiting == "Y", na.rm = T),
            not_fruiting = sum(fruiting == "N", na.rm = T)) %>%
  mutate(total = is_fruiting + not_fruiting) 

#turning the contingency table into a matrix suitable for the chi square testing function
chi_matrix <- cont_table %>%
  column_to_rownames(var = "locality") %>%
  as.matrix()

# Performing the Chi-Square Test
chi_result <- chisq.test(chi_matrix)
print(chi_result)

# View expected frequencies (check if > 5)
chi_result$expected

# View observed frequencies
chi_result$observed

# Running post hoc test with Bonferroni correction
library(chisq.posthoc.test)
chisq.posthoc.test(chi_matrix, method = "bonferroni")

ggplot(fixed_field_data_processed_sf_transformed) +
  geom_bar(aes(x = locality, fill = fruiting), position = "fill")

#looking at the pairwise differences
chisq.test(chi_matrix[c(1, 2),], correct = FALSE)
chisq.test(chi_matrix[c(1, 3),], correct = FALSE)
chisq.test(chi_matrix[c(2, 3),], correct = FALSE)



