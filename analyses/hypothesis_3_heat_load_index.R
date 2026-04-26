# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Investigating to see if Quercus brandegeei tree shapes are influenced by Water Availability Proxies%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by their access to water using generalized additive models (GAMs)

# We used Generalized Additive Models (GAMs) after having trouble with multiple linear 
# regressions because of issues with the normality conditions and some nervousness about linearity. 
# We used four water availability proxies (distance to river, TWI, heat load index, and elevation).

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
# Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
# Extracting and processing slope, elevation, and aspect (4 and 8 cardinal directions) data using 15 m res rasters,
# Extracting the distance to the river of each tree for each population,
# 2) Creating Generalized Additive Models to look at the relationships between the shape/size of trees and their aspect, elevation, slope, and distance to river

# NOTE: Uncomment and run line 45, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

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
library(raster) #to use crop
library(starsExtra) #to use dist_to_nearest
library(geostatsp) 
library(tmaptools)
library(mgcv) #to use GAM function 
library(plotly) #to 3d plot variables
library(MuMIn) #to be able to use dredge
library(visreg) # to be able to plot Aspect/categorical variables with GAM
library(performance) #for converting concurvity values to VIF
library(spdep) #for checking the moran's I

# loading in the processed tree data 
# NOTE: Uncomment and run line 45, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")

#### Generalized Additive Models ####

# For each population/all populations and each size/shape metric we created
#generalized additive models whereby...

# a) We removed NAs from the explanatory and response variables to be able to run the GAMs
# b) We then created a base GAM model with smoothing splines on all of the quantitative variables (slope, distance, and elevation)
#and not with aspect because categorical variables cannot be smoothed in GAMs.
# c) We then looked at the summary of the base model (looking at significant varibales) and used the dredge function to create subsequent models 
#using different combinations of significant variables and smoothing splines to compare the models with each other to see which one performs the best. 
# d) We compared the models by using ANOVA f-tests (only works if one models variables are entirely present in the other model)
#and by comparing the models' AIC values to see which one had the lowest, supporting which model fits the data the best
# e) Once we felt like we had narrowed down our choices to the best option, we checked if the conditions of the GAM were met well (normal distribution and equal variance 
#of the residuals since we are using a GAM with a Gaussian distribution)
# h) We then plotted the GAMs in 2D and 3D to be able to describe the relationship between the explanatory variables and the size/shape characteristics 
# i) We then looked for interactions and compared the interaction model to the previous best model and if we have a significant interaction term we have to consider
# j) We then plotted the interactions both to be able to describe the relationships


#We used only the 8 categories to get a more specific idea of how direction may be influencing size 


### LM ###

# removing the NAs from the explanatory and response variables to avoid issues while making the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA <- LM_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.FIXED)) %>% #Elevation NAs removed
  filter(!is.na(LM_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(LM_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(LM_Eastness)) %>% #eastness NAs removed
  filter(!is.na(LM_Northness)) %>% #northness NAs removed
  filter(!is.na(LM_TWI_values)) %>% #TWI NAs removed
  filter(!is.na(heat.load)) %>% #HLI NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #crown spread NAs removed
  filter(!is.na(DBH_ag)) %>% #DBH NAs removed
  filter(!is.na(Canopy_short_lg)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long_lg)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area_lg)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread_lg)) %>% #crown spread NAs removed
  filter(!is.na(DBH_ag_lg)) %>% #DBH NAs removed
  filter(!is.na(Canopy_short_inv)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long_inv)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area_inv)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread_inv)) %>% #crown spread NAs removed
  filter(!is.na(DBH_ag_inv)) #DBH NAs removed

## SCA ##


#removing the spatial geometry to be able to use the GAM function
LM_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LM_fixed_field_data_processed_terrain_dist_no_NA)

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short_lg <- log1p(LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short) #removes the infinites allowing us to run the GAM

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable and that we added a control for spatial autocorrelation with the s(x.1, y.1)
LM_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_TWI_values) + s(heat.load) + s(X.1, Y.1),    #s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness)
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#northness, twi has significant non-linear fit and heat load is borderline

#checking for concurvity (equivalent to co-linearity but for models with non-linear functions)
concurvity(LM_add.gam_SCA.terrain_dist, full = T)

library(performance)
check_concurvity(LM_add.gam_SCA.terrain_dist)

#high concurvity with all terms currently

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(elevation), s(northness), s(slope), s(TWI) as allowing for the best model

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge) 
#elevation is significant and twi is borderline

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.2 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_TWI_values) + s(X.1, Y.1), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.2) 
#twi is borderline

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.3 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_TWI_values) 
                                            + s(heat.load), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.3) 
#elevation is significant and twi is borderline

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.4 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_TWI_values) 
                                            + s(heat.load)  + s(X.1, Y.1), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.4)
#twi is borderline

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.5 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_TWI_values) + 
                                            s(heat.load) + s(d), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation has significant non-linear fit and TWI are borderline

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.6 <- gam(Canopy_short ~ s(Elevation..m.FIXED), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.6) 

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.7 <- gam(Canopy_short ~ s(Elevation..m.FIXED)  + s(X.1, Y.1), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.7) 

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.8 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(heat.load) , 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.8) 
#elevation is significant

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge.9 <- gam(Canopy_short ~ s(Elevation..m.FIXED)  + s(heat.load)  + s(X.1, Y.1), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.9) 
#none are signigicant

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LM_add.gam_SCA.terrain_dist.dredge.9) #storing the residuals of the GAM
coords <- data.frame(
  x = LM_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LM_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 15)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LM_add.gam_CA.terrain_dist controls for spatial autocorrelation
#LM_add.gam_CA.terrain_dist.dredge controls for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.2 controls for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.3 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.4 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.5 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.6 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.7 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.8 control for spatial autocorrelation
#LM_add.gam_SCA.terrain_dist.dredge.9 control for spatial autocorrelation

# checking concurvity
check_concurvity(LM_add.gam_SCA.terrain_dist.dredge.6)

#LM_add.gam_SCA.terrain_dist high concurvity across all terms
#LM_add.gam_SCA.terrain_dist.dredge elevation and TWI have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.2 elevation, TWI, location have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.3 elevation and twi have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.4 elevation, twi, heat load, location have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.5 elevation, twi, distance have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.6 low concurvity
#LM_add.gam_SCA.terrain_dist.dredge.7 elevation, location have high concurvity
#LM_add.gam_SCA.terrain_dist.dredge.8 low concurvity
#LM_add.gam_SCA.terrain_dist.dredge.9 elevation, heat load, location have high concurvity

#checking normality of residuals
shapiro.test(LM_add.gam_SCA.terrain_dist.dredge.9$residuals)

#LM_add.gam_SCA.terrain_dist not normal
#LM_add.gam_SCA.terrain_dist.dredge not normal
#LM_add.gam_SCA.terrain_dist.dredge.2 not normal, but the least not normal
#LM_add.gam_SCA.terrain_dist.dredge.3 not normal
#LM_add.gam_SCA.terrain_dist.dredge.4 not normal
#LM_add.gam_SCA.terrain_dist.dredge.5 not normal
#LM_add.gam_SCA.terrain_dist.dredge.6 not normal
#LM_add.gam_SCA.terrain_dist.dredge.7 not normal
#LM_add.gam_SCA.terrain_dist.dredge.8 not normal
#LM_add.gam_SCA.terrain_dist.dredge.9 not normal

#model 6 and 8 are the best supported
AIC(LM_add.gam_SCA.terrain_dist.dredge.6, LM_add.gam_SCA.terrain_dist.dredge.8)
anova(LM_add.gam_SCA.terrain_dist.dredge.6, LM_add.gam_SCA.terrain_dist.dredge.8)

#I will go with LM_add.gam_SCA.terrain_dist.dredge.8 because it is a fuller model and not significantly worse than model 6


#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_SCA.terrain_dist.dredge.8 is the best supported model


#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_SCA.terrain_dist.dredge 
summary(LM_add.gam_SCA.terrain_dist.dredge.8)
check_concurvity(LM_add.gam_SCA.terrain_dist.dredge.8)
#elevation and northness are significant

#while the original dredge model has the lowest AIC the dredge 3 model has the lowest AIC without concurvity

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.terrain_dist.dredge.8) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LM_add.gam_SCA.terrain_dist.dredge.8$residuals)

#looking at significance
summary(LM_add.gam_SCA.terrain_dist.dredge.8)

#Chosen model: LM_add.gam_SCA.terrain_dist.dredge.3

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge.8, select=2, 
         all.terms=T, xlab = "Heat Load", 
         ylab = expression(f[1]*'Heat Load'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.terrain_dist.dredge.8, select=1, 
         all.terms=T, xlab = "Elevation (m)", ylab = expression(f[1]*'Elevation (m)'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_Northness, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, heat.load), 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):Northness'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_SCA.inter, LM_add.gam_SCA.terrain_dist.dredge.8)

#overall best model:LM_add.gam_SCA.terrain_dist.dredge.8


## LCA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short_lg <- log1p(LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short) #removes the infinites allowing us to run the GAM

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + 
                                     s(LM_Eastness) + s(LM_TWI_values) + s(heat.load) + s(X.1, Y.1),   
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#TWI and elevations showing significant smoothness

#checking for concurvity (equivalent to co-linearity but for models with non-linear functions)
check_concurvity(LM_add.gam_LCA.terrain_dist)

#high concurvity across all terms

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(elevation), s(Northness), s(slope), s(TWI) as allowing for the best model

# Checking a GAM after dredge
LM_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_Northness) + 
                                            s(LM_slope_raster_15_data_pts) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(elevation) and s(TWI) are significant, northness is borderline

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.dredge) # AICs
anova(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

# checking concurvity
check_concurvity(LM_add.gam_LCA.terrain_dist.dredge) #high concurvity between elevation and TWI

# Checking a GAM without TWI
LM_add.gam_LCA.terrain_dist.dredge.2 <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_Northness) + 
                                              s(LM_slope_raster_15_data_pts),
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(elevation) and s(northness) are significant, northness is borderline

# Checking a GAM without Elevation
LM_add.gam_LCA.terrain_dist.dredge.3 <- gam(Canopy_long ~ s(LM_Northness) + 
                                              s(LM_slope_raster_15_data_pts) + s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(northness) and s(TWI) are significant, northness is borderline

#checking AIC for these models
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.dredge, LM_add.gam_LCA.terrain_dist.dredge.2, LM_add.gam_LCA.terrain_dist.dredge.3)
#LM_add.gam_LCA.terrain_dist has the lowest AIC

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LM_add.gam_LCA.terrain_dist.dredge.3) #storing the residuals of the GAM
coords <- data.frame(
  x = LM_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LM_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 15)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LM_add.gam_LCA.terrain_dist controls for spatial autocorrelation
#LM_add.gam_LCA.terrain_dist.dredge controls for spatial autocorrelation
#LM_add.gam_LCA.terrain_dist.dredge.2 controls for spatial autocorrelation
#LM_add.gam_LCA.terrain_dist.dredge.3 controls for spatial autocorrelation

# checking concurvity
check_concurvity(LM_add.gam_LCA.terrain_dist.dredge.3)

#LM_add.gam_LCA.terrain_dist high concurvity across all terms
#LM_add.gam_LCA.terrain_dist.dredge elevation and TWI have high concurvity
#LM_add.gam_LCA.terrain_dist.dredge.2 low concurvity across all terms
#LM_add.gam_LCA.terrain_dist.dredge.3 low concurvity across all terms

#checking normality of residuals
shapiro.test(LM_add.gam_LCA.terrain_dist.dredge.3$residuals)

#LM_add.gam_LCA.terrain_dist not normal
#LM_add.gam_LCA.terrain_dist.dredge not normal
#LM_add.gam_LCA.terrain_dist.dredge.2 not normal, but the least not normal
#LM_add.gam_LCA.terrain_dist.dredge.3 not normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_LCA.terrain_dist.dredge.3 is the best supported model

#checking if the model without location meets spatial autocorrelation 

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_LCA.terrain_dist.dredge.3 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_LCA.terrain_dist.dredge 
summary(LM_add.gam_LCA.terrain_dist.dredge.3)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.terrain_dist.dredge.3) #pretty normal residuals and no heteroscedasticity 

#looking at significance
summary(LM_add.gam_LCA.terrain_dist.dredge.3)

#Chosen model: LM_add.gam_LCA.terrain_dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = 'Slope (º)', ylab = expression(f[1]*'Slope (º)'))
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = "LM_Northness", 
         ylab = expression(f[1]*'(LM_Northness)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.3, select=3, 
         all.terms=T, xlab = "TWI", ylab = expression(f[1]*'TWI'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_TWI_values, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_Northness, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_LCA.inter <- gam(Canopy_long ~ ti(LM_TWI_values, LM_Northness, LM_slope_raster_15_data_pts), 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):Northness'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_LCA.inter, LM_add.gam_LCA.terrain_dist)

#overall best model: LM_add.gam_LCA.terrain_dist

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_area_lg <- log1p(LM_fixed_field_data_processed_terrain_dist_no_NA$Canopy_area) #removes the infinites allowing us to run the GAM

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness) + s(LM_TWI_values) + s(X.1, Y.1), 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation, eastness, and TWI showing potential for significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(northness), s(slope), s(TWI)

# Checking a GAM without aspect
LM_add.gam_CA.terrain_dist.dredge <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + 
                                           s(LM_Northness) + s(LM_TWI_values), 
                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation), s(northness), and somewhat s(TWI) have significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.dredge) # AICs
anova(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.dredge)  #ANOVA F-Test
#the dredged model is preferred

#checking concurvity
check_concurvity(LM_add.gam_CA.terrain_dist.dredge) #low concurvity

#elevation and twi have high concurvity

# model without TWI
LM_add.gam_CA.terrain_dist.dredge.2 <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness), 
                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) and somewhat s(northness) have significant non-linear functions

# model without elevation
LM_add.gam_CA.terrain_dist.dredge.3 <- gam(Canopy_area_lg ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(northness) has significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.dredge, LM_add.gam_CA.terrain_dist.dredge.2, LM_add.gam_CA.terrain_dist.dredge.3) # AICs

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LM_add.gam_CA.terrain_dist.dredge.3) #storing the residuals of the GAM
coords <- data.frame(
  x = LM_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LM_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 15)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LM_add.gam_CA.terrain_dist controls for spatial autocorrelation
#LM_add.gam_CA.terrain_dist.dredge controls for spatial autocorrelation
#LM_add.gam_CA.terrain_dist.dredge.2 controls for spatial autocorrelation
#LM_add.gam_CA.terrain_dist.dredge.3 controls for spatial autocorrelation

# checking concurvity
check_concurvity(LM_add.gam_CA.terrain_dist.dredge.3)

#LM_add.gam_CA.terrain_dist high concurvity across all terms
#LM_add.gam_CA.terrain_dist.dredge elevation and TWI have high concurvity
#LM_add.gam_CA.terrain_dist.dredge.2 low concurvity across all terms
#LM_add.gam_CA.terrain_dist.dredge.3 low concurvity across all terms

#checking normality of residuals
shapiro.test(LM_add.gam_CA.terrain_dist.dredge.3$residuals)

#LM_add.gam_CA.terrain_dist not normal
#LM_add.gam_CA.terrain_dist.dredge not normal
#LM_add.gam_CA.terrain_dist.dredge.2 not normal, but the least not normal
#LM_add.gam_CA.terrain_dist.dredge.3 not normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_CA.terrain_dist.dredge.2 is the best supported model


#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CA.terrain_dist.dredge.2 
summary(LM_add.gam_CA.terrain_dist.dredge.2)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.terrain_dist.dredge.2) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_CA.terrain_dist.dredge.2)

#Chosen model: LM_add.gam_CA.terrain_dist.no.twi

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.terrain_dist.dredge.2, select=2, 
         all.terms=T, xlab = 'Slope (º)', ylab = expression(f[1]*'Slope (º)'))
plot.gam(LM_add.gam_CA.terrain_dist.dredge.2, select=3, 
         all.terms=T, xlab = "LM_Northness", 
         ylab = expression(f[1]*'(LM_Northness)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.terrain_dist.dredge.2, select=1, 
         all.terms=T, xlab = "Elevation (m)", ylab = expression(f[1]*'Elevation (m)'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_Northness, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m.FIXED, LM_Northness, LM_slope_raster_15_data_pts), 
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_CA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_CA.inter, LM_add.gam_CA.terrain_dist.no.aspect)

#overall best model: LM_add.gam_CA.terrain_dist.no.aspect

## CS ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness) + s(LM_TWI_values) + s(X.1, Y.1), 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#northness and TWI showing signs of significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM after the dredge
LM_add.gam_CS.terrain_dist.dredge <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + 
                                           s(LM_Northness) + s(LM_TWI_values), 
                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) is significant, signs that northness and twi are close to significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.dredge) # AICs
anova(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preffered

#checking concurvity
check_concurvity(LM_add.gam_CS.terrain_dist.dredge)
#elevation and TWI have high concurvity

# model without TWI
LM_add.gam_CS.terrain_dist.dredge.2 <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness), 
                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) and s(northness) have significant non-linear functions

# model without elevation
LM_add.gam_CS.terrain_dist.dredge.3 <- gam(Crown_spread ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(northness) and s(TWI) have significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.dredge, LM_add.gam_CS.terrain_dist.dredge.2, LM_add.gam_CS.terrain_dist.dredge.3) # AICs

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LM_add.gam_CS.terrain_dist.dredge.3) #storing the residuals of the GAM
coords <- data.frame(
  x = LM_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LM_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 15)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LM_add.gam_CS.terrain_dist controls for spatial autocorrelation
#LM_add.gam_CS.terrain_dist.dredge controls for spatial autocorrelation
#LM_add.gam_CS.terrain_dist.dredge.2 controls for spatial autocorrelation
#LM_add.gam_CS.terrain_dist.dredge.3 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(LM_add.gam_CS.terrain_dist.dredge.3)

#LM_add.gam_CS.terrain_dist high concurvity across all terms
#LM_add.gam_CS.terrain_dist.dredge elevation and TWI have high concurvity
#LM_add.gam_CS.terrain_dist.dredge.2 low concurvity across all terms
#LM_add.gam_CS.terrain_dist.dredge.3 low concurvity across all terms

#checking normality of residuals
shapiro.test(LM_add.gam_CS.terrain_dist.dredge.3$residuals)

#LM_add.gam_CS.terrain_dist not normal
#LM_add.gam_CS.terrain_dist.dredge not normal
#LM_add.gam_CS.terrain_dist.dredge.2 not normal, but the least not normal
#LM_add.gam_CS.terrain_dist.dredge.3 not normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_CS.terrain_dist.dredge.2 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CA.terrain_dist.dredge.2 
summary(LM_add.gam_CS.terrain_dist.dredge.2)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.terrain_dist.dredge.2) #pretty normal residuals and no Heteroscedasticity

#looking at significance
summary(LM_add.gam_CS.terrain_dist.dredge.2)

#Chosen model: LM_add.gam_CS.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CS.terrain_dist.no.twi, select=3, 
         all.terms=T, xlab = 'Northness', ylab = expression(f[1]*'Northness'))
plot.gam(LM_add.gam_CS.terrain_dist.no.twi, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CS.terrain_dist.no.twi, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_Northness, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, LM_Northness, LM_slope_raster_15_data_pts), 
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail, method = "REML")
summary(LM_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_CS.inter, LM_add.gam_CS.terrain_dist.no.aspect)

#overall best model: LM_add.gam_CS.terrain_dist.no.aspect

## DBH ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA$DBH_ag_lg <- log1p(LM_fixed_field_data_processed_terrain_dist_no_NA$DBH_ag) #removes the infinites allowing us to run the GAM

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_DBH.terrain_dist <- gam(DBH_ag_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness) + s(LM_TWI_values) + s(X.1, Y.1), 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only elevation is showing signs of significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(eastness), s(twi)

# Checking a GAM after the dredge
LM_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~ s(LM_Eastness) +
                                            s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness and twi showing somewhat significance

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.dredge) # AICs
anova(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.dredge)  #ANOVA F-Test
#the more complex, original model is preferred

# Checking a GAM after the dredge
LM_add.gam_DBH.terrain_dist.dredge.2 <- gam(DBH_ag_lg ~ s(LM_Eastness) + s(Elevation..m.FIXED) + s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#twi showing somewhat significance

# Checking a GAM after the dredge
LM_add.gam_DBH.terrain_dist.dredge.3 <- gam(DBH_ag_lg ~ s(LM_Eastness) + s(LM_slope_raster_15_data_pts) + 
                                              s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness showing somewhat significance, twi is borderline

#comparing AIC values
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.dredge, LM_add.gam_DBH.terrain_dist.dredge.2, LM_add.gam_DBH.terrain_dist.dredge.3)
#the second and third dredge models are preferred

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LM_add.gam_DBH.terrain_dist.dredge.3) #storing the residuals of the GAM
coords <- data.frame(
  x = LM_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LM_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 15)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LM_add.gam_DBH.terrain_dist controls for spatial autocorrelation
#LM_add.gam_DBH.terrain_dist.dredge controls for spatial autocorrelation
#LM_add.gam_DBH.terrain_dist.dredge.2 controls for spatial autocorrelation
#LM_add.gam_DBH.terrain_dist.dredge.3 controls for spatial autocorrelation

# checking concurvity
check_concurvity(LM_add.gam_DBH.terrain_dist)

#LM_add.gam_DBH.terrain_dist high concurvity across all terms
#LM_add.gam_DBH.terrain_dist.dredge eastness and TWI have high concurvity
#LM_add.gam_DBH.terrain_dist.dredge.2 elevation and TWI have high concurvity
#LM_add.gam_DBH.terrain_dist.dredge.3 low concurvity across all terms

#checking normality of residuals
shapiro.test(LM_add.gam_DBH.terrain_dist$residuals)

#LM_add.gam_DBH.terrain_dist not normal
#LM_add.gam_DBH.terrain_dist.dredge not normal
#LM_add.gam_DBH.terrain_dist.dredge.2 not normal
#LM_add.gam_DBH.terrain_dist.dredge.3 not normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LM_add.gam_DBH.terrain_dist.dredge.3 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LM_add.gam_DBH.terrain_dist.dredge.3)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.terrain_dist.dredge.3) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_DBH.terrain_dist.dredge.3)

#Chosen model: LM_add.gam_DBH.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = 'Eastness', ylab = expression(f[1]*'(Eastness)'))
plot.gam(LM_add.gam_DBH.terrain_dist.dredge.3, select=3, 
         all.terms=T, xlab = "TWI", 
         ylab = expression(f[1]*'(TWI)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_Eastness, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        type="scatter3d", mode="markers")


#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m.FIXED, LM_Eastness, LM_slope_raster_15_data_pts), 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):Elevation'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_DBH.inter, LM_add.gam_DBH.terrain_dist)
#while there is no significant interaction, the model with the interaction performs significantly better

#overall best model: LM_add.gam_DBH.terrain_dist.dredge.no.smooth (a multiple linear regression)

### LC ###

# removing the NAs from the explanatory and response variables to avoid issues while making the GAMs
LC_fixed_field_data_processed_terrain_dist_no_NA <- LC_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.FIXED)) %>% #Elevation NAs removed
  filter(!is.na(LC_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(LC_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(LC_Eastness)) %>% #eastness NAs removed
  filter(!is.na(LC_Northness)) %>% #northness NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed

## SCA ##

#removing the spatial geometry to be able to use the GAM function
LC_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LC_fixed_field_data_processed_terrain_dist_no_NA)

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs
LC_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short_lg <- log1p(LC_fixed_field_data_processed_terrain_dist_no_NA$Canopy_short) #removes the infinites allowing us to run the GAM

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_SCA.terrain_dist <- gam(Canopy_short_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation) and aspect

# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short_lg ~ s(d) + s(LC_Eastness) + s(X.1, Y.1), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has a significant non-linear fit

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.dredge) # AICs
anova(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants all three variables

#checking concurvity
check_concurvity(LC_add.gam_SCA.terrain_dist.dredge) #low concurvity

# Checking a GAM after dredge
LC_add.gam_SCA.terrain_dist.dredge.2 <- gam(Canopy_short_lg ~ s(d) + s(LC_Eastness), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness has a significant non-linear fit

# Checking a GAM after dredge
LC_add.gam_SCA.terrain_dist.dredge.3 <- gam(Canopy_short_lg ~ s(d) + s(LC_slope_raster_15_data_pts) + s(X.1,Y.1), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has a significant non-linear fit

# Checking a GAM after dredge
LC_add.gam_SCA.terrain_dist.dredge.4 <- gam(Canopy_short_lg ~ s(d) + s(LC_slope_raster_15_data_pts), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#distance and slope have significant non-linear fits

#checking the AICs
AIC(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.dredge, LC_add.gam_SCA.terrain_dist.dredge.2, 
    LC_add.gam_SCA.terrain_dist.dredge.3, LC_add.gam_SCA.terrain_dist.dredge.4) # AICs

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_SCA.terrain_dist) #storing the residuals of the GAM
coords <- data.frame(
  x = LC_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 10)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LC_add.gam_SCA.terrain_dist controls for spatial autocorrelation
#LC_add.gam_SCA.terrain_dist.dredge controls for spatial autocorrelation
#LC_add.gam_SCA.terrain_dist.dredge.2 does not control for spatial autocorrelation
#LC_add.gam_SCA.terrain_dist.dredge.3 controls for spatial autocorrelation
#LC_add.gam_SCA.terrain_dist.dredge.4 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(LC_add.gam_SCA.terrain_dist)

#LC_add.gam_SCA.terrain_dist high concurvity across all terms
#LC_add.gam_SCA.terrain_dist.dredge low concurvity across all terms
#LC_add.gam_SCA.terrain_dist.dredge.2 low concurvity across all terms
#LC_add.gam_SCA.terrain_dist.dredge.3 high concurvity with slope and location
#LC_add.gam_SCA.terrain_dist.dredge.4 low concurvity across all terms

#checking normality of residuals
shapiro.test(LC_add.gam_SCA.terrain_dist$residuals)

#LC_add.gam_SCA.terrain_dist normal
#LC_add.gam_SCA.terrain_dist.dredge normal
#LC_add.gam_SCA.terrain_dist.dredge.2 normal
#LC_add.gam_SCA.terrain_dist.dredge.3 normal
#LC_add.gam_SCA.terrain_dist.dredge.4 normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
# LC_add.gam_SCA.terrain_dist.dredge is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_SCA.terrain_dist.dredge 
summary(LC_add.gam_SCA.terrain_dist.dredge)
#just location is significant

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_SCA.terrain_dist.dredge$residuals)

#looking at significance
summary(LC_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LC_add.gam_SCA.terrain_dist.dredge

#plotting all of the variables for observation
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = 'Eastness', ylab = expression(f[1]*'Eastness'))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = "distance to river", 
         ylab = expression(f[1]*'(distance to river)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.terrain_dist.dredge, select=3, 
         all.terms=T, xlab = "Location", ylab = expression(f[1]*'Location'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_Eastness, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$X.1, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_SCA.inter <- gam(Canopy_short ~ ti(d, LC_Eastness, X.1, Y.1), 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail, method = "REML")
summary(LC_add.gam_SCA.inter)
#there was a significant interaction term 

#interaction plots
plot.gam(LC_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Eastness:X:Y'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_SCA.inter, LC_add.gam_SCA.terrain_dist.dredge.3)

#overall best model: LC_add.gam_SCA.terrain_dist.dredge.3

## LCA ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_LCA.terrain_dist <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has borderline significance 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation)

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long_lg ~ s(d) + 
                                            s(LC_Eastness) + s(X.1,Y.1), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant, distance is borderline

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge) # AICs
anova(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.dredge.2 <- gam(Canopy_long_lg ~ s(d) + 
                                              s(LC_Eastness), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness is significant, distance is borderling

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.dredge.3 <- gam(Canopy_long_lg ~ s(d) + 
                                              s(LC_Eastness) + s(Elevation..m.FIXED) + s(LC_Northness), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and eastness is significant, northness is borderline

# Checking a GAM with eastness
LC_add.gam_LCA.terrain_dist.dredge.4 <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_Northness), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation is significant, northness is borderline

# Checking a GAM without northness
LC_add.gam_LCA.terrain_dist.dredge.5 <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_Eastness), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and eastness are significant

#checking AIC
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge, LC_add.gam_LCA.terrain_dist.dredge.2,
    LC_add.gam_LCA.terrain_dist.dredge.3, LC_add.gam_LCA.terrain_dist.dredge.4, LC_add.gam_LCA.terrain_dist.dredge.5)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_LCA.terrain_dist) #storing the residuals of the GAM
coords <- data.frame(
  x = LC_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 10)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LC_add.gam_LCA.terrain_dist controls for spatial autocorrelation
#LC_add.gam_LCA.terrain_dist.dredge controls for spatial autocorrelation
#LC_add.gam_LCA.terrain_dist.dredge.2 does not control for spatial autocorrelation
#LC_add.gam_LCA.terrain_dist.dredge.3 controls for spatial autocorrelation
#LC_add.gam_LCA.terrain_dist.dredge.4 controls for spatial autocorrelation
#LC_add.gam_LCA.terrain_dist.dredge.5 controls for spatial autocorrelation

# checking concurvity
check_concurvity(LC_add.gam_LCA.terrain_dist.dredge.2)

#LC_add.gam_LCA.terrain_dist high concurvity across all terms
#LC_add.gam_LCA.terrain_dist.dredge low concurvity across all terms
#LC_add.gam_LCA.terrain_dist.dredge.2 low concurvity across all terms
#LC_add.gam_LCA.terrain_dist.dredge.3 high concurvity with eastness and northness
#LC_add.gam_LCA.terrain_dist.dredge.4 low concurvity across all terms
#LC_add.gam_LCA.terrain_dist.dredge.5 low concurvity across all terms

#checking normality of residuals
shapiro.test(LC_add.gam_LCA.terrain_dist$residuals)

#LC_add.gam_LCA.terrain_dist normal
#LC_add.gam_LCA.terrain_dist.dredge normal
#LC_add.gam_LCA.terrain_dist.dredge.2 normal
#LC_add.gam_LCA.terrain_dist.dredge.3 normal
#LC_add.gam_LCA.terrain_dist.dredge.4 normal
#LC_add.gam_LCA.terrain_dist.dredge.5 normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although LC_add.gam_LCA.terrain_dist.dredge has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# LC_add.gam_LCA.terrain_dist.dredge.5 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_LCA.terrain_dist.dredge 
summary(LC_add.gam_LCA.terrain_dist.dredge.5)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.terrain_dist.dredge.5) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_LCA.terrain_dist.dredge.5$residuals)

#looking at significance
summary(LC_add.gam_LCA.terrain_dist.dredge.5)

#Chosen model: LC_add.gam_LCA.terrain_dist.dredge.5

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.5, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.5, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.5, select=3, 
         all.terms=T, xlab = "Eastness", ylab = expression(f[1]*'Eastness'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_LCA.inter, LC_add.gam_LCA.terrain_dist.dredge.5)

#overall best model: LC_add.gam_LCA.terrain_dist.dredge.5

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(location), s(eastness)

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge <- gam(Canopy_area_lg ~ s(d) + s(X.1, Y.1) + s(LC_Eastness), 
                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.dredge) # AICs
anova(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.2 <- gam(Canopy_area_lg ~ s(X.1, Y.1) + s(d), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.2)
#location is significant, distance is borderline

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.2 <- gam(Canopy_area_lg ~ s(d), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.2)
#distance is significant

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.3 <- gam(Canopy_area_lg ~ s(d) + s(LC_Eastness) + s(X.1,Y.1), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.3)
#location is significant and distance is borderline

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge.4 <- gam(Canopy_area_lg ~ s(d) + s(LC_Eastness), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness is significant and distance is borderling

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge.5 <- gam(Canopy_area_lg ~ s(d) + s(LC_slope_raster_15_data_pts) + s(X.1,Y.1), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge.6 <- gam(Canopy_area_lg ~ s(d) + s(LC_slope_raster_15_data_pts), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#distance is significant and slope is borderline

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.dredge,
    LC_add.gam_CA.terrain_dist.dredge.2, LC_add.gam_CA.terrain_dist.dredge.3, 
    LC_add.gam_CA.terrain_dist.dredge.4, LC_add.gam_CA.terrain_dist.dredge.5,
    LC_add.gam_CA.terrain_dist.dredge.6) # AICs
#the dredge model with s(d), s(eastness), s(location) is preferred

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_CA.terrain_dist) #storing the residuals of the GAM
coords <- data.frame(
  x = LC_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 10)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LC_add.gam_CA.terrain_dist controls for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge controls for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge.2 does not control for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge.3 controls for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge.4 does not control for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge.5 controls for spatial autocorrelation
#LC_add.gam_CA.terrain_dist.dredge.6 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(LC_add.gam_CA.terrain_dist)

#LC_add.gam_CA.terrain_dist high concurvity across all terms
#LC_add.gam_CA.terrain_dist.dredge low concurvity across all terms
#LC_add.gam_CA.terrain_dist.dredge.2 low concurvity across all terms
#LC_add.gam_CA.terrain_dist.dredge.3 low concurvity across all terms
#LC_add.gam_CA.terrain_dist.dredge.4 low concurvity across all terms
#LC_add.gam_CA.terrain_dist.dredge.5 High concurvity with slope and location
#LC_add.gam_CA.terrain_dist.dredge.6 low concurvity across all terms

#checking normality of residuals
shapiro.test(LC_add.gam_CA.terrain_dist.dredge$residuals)

#LC_add.gam_CA.terrain_dist normal
#LC_add.gam_CA.terrain_dist.dredge normal
#LC_add.gam_CA.terrain_dist.dredge.2 normal
#LC_add.gam_CA.terrain_dist.dredge.3 normal
#LC_add.gam_CA.terrain_dist.dredge.4 normal
#LC_add.gam_CA.terrain_dist.dredge.5 normal
#LC_add.gam_CA.terrain_dist.dredge.6 normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although LC_add.gam_CA.terrain_dist.dredge has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# LC_add.gam_CA.terrain_dist.dredge is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CA.terrain_dist.dredge 
summary(LC_add.gam_CA.terrain_dist.dredge)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_CA.terrain_dist.dredge$residuals)

#looking at significance
summary(LC_add.gam_CA.terrain_dist.dredge)

#Chosen model: LC_add.gam_CA.terrain_dist.dredge

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CA.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = "Location", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.terrain_dist.dredge, select=3, 
         all.terms=T, xlab = "Eastness", ylab = expression(f[1]*'Eastness'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_add.gam_CA.terrain_dist.dredge$Elevation..m.FIXED, 
        y=LC_add.gam_CA.terrain_dist.dredge$d, 
        z=LC_add.gam_CA.terrain_dist.dredge$LC_slope_raster_15_data_pts, 
        color = LC_add.gam_CA.terrain_dist.dredge$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_CA.inter <- gam(Canopy_area ~ ti(d,X.1,Y.1, LC_Eastness), 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CA.inter)
#there was is a significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_CA.inter, LC_add.gam_CA.terrain_dist.dredge)

#overall best model: LC_add.gam_CA.terrain_dist.no.aspect

## CS ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_CS.terrain_dist <- gam(Crown_spread_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + 
                                    s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has significant non-linear functions, distance is borderline

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge <- gam(Crown_spread_lg ~ s(d) + s(X.1, Y.1), 
                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(distance) + s(location) is borderline significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.dredge) # AICs
anova(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.2 <- gam(Crown_spread_lg ~ s(d) + s(LC_Eastness) + s(X.1, Y.1), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant, distance is borderline

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.3 <- gam(Crown_spread_lg ~ s(d) + s(LC_Eastness), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(eastness) is  significant, distance is borderline

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.4 <- gam(Crown_spread_lg ~ s(d) + s(LC_slope_raster_15_data_pts) + s(X.1, Y.1), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
# location is significant, distance is borderline

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.5 <- gam(Crown_spread_lg ~ s(d) + s(LC_slope_raster_15_data_pts), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(distance) is significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.dredge,
    LC_add.gam_CS.terrain_dist.dredge.2, LC_add.gam_CS.terrain_dist.dredge.3, 
    LC_add.gam_CS.terrain_dist.dredge.4, LC_add.gam_CS.terrain_dist.dredge.5) # AICs
#the dredge model with s(d), s(eastness), s(location) is preferred

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_CS.terrain_dist) #storing the residuals of the GAM
coords <- data.frame(
  x = LC_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 10)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LC_add.gam_CS.terrain_dist controls for spatial autocorrelation
#LC_add.gam_CS.terrain_dist.dredge controls for spatial autocorrelation
#LC_add.gam_CS.terrain_dist.dredge.2 controls for spatial autocorrelation
#LC_add.gam_CS.terrain_dist.dredge.3 does not control for spatial autocorrelation
#LC_add.gam_CS.terrain_dist.dredge.4 controls for spatial autocorrelation
#LC_add.gam_CS.terrain_dist.dredge.5 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(LC_add.gam_CS.terrain_dist)

#LC_add.gam_CS.terrain_dist high concurvity across all terms
#LC_add.gam_CS.terrain_dist.dredge low concurvity across all terms
#LC_add.gam_CS.terrain_dist.dredge.2 low concurvity across all terms
#LC_add.gam_CS.terrain_dist.dredge.3 low concurvity across all terms
#LC_add.gam_CS.terrain_dist.dredge.4 High concurvity with slope and location
#LC_add.gam_CS.terrain_dist.dredge.5 low concurvity across all terms

#checking normality of residuals
shapiro.test(LC_add.gam_CS.terrain_dist$residuals)

#LC_add.gam_CS.terrain_dist normal
#LC_add.gam_CS.terrain_dist.dredge normal
#LC_add.gam_CS.terrain_dist.dredge.2 normal
#LC_add.gam_CS.terrain_dist.dredge.3 normal
#LC_add.gam_CS.terrain_dist.dredge.4 normal
#LC_add.gam_CS.terrain_dist.dredge.5 normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although LC_add.gam_CA.terrain_dist.dredge has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# LC_add.gam_CS.terrain_dist.dredge is the best supported model

#LC_add.gam_CS.terrain_dist.dredge is currently the best model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CS.terrain_dist.no.aspect 
summary(LC_add.gam_CS.terrain_dist.dredge)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_CS.terrain_dist.dredge)

#Chosen model: LC_add.gam_CS.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CS.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CS.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = "Location", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")



# 3d plotting in plotly and with gg3D
plot_ly(x=LC_add.gam_CS.terrain_dist.dredge$Y.1, 
        y=LC_add.gam_CS.terrain_dist.dredge$d, 
        z=LC_add.gam_CS.terrain_dist.dredge$X.1, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_CS.inter <- gam(Crown_spread ~ ti(X.1, d, Y.1), 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_CS.inter, LC_add.gam_CS.terrain_dist.dredge)

#overall best model: LC_add.gam_CS.terrain_dist.dredge

## DBH ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_DBH.terrain_dist <- gam(DBH_ag_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + 
                                     s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only location has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(location)

# Checking a GAM without aspect
LC_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~ s(d) + s(X.1,Y.1), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant 

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.dredge) # AICs
anova(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.2 <- gam(DBH_ag_lg ~ s(d),
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.2)
#distance is significant

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.3 <- gam(DBH_ag_lg ~ s(d) + s(X.1,Y.1) + 
                                              s(LC_Eastness) + s(LC_TWI_values), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.3)
#location is significant 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.4 <- gam(DBH_ag_lg ~ s(d) + s(LC_Eastness) + s(LC_TWI_values), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.4)
#s(eastiness) is significantly non-linear 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.5 <- gam(DBH_ag_lg ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) 
                                            + s(X.1,Y.1), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.5)
#s(eastiness) is significantly non-linear 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.6 <- gam(DBH_ag_lg ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.6)
#s(eastiness) is significantly non-linear 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.7 <- gam(DBH_ag_lg ~ s(LC_slope_raster_15_data_pts) 
                                            + s(X.1,Y.1), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.7)
#s(eastiness) is significantly non-linear 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.8 <- gam(DBH_ag_lg ~ s(LC_slope_raster_15_data_pts), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.8)
#s(eastiness) is significantly non-linear 

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.9 <- gam(DBH_ag_lg ~ s(Elevation..m.FIXED) + s(X.1,Y.1), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.9)
#s(eastiness) is significantly non-linear 


#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.10 <- gam(DBH_ag_lg ~ s(Elevation..m.FIXED), 
                                             data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.10)
#s(eastiness) is significantly non-linear 


#comparing the models 
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.dredge, LC_add.gam_DBH.terrain_dist.dredge.2, 
    LC_add.gam_DBH.terrain_dist.dredge.3, LC_add.gam_DBH.terrain_dist.dredge.4,
    LC_add.gam_DBH.terrain_dist.dredge.5, LC_add.gam_DBH.terrain_dist.dredge.6,
    LC_add.gam_DBH.terrain_dist.dredge.7, LC_add.gam_DBH.terrain_dist.dredge.8,
    LC_add.gam_DBH.terrain_dist.dredge.9, LC_add.gam_DBH.terrain_dist.dredge.10)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_DBH.terrain_dist.dredge.10) #storing the residuals of the GAM
coords <- data.frame(
  x = LC_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = LC_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 10)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#LC_add.gam_DBH.terrain_dist controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.2 does not control for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.3 controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.4 does not control for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.5 controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.6 does not control for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.7 controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.8 does not control for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.9 controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.10 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(LC_add.gam_DBH.terrain_dist.dredge.10)

#LC_add.gam_DBH.terrain_dist high concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge low concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge.2 low concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge.3 moderate concurvity with location, eastness, TWI
#LC_add.gam_DBH.terrain_dist.dredge.4 low concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge.5 high concurvity elevation, slope, and location
#LC_add.gam_DBH.terrain_dist.dredge.6 low concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge.7 moderate concurvity slope and location
#LC_add.gam_DBH.terrain_dist.dredge.8 low concurvity across all terms
#LC_add.gam_DBH.terrain_dist.dredge.9 high concurvity elevation and location
#LC_add.gam_DBH.terrain_dist.dredge.10 low concurvity across all terms

#checking normality of residuals
shapiro.test(LC_add.gam_DBH.terrain_dist.dredge.10$residuals)

#LC_add.gam_DBH.terrain_dist normal
#LC_add.gam_DBH.terrain_dist.dredge normal
#LC_add.gam_DBH.terrain_dist.dredge.2 normal
#LC_add.gam_DBH.terrain_dist.dredge.3 normal
#LC_add.gam_DBH.terrain_dist.dredge.4 normal
#LC_add.gam_DBH.terrain_dist.dredge.5 normal
#LC_add.gam_DBH.terrain_dist.dredge.6 not normal
#LC_add.gam_DBH.terrain_dist.dredge.7 normal
#LC_add.gam_DBH.terrain_dist.dredge.8 normal
#LC_add.gam_DBH.terrain_dist.dredge.9 normal
#LC_add.gam_DBH.terrain_dist.dredge.10 not normal

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although LC_add.gam_DBH.terrain_dist.dredge has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# LC_add.gam_DBH.terrain_dist.dredge is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_DBH.terrain_dist.dredge.3 
summary(LC_add.gam_DBH.terrain_dist.dredge)
#the model with just distance smooth is marginally better, so I will use the one with smoothing splines on all of the quantitative variables

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.terrain_dist.dredge.3) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_DBH.terrain_dist.dredge.3$residuals)

#looking at significance
summary(LC_add.gam_DBH.terrain_dist.dredge.3)

#Chosen model: LC_add.gam_DBH.terrain_dist.dredge.3

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_DBH.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_DBH.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = "Location", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.terrain_dist.dredge.3, select=3, 
         all.terms=T, xlab = "Eastness", 
         ylab = expression(f[1]*'(Eastness)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.terrain_dist.dredge.3, select=4, 
         all.terms=T, xlab = "TWI", 
         ylab = expression(f[1]*'(TWI)'), se = TRUE , col = "black")


# # 3d plotting in plotly and with gg3D
plot_ly(x=LC_add.gam_DBH.terrain_dist.dredge$d,
        y=LC_add.gam_DBH.terrain_dist.dredge$X.1,
        z=LC_add.gam_DBH.terrain_dist.dredge.3$Y.1,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_DBH.inter <- gam(DBH_ag ~ ti(X.1, Y.1, d), 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_DBH.inter, LC_add.gam_DBH.terrain_dist.dredge)
#while there is significant interaction, the model with the interaction performs significantly better

#overall best model: LC_add.gam_DBH.inter

### SD ###

# removing the NAs from the explanatory and response variables to avoid issues while making the GAMs
SD_fixed_field_data_processed_terrain_dist_no_NA <- SD_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.)) %>% #Elevation NAs removed
  filter(!is.na(SD_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(SD_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(SD_Eastness)) %>% #eastness NAs removed
  filter(!is.na(SD_Northness)) %>% #northness NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) %>% #DBH NAs removed
  mutate(Elevation..m. = as.numeric(Elevation..m.))

## SCA ##

#removing the spatial geometry to be able to use the GAM function
SD_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(SD_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_SCA.terrain_dist <- gam(Canopy_short_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + s(SD_Northness) + 
                                     s(SD_Eastness) + s(SD_TWI_values) + s(X.1, Y.1), 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location have significant non-linear function

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short_lg ~ s(Elevation..m.) + s(SD_Northness) +
                                            s(SD_slope_raster_15_data_pts), 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation, northness, and slope are significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.dredge) # AICs
anova(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.dredge)  #ANOVA F-Test
#the original model is preferred over the dredged one

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.2 <- gam(Canopy_short_lg ~ s(d) + s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has significant non-linear fit

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.3 <- gam(Canopy_short_lg ~ s(d) + s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and slope of the linear fits are significant

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.4 <- gam(Canopy_short_lg ~ s(Elevation..m.) + s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has significant non-linear fit, borderline with elevation

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.5 <- gam(Canopy_short_lg ~ s(Elevation..m.) + s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation has significant non-linear fit, borderline with TWI

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.6 <- gam(Canopy_short_lg ~ s(d) + s(SD_slope_raster_15_data_pts) + s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has significant non-linear fit

# Checking a GAM after dredge
SD_add.gam_SCA.terrain_dist.dredge.7 <- gam(Canopy_short_lg ~ s(d) + s(SD_slope_raster_15_data_pts) + s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.dredge.7) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#slope has significant non-linear fit

#comparing the models 
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.dredge,  
    SD_add.gam_SCA.terrain_dist.dredge.2, SD_add.gam_SCA.terrain_dist.dredge.3,
    SD_add.gam_SCA.terrain_dist.dredge.4, SD_add.gam_SCA.terrain_dist.dredge.5,
    SD_add.gam_SCA.terrain_dist.dredge.6, SD_add.gam_SCA.terrain_dist.dredge.7)

#SD_add.gam_SCA.terrain_dist.dredge.4 has the lowest AIC

# checking to see if we need to control for spatial autocorrelation

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(SD_add.gam_SCA.terrain_dist.dredge.7) #storing the residuals of the GAM
coords <- data.frame(
  x = SD_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 16) #K chosen based on what removes neighbor sub-graphs
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#SD_add.gam_SCA.terrain_dist controls for  spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge controls for  spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.2 controls for spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.3 does not control for spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.4 controls for spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.5 does not control for spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.6 controls for spatial autocorrelation
#SD_add.gam_SCA.terrain_dist.dredge.7 does not control for spatial autocorrelation

shapiro.test(SD_add.gam_SCA.terrain_dist.dredge.7$residuals) #.dredge

##SD_add.gam_SCA.terrain_dist meets normality of residuals
#SD_add.gam_DBH.terrain_dist.dredge does not meet normality of residuals
#SD_add.gam_SCA.terrain_dist.dredge.2 meets normality of residuals
#SD_add.gam_SCA.terrain_dist.dredge.3 does not meet normality of residuals
#SD_add.gam_DBH.terrain_dist.dredge.4 does not meet normality of residuals
#SD_add.gam_SCA.terrain_dist.dredge.5 meets normality of residuals
#SD_add.gam_SCA.terrain_dist.dredge.6 does not meet normality of residuals
#SD_add.gam_SCA.terrain_dist.dredge.7 does not meet normality of residuals

#only SD_add.gam_SCA.terrain_dist and SD_add.gam_SCA.terrain_dist.dredge.2 met normal residuals and control for spatial autocorrelation

#ANOVA to compare models with normal residuals and that control for spatial autocorrelation
anova(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.dredge.2)
#dredge model is fine

#Based on the comparisons (AIC/Anova) of these models and which ones control for spatial autoccorlation, the best model seems to be: SD_add.gam_SCA.terrain_dist.dredge.2 
summary(SD_add.gam_SCA.terrain_dist.dredge.2)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.terrain_dist.dredge.2) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(SD_add.gam_SCA.terrain_dist.dredge.5$residuals)

#looking at significance
summary(SD_add.gam_SCA.terrain_dist.dredge.2)

#Chosen model: SD_add.gam_SCA.terrain_dist.dredge.2

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_SCA.terrain_dist.dredge.2, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_SCA.terrain_dist.dredge.2, select=2, 
         all.terms=T, xlab = "TWI", 
         ylab = expression(f[1]*'(TWI)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_SCA.terrain_dist.dredge.2, select=3, 
         all.terms=T, xlab = "Location", ylab = expression(f[1]*'Location'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_TWI_values, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$X.1, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_SCA.inter <- gam(Canopy_short ~ ti(SD_TWI_values, d, X.1, Y.1), 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_SCA.inter, SD_add.gam_SCA.terrain_dist.no.aspect)

#overall best model:SD_add.gam_SCA.terrain_dist.no.aspect

## LCA ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + s(SD_Northness) + 
                                     s(SD_Eastness) + s(SD_TWI_values) + s(X.1, Y.1),  
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location is significant, and TWI is borderline

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(slope) and aspect

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long ~  s(SD_TWI_values) + s(X.1,Y.1, k = 13), 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long ~  s(SD_TWI_values), 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.2 <- gam(Canopy_long ~  s(Elevation..m.) + s(SD_TWI_values) + s(X.1,Y.1, k = 13), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.3 <- gam(Canopy_long ~  s(Elevation..m.) + s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.4 <- gam(Canopy_long ~  s(SD_slope_raster_15_data_pts) + s(SD_TWI_values) + s(X.1,Y.1, k = 13), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.5 <- gam(Canopy_long ~  s(SD_slope_raster_15_data_pts) + s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.6 <- gam(Canopy_long ~  s(SD_Eastness) + s(SD_slope_raster_15_data_pts) + 
                                              s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.7 <- gam(Canopy_long ~  s(SD_Eastness) + s(SD_slope_raster_15_data_pts) + 
                                              s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.7) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.8 <- gam(Canopy_long ~  s(SD_Eastness) + s(SD_slope_raster_15_data_pts) 
                                            + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.8) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.dredge.9 <- gam(Canopy_long ~  s(SD_Eastness) + s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.9) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant

#checking AIC
AIC(SD_add.gam_LCA.terrain_dist, SD_add.gam_LCA.terrain_dist.dredge, 
    SD_add.gam_LCA.terrain_dist.dredge.2, SD_add.gam_LCA.terrain_dist.dredge.3,
    SD_add.gam_LCA.terrain_dist.dredge.4, SD_add.gam_LCA.terrain_dist.dredge.5,
    SD_add.gam_LCA.terrain_dist.dredge.6, SD_add.gam_LCA.terrain_dist.dredge.7,
    SD_add.gam_LCA.terrain_dist.dredge.8, SD_add.gam_LCA.terrain_dist.dredge.9)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(SD_add.gam_LCA.terrain_dist.dredge.7) #storing the residuals of the GAM
coords <- data.frame(
  x = SD_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 17)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#SD_add.gam_LCA.terrain_dist controls for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge does not control for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.2 controls for spatial autocorrelation, barely
#SD_add.gam_LCA.terrain_dist.dredge.3 does not control for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.4 does not control for spatial autocorrelation, barely
#SD_add.gam_LCA.terrain_dist.dredge.5 does not control for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.6 controls for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.7 does not control for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.8 controls for spatial autocorrelation
#SD_add.gam_LCA.terrain_dist.dredge.9 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(SD_add.gam_LCA.terrain_dist.dredge.9)

#SD_add.gam_LCA.terrain_dist high concurvity across all terms
#SD_add.gam_LCA.terrain_dist.dredge low concurvity across all terms
#SD_add.gam_LCA.terrain_dist.dredge.2 high concurvity with elevation and location
#SD_add.gam_LCA.terrain_dist.dredge.3 low concurvity across all terms
#SD_add.gam_LCA.terrain_dist.dredge.4 moderate concurvity with slope and location
#SD_add.gam_LCA.terrain_dist.dredge.5 low concurvity across all terms
#SD_add.gam_LCA.terrain_dist.dredge.6 high concurvity on eastness, slope, and location
#SD_add.gam_LCA.terrain_dist.dredge.7 low concurvity across all terms
#SD_add.gam_LCA.terrain_dist.dredge.8 high concurvity eastness, slope, location
#SD_add.gam_LCA.terrain_dist.dredge.9 low concurvity across all terms

#checking normality of residuals
shapiro.test(SD_add.gam_LCA.terrain_dist.dredge.9$residuals)

#SD_add.gam_LCA.terrain_dist not normal
#SD_add.gam_LCA.terrain_dist.dredge not normal
#SD_add.gam_LCA.terrain_dist.dredge.2 not normal
#SD_add.gam_LCA.terrain_dist.dredge.3 not normal
#SD_add.gam_LCA.terrain_dist.dredge.4 not normal, closest to normal
#SD_add.gam_LCA.terrain_dist.dredge.5 not normal
#SD_add.gam_LCA.terrain_dist.dredge.6 not normal
#SD_add.gam_LCA.terrain_dist.dredge.7 not normal
#SD_add.gam_LCA.terrain_dist.dredge.8 not normal
#SD_add.gam_LCA.terrain_dist.dredge.9 not normal

#no good model, SD_add.gam_LCA.terrain_dist.dredge.4 is best of the bad models

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although SD_add.gam_LCA.terrain_dist.dredge.3 has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# SD_add.gam_LCA.terrain_dist.dredge.3 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist 
summary(SD_add.gam_LCA.terrain_dist.dredge.4)
#location is significant

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.terrain_dist.dredge.4) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(SD_add.gam_LCA.terrain_dist.dredge.4$residuals)

#looking at model significance
summary(SD_add.gam_LCA.terrain_dist.dredge.4)

#Chosen model: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_LCA.terrain_dist.dredge.4, select=2, 
         all.terms=T, xlab = 'TWI', ylab = expression(f[1]*'(TWI)'))
plot.gam(SD_add.gam_LCA.terrain_dist.dredge.4, select=3, 
         all.terms=T, xlab = "Location", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_LCA.terrain_dist.dredge.4, select=1, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$X.1, 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_LCA.inter <- gam(Canopy_long ~ ti(X.1, Y.1, d, SD_slope_raster_15_data_pts), 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_LCA.inter, SD_add.gam_LCA.terrain_dist.dredge.4)

#overall best model: SD_add.gam_LCA.inter

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + s(SD_Northness) + 
                                    s(SD_Eastness) + s(SD_TWI_values) + s(X.1, Y.1), 
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location is significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge <- gam(Canopy_area_lg ~ s(SD_TWI_values) + s(X.1, Y.1), 
                                         data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.dredge) # AICs
anova(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.dredge)  #ANOVA F-Test
#the full model is better

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.2 <- gam(Canopy_area_lg ~ s(SD_TWI_values), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#no signficance

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.3 <- gam(Canopy_area_lg ~ s(SD_slope_raster_15_data_pts) + s(SD_TWI_values) + s(X.1,Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.4 <- gam(Canopy_area_lg ~ s(SD_slope_raster_15_data_pts) + s(SD_TWI_values), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#slope is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.5 <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(SD_TWI_values) + s(X.1,Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant, twi is borderline

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.6 <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(SD_TWI_values), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation is sigificant 

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.7 <- gam(Canopy_area_lg ~ s(SD_slope_raster_15_data_pts) + s(X.1,Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.7) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant and slope is borderline

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.8 <- gam(Canopy_area_lg ~ s(SD_slope_raster_15_data_pts), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.8) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#slope is significant 

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.9 <- gam(Canopy_area_lg ~ s(SD_TWI_values) + s(X.1,Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.9) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.10 <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_Northness) +
                                              s(SD_slope_raster_15_data_pts) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.10) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.11 <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_Northness) +
                                              s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.11) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.12 <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_Eastness) +
                                              s(SD_slope_raster_15_data_pts) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.12) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

# Checking a GAM post dredge
SD_add.gam_CA.terrain_dist.dredge.13 <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_Eastness) +
                                              s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.13) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location is significant

#checking AIC
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.dredge,
    SD_add.gam_CA.terrain_dist.dredge.2, SD_add.gam_CA.terrain_dist.dredge.3,
    SD_add.gam_CA.terrain_dist.dredge.4, SD_add.gam_CA.terrain_dist.dredge.5,
    SD_add.gam_CA.terrain_dist.dredge.6, SD_add.gam_CA.terrain_dist.dredge.7, 
    SD_add.gam_CA.terrain_dist.dredge.8, SD_add.gam_CA.terrain_dist.dredge.9,
    SD_add.gam_CA.terrain_dist.dredge.10, SD_add.gam_CA.terrain_dist.dredge.11,
    SD_add.gam_CA.terrain_dist.dredge.12, SD_add.gam_CA.terrain_dist.dredge.13)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(SD_add.gam_CA.terrain_dist.dredge.13) #storing the residuals of the GAM
coords <- data.frame(
  x = SD_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 17)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#SD_add.gam_CA.terrain_dist controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.2 does not control for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.3 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.4 does not control for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.5 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.6 does not control for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.7 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.8 does not control for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.9 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.10 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.11 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.12 controls for spatial autocorrelation
#SD_add.gam_CA.terrain_dist.dredge.13 controls for spatial autocorrelation

# checking concurvity
check_concurvity(SD_add.gam_CA.terrain_dist.dredge)

#SD_add.gam_CA.terrain_dist high concurvity across all terms
#SD_add.gam_CA.terrain_dist.dredge moderate concurvity with twi and locatin
#SD_add.gam_CA.terrain_dist.dredge.2 low concurvity across all terms
#SD_add.gam_CA.terrain_dist.dredge.3 high concurvity with slope and location
#SD_add.gam_CA.terrain_dist.dredge.4 low concurvity across all terms
#SD_add.gam_CA.terrain_dist.dredge.5 high concurvity with elevation and location
#SD_add.gam_CA.terrain_dist.dredge.6 low concurvity across all terms
#SD_add.gam_CA.terrain_dist.dredge.7 high concurvity wih slope and location
#SD_add.gam_CA.terrain_dist.dredge.8 low concurvity across all terms
#SD_add.gam_CA.terrain_dist.dredge.9 moderate concurvity with twi and location
#SD_add.gam_CA.terrain_dist.dredge.10 high concurvity d, elevation, northness, slope, location
#SD_add.gam_CA.terrain_dist.dredge.11 moderate concurvity with elevation, northness, slope
#SD_add.gam_CA.terrain_dist.dredge.12 high concurvity with distance, elevation, eastness, slope, location
#SD_add.gam_CA.terrain_dist.dredge.13 moderate concurvity with eastness and slope

#checking normality of residuals
shapiro.test(SD_add.gam_CA.terrain_dist.dredge.13$residuals)

#SD_add.gam_CA.terrain_dist not normal
#SD_add.gam_CA.terrain_dist.dredge not normal
#SD_add.gam_CA.terrain_dist.dredge.2 not normal
#SD_add.gam_CA.terrain_dist.dredge.3 not normal
#SD_add.gam_CA.terrain_dist.dredge.4 not normal, closest to normal
#SD_add.gam_CA.terrain_dist.dredge.5 not normal
#SD_add.gam_CA.terrain_dist.dredge.6 normal
#SD_add.gam_CA.terrain_dist.dredge.7 not normal
#SD_add.gam_CA.terrain_dist.dredge.8 not normal, but barely
#SD_add.gam_CA.terrain_dist.dredge.9 not normal
#SD_add.gam_CA.terrain_dist.dredge.10 not normal
#SD_add.gam_CA.terrain_dist.dredge.11 not normal
#SD_add.gam_CA.terrain_dist.dredge.12 not normal
#SD_add.gam_CA.terrain_dist.dredge.13 not normal

#no good model, SD_add.gam_CA.terrain_dist.dredge is best of the bad models

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although SD_add.gam_CA.terrain_dist.dredge has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# SD_add.gam_CA.terrain_dist.dredge is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CA.terrain_dist.no.aspect 
summary(SD_add.gam_CA.terrain_dist.dredge)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_CA.terrain_dist.dredge)

#Chosen model: SD_add.gam_CA.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'TWI', ylab = expression(f[1]*'(TWI)'))
plot.gam(SD_add.gam_CA.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = "Location (m)", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$X.1, 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_TWI_values, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_CA.inter <- gam(Canopy_area ~ ti(X.1, Y.1, SD_TWI_values), 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_CA.inter, SD_add.gam_CA.terrain_dist.dredge)

#overall best model: SD_add.gam_CA.terrain_dist.dredge

## CS ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + s(SD_Northness) + 
                                    s(SD_Eastness) + s(SD_TWI_values) + s(X.1, Y.1),  
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location is significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope), and aspect

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge <- gam(Crown_spread ~ s(SD_slope_raster_15_data_pts) + s(X.1, Y.1), 
                                         data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only location is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.2 <- gam(Crown_spread ~ s(SD_slope_raster_15_data_pts), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only slope is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.3 <- gam(Crown_spread ~ s(SD_TWI_values) + s(X.1, Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only location is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.4 <- gam(Crown_spread ~ s(SD_TWI_values), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#no significance

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.5 <- gam(Crown_spread ~ s(SD_slope_raster_15_data_pts) + s(SD_TWI_values)  + s(X.1, Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only location is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.6 <- gam(Crown_spread ~ s(SD_slope_raster_15_data_pts) + s(SD_TWI_values), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only slope is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.7 <- gam(Crown_spread ~ s(d) + s(SD_slope_raster_15_data_pts) + s(X.1, Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.7) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only location is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.8 <- gam(Crown_spread ~ s(d) + s(SD_slope_raster_15_data_pts), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.8) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only slope is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.9 <- gam(Crown_spread ~ s(SD_Eastness) + s(SD_slope_raster_15_data_pts) + s(X.1, Y.1), 
                                           data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.9) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only location is significant

# Checking a GAM after dredge
SD_add.gam_CS.terrain_dist.dredge.10 <- gam(Crown_spread ~ s(SD_Eastness) + s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.dredge.10) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only slope is significant

#checking AIC
AIC(SD_add.gam_CS.terrain_dist, SD_add.gam_CS.terrain_dist.dredge,
    SD_add.gam_CS.terrain_dist.dredge.2, SD_add.gam_CS.terrain_dist.dredge.3,
    SD_add.gam_CS.terrain_dist.dredge.4, SD_add.gam_CS.terrain_dist.dredge.5, 
    SD_add.gam_CS.terrain_dist.dredge.6, SD_add.gam_CS.terrain_dist.dredge.7,
    SD_add.gam_CS.terrain_dist.dredge.8, SD_add.gam_CS.terrain_dist.dredge.9,
    SD_add.gam_CS.terrain_dist.dredge.10)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(SD_add.gam_CS.terrain_dist.dredge.3) #storing the residuals of the GAM
coords <- data.frame(
  x = SD_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 17)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#SD_add.gam_CS.terrain_dist controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.2 does not control for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.3 controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.4 does not control for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.5 controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.6 does not control for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.7 controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.8 does not control for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.9 controls for spatial autocorrelation
#SD_add.gam_CS.terrain_dist.dredge.10 does not control for spatial autocorrelation

# checking concurvity
check_concurvity(SD_add.gam_CS.terrain_dist.dredge.3)

#SD_add.gam_CS.terrain_dist high concurvity across all terms
#SD_add.gam_CS.terrain_dist.dredge high concurvity with slope and location
#SD_add.gam_CS.terrain_dist.dredge.2 low concurvity across all terms
#SD_add.gam_CS.terrain_dist.dredge.3 moderate concurvity with twi and location
#SD_add.gam_CS.terrain_dist.dredge.4 low concurvity across all terms
#SD_add.gam_CS.terrain_dist.dredge.5 high concurvity with slope and location
#SD_add.gam_CS.terrain_dist.dredge.6 low concurvity across all terms
#SD_add.gam_CS.terrain_dist.dredge.7 high concurvity wih slope, distance, and location
#SD_add.gam_CS.terrain_dist.dredge.8 low concurvity across all terms
#SD_add.gam_CS.terrain_dist.dredge.9 high concurvity with eastness, slope, and location
#SD_add.gam_CS.terrain_dist.dredge.10 low concurvity across all terms

#checking normality of residuals
shapiro.test(SD_add.gam_CS.terrain_dist.dredge.3$residuals)

#SD_add.gam_CS.terrain_dist not normal
#SD_add.gam_CS.terrain_dist.dredge not normal
#SD_add.gam_CS.terrain_dist.dredge.2 not normal
#SD_add.gam_CS.terrain_dist.dredge.3 not normal
#SD_add.gam_CS.terrain_dist.dredge.4 not normal
#SD_add.gam_CS.terrain_dist.dredge.5 not normal
#SD_add.gam_CS.terrain_dist.dredge.6 not normal
#SD_add.gam_CS.terrain_dist.dredge.7 normal
#SD_add.gam_CS.terrain_dist.dredge.8 not normal
#SD_add.gam_CS.terrain_dist.dredge.9 not normal
#SD_add.gam_CS.terrain_dist.dredge.10 not normal

#no good model, SD_add.gam_CS.terrain_dist.dredge.3 is best of the bad models

#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although SD_add.gam_CS.terrain_dist.dredge.3 has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# SD_add.gam_CS.terrain_dist.dredge.3 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CS.terrain_dist.no.aspect 
summary(SD_add.gam_CS.terrain_dist.dredge.3)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.terrain_dist.dredge.3) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_CS.terrain_dist.dredge.3)

#Chosen model: SD_add.gam_CS.terrain_dist.dredge.3

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CS.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = 'TWI', ylab = expression(f[1]*'(TWI)'))
plot.gam(SD_add.gam_CS.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = "Location (m)", 
         ylab = expression(f[1]*'(Location)'), se = TRUE , col = "black")


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_TWI_values, 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$X.1, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_CS.inter <- gam(Crown_spread ~ ti(SD_TWI_values, X.1, Y.1), 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_CS.inter, SD_add.gam_CS.terrain_dist.dredge.3)

#overall best model: SD_add.gam_CS.inter


## DBH ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_DBH.terrain_dist <- gam(DBH_ag_lg ~ s(d) + s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + s(SD_Northness) + 
                                     s(SD_Eastness) + s(SD_TWI_values) + s(X.1, Y.1),  
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#twi, northness, location have significant non-linear function , eastness is borderline

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_Eastness) + 
                                            s(SD_TWI_values) + s(X.1,Y.1), 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation, eastness, location is significant

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge.2 <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_Eastness) + 
                                              s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation is significant

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge.3 <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_Northness) + 
                                              s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
# twi, location is significant, elevation is borderline

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge.4 <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_Northness) + 
                                              s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#northness and elevation are significant 

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge.5 <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + 
                                              s(SD_TWI_values) + s(X.1,Y.1), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.5) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and location is significant, twi is borderline

# Checking a GAM post dredge
SD_add.gam_DBH.terrain_dist.dredge.6 <- gam(DBH_ag_lg ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + 
                                              s(SD_TWI_values), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail, method = "REML") #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.6) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation, slope is significant

#checking AICs
AIC(SD_add.gam_DBH.terrain_dist, SD_add.gam_DBH.terrain_dist.dredge,
    SD_add.gam_DBH.terrain_dist.dredge.2, SD_add.gam_DBH.terrain_dist.dredge.3,
    SD_add.gam_DBH.terrain_dist.dredge.4, SD_add.gam_DBH.terrain_dist.dredge.5,
    SD_add.gam_DBH.terrain_dist.dredge.6)

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(SD_add.gam_DBH.terrain_dist.dredge.6) #storing the residuals of the GAM
coords <- data.frame(
  x = SD_fixed_field_data_processed_terrain_dist_no_NA$X.1,
  y = SD_fixed_field_data_processed_terrain_dist_no_NA$Y.1
)
coords_mat <- as.matrix(coords)
#k-nearest neighbors 
knn <- knearneigh(coords_mat, k = 17)
nb <- knn2nb(knn)
#spatial weights
lw <- nb2listw(nb, style = "W")
#running the moran's I test
moran.test(res, lw)

#SD_add.gam_DBH.terrain_dist controls for spatial autocorrelation
#SD_add.gam_DBH.terrain_dist.dredge controls for spatial autocorrelation
#SD_add.gam_DBH.terrain_dist.dredge.2 controls for spatial autocorrelation, borderline
#SD_add.gam_DBH.terrain_dist.dredge.3 controls for spatial autocorrelation
#SD_add.gam_DBH.terrain_dist.dredge.4 controls for spatial autocorrelation
#SD_add.gam_DBH.terrain_dist.dredge.5 controls for spatial autocorrelation
#SD_add.gam_DBH.terrain_dist.dredge.6 controls for spatial autocorrelation

# checking concurvity
check_concurvity(SD_add.gam_DBH.terrain_dist.dredge.6)

#SD_add.gam_DBH.terrain_dist high concurvity across all terms
#SD_add.gam_DBH.terrain_dist.dredge high concurvity with elevation, eastness and location
#SD_add.gam_DBH.terrain_dist.dredge.2 low concurvity across all terms
#SD_add.gam_DBH.terrain_dist.dredge.3 high concurvity with elevation, northness, and location
#SD_add.gam_DBH.terrain_dist.dredge.4 low concurvity across all terms
#SD_add.gam_DBH.terrain_dist.dredge.5 high concurvity with elevation, slope and location
#SD_add.gam_DBH.terrain_dist.dredge.6 low concurvity across all terms

#checking normality of residuals
shapiro.test(SD_add.gam_DBH.terrain_dist.dredge.4$residuals)

#SD_add.gam_DBH.terrain_dist normal
#SD_add.gam_CS.terrain_dist.dredge not normal
#SD_add.gam_DBH.terrain_dist.dredge.2 normal
#SD_add.gam_DBH.terrain_dist.dredge.3 normal
#SD_add.gam_DBH.terrain_dist.dredge.4 normal
#SD_add.gam_DBH.terrain_dist.dredge.5 normal
#SD_add.gam_DBH.terrain_dist.dredge.6 normal

#SD_add.gam_DBH.terrain_dist.dredge.6 is best of the bad models


#based on the models with lower AIC, normality of residuals, low concurvity, amd control for spatial autocorrelation
#although SD_add.gam_DBH.terrain_dist.dredge.6 has the lowest AIC, if we do not need x.1, y.1 to control for spatial autocorrelation we should not include it
# SD_add.gam_DBH.terrain_dist.dredge.6 is the best supported model

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(SD_add.gam_DBH.terrain_dist.dredge.6)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.terrain_dist.dredge.6) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_DBH.terrain_dist.dredge.6)

#Chosen model: SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.6, select=3, 
         all.terms=T, xlab = 'TWI (m)', ylab = expression(f[1]*'(TWI)'))
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.6, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.6, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_TWI_values, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m., SD_TWI_values, SD_slope_raster_15_data_pts), 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_DBH.inter, SD_add.gam_DBH.terrain_dist.dredge.6)
#the interaction is not significant

#overall best model: SD_add.gam_DBH.terrain_dist.dredge.6 




