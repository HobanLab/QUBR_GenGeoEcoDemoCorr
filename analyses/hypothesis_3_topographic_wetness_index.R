# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Investigating to see if Quercus brandegeei tree shapes are influenced by Water Availability Proxies%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by their access to water using generalized additive models (GAMs)

# We used Generalized Additive Models (GAMs) after having trouble with multiple linear 
# regressions because of issues with the normality conditions and some nervousness about linearity. 
# We used four water availability proxies (distance to river, aspect, slope, and elevation).

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
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #crown spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed

## SCA ##

hist(LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts)

#removing the spatial geometry to be able to use the GAM function
LM_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LM_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable and that we added a control for spatial autocorrelation with the s(x.1, y.1)
LM_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness) + s(LM_TWI_values) + s(X.1, Y.1),   
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#checking for concurvity (equivalent to co-linearity but for models with non-linear functions)
concurvity(LM_add.gam_SCA.terrain_dist, full = T)

library(performance)
check_concurvity(LM_add.gam_SCA.terrain_dist)

#elevation has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(elevation), s(northness), s(slope), s(TWI) as allowing for the best model

# Checking a GAM with dredged variables
LM_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# it wants all four variables

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.dredge) # AICs
anova(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.dredge)  #ANOVA F-Test
#the less complex model is prefrered,  not significant ANOVA

check_concurvity(LM_add.gam_SCA.terrain_dist.dredge)

# Checking a GAM with significant smooth terms maintained
LM_add.gam_SCA.terrain_dist.dredge.2 <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_Northness + LM_TWI_values, 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist.dredge.2) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# it wants only s(elevation) variables

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.2) # AICs
anova(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.2)  #ANOVA F-Test
#the more complex model is preferred, significant ANOVA

#making a dredge model with just smooth elevation
LM_add.gam_SCA.terrain_dist.dredge.just.elev <- gam(Canopy_short ~ s(Elevation..m.FIXED), 
                                                    data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.just.elev)
#the linear fits for distance and slope were not significant

#comparing the models with smoothed distance, elevation, and slope to one with just smoothed elevation and the other variables, and one with only s(elevation)
AIC(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.just.elev)
anova(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.just.elev)
#the more complex model is preferred, significant ANOVA

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_SCA.terrain_dist.dredge 
summary(LM_add.gam_SCA.terrain_dist.dredge)
check_concurvity(LM_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler, so it could be an alternative choice

#making a dredge model with just smooth elevation
LM_add.gam_SCA.terrain_dist.dredge.3 <- gam(Canopy_short ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(Elevation..m.FIXED), 
                                                    data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.3)
#the linear fits for distance and slope were not significant

#making a dredge model with just smooth elevation
LM_add.gam_SCA.terrain_dist.dredge.4 <- gam(Canopy_short ~  LM_slope_raster_15_data_pts + s(LM_Northness) + s(Elevation..m.FIXED), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.4)
#the linear fits for distance and slope were not significant

AIC(LM_add.gam_SCA.terrain_dist.dredge.3, LM_add.gam_SCA.terrain_dist.dredge.4)
AIC(LM_add.gam_SCA.terrain_dist.dredge.3, LM_add.gam_SCA.terrain_dist.dredge)

#while the original dredge model has the lowest AIC the dredge 3 model has the lowest AIC without concurvity

# model without TWI
LM_add.gam_SCA.terrain_dist.dredge.3 <- gam(Canopy_short ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(Elevation..m.FIXED), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.3)
#the linear fits for distance and slope were not significant

# model without elevation
LM_add.gam_SCA.terrain_dist.dredge.no_elev <- gam(Canopy_short ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.no_elev)
#the linear fits for distance and slope were not significant


#comparing the model with TWI vs. with elevation 
AIC(LM_add.gam_SCA.terrain_dist.dredge.3, LM_add.gam_SCA.terrain_dist.dredge.no_elev)
#the model elevation is preferred

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.terrain_dist.dredge.3) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_SCA.terrain_dist.dredge.3)

#Chosen model: LM_add.gam_SCA.terrain_dist.dredge.3

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = 'Slope (º)', ylab = expression(f[1]*'Slope (º)'))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = "LM_Northness", 
         ylab = expression(f[1]*'(LM_Northness)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.terrain_dist.dredge.3, select=3, 
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
LM_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, LM_Northness, LM_slope_raster_15_data_pts), 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):Northness'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_SCA.inter, LM_add.gam_SCA.terrain_dist.dredge)

#overall best model:LM_add.gam_SCA.terrain_dist.dredge.3


## LCA ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_Eastness) + s(LM_TWI_values) + s(X.1, Y.1),   
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#TWI and elevations showing significant smoothness

#checking for concurvity (equivalent to co-linearity but for models with non-linear functions)
concurvity(LM_add.gam_LCA.terrain_dist, full = T)

library(performance)
check_concurvity(LM_add.gam_LCA.terrain_dist)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(elevation), s(Northness), s(slope), s(TWI) as allowing for the best model

# Checking a GAM without aspect
LM_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_Northness) + s(LM_slope_raster_15_data_pts) + s(LM_TWI_values), 
                                             data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(elevation) and s(TWI) are significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.dredge) # AICs
anova(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#checking for concurvity (equivalent to co-linearity but for models with non-linear functions)
concurvity(LM_add.gam_LCA.terrain_dist.dredge, full = T)

check_concurvity(LM_add.gam_LCA.terrain_dist.dredge)

#elevation and TWI show high concurvity

# model without TWI
LM_add.gam_LCA.terrain_dist.dredge.no.twi <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_Northness) + s(LM_slope_raster_15_data_pts), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.no.twi) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(elevation) and s(northness) are significant

# model without elevation
LM_add.gam_LCA.terrain_dist.dredge.no.elev <- gam(Canopy_long ~ s(LM_Northness) + s(LM_slope_raster_15_data_pts) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.no.elev) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(northness) and s(TWI) are significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist.dredge.no.twi, LM_add.gam_LCA.terrain_dist.dredge.no.elev) # AICs

#there is not a significant difference between the model with elevation or with twi

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.terrain_dist.dredge.no.elev) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(d), s(elevation), s(slope) as allowing for the best model

# model without elevation
LM_add.gam_LCA.terrain_dist.dredge.2 <- gam(Canopy_long ~ s(LM_Northness) + s(LM_TWI_values), 
                                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(northness) and s(TWI) are significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist.dredge.no.elev, LM_add.gam_LCA.terrain_dist.dredge.2) # AICs
anova(LM_add.gam_LCA.terrain_dist.dredge.no.elev, LM_add.gam_LCA.terrain_dist.dredge.2)  #ANOVA F-Test
#the more complex model is not preferred, dredge model is preferred, but not significantly

#checking concurvity
check_concurvity(LM_add.gam_LCA.terrain_dist.dredge.no.elev) #low concurvity
check_concurvity(LM_add.gam_LCA.terrain_dist.dredge.2)  #low concurvity

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_LCA.terrain_dist.dredge 
summary(LM_add.gam_LCA.terrain_dist.dredge.no.elev)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.terrain_dist.dredge.no.elev) #pretty normal residuals and no heteroscedasticity 

#looking at significance
summary(LM_add.gam_LCA.terrain_dist.dredge.no.elev)

#Chosen model: LM_add.gam_LCA.terrain_dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.no.elev, select=2, 
         all.terms=T, xlab = 'Slope (º)', ylab = expression(f[1]*'Slope (º)'))
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.no.elev, select=1, 
         all.terms=T, xlab = "LM_Northness", 
         ylab = expression(f[1]*'(LM_Northness)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.terrain_dist.dredge.no.elev, select=3, 
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
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation, eastness, and TWI showing potential for significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(northness), s(slope), s(TWI)

# Checking a GAM without aspect
LM_add.gam_CA.terrain_dist.dredge <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation), s(northness), and somewhat s(TWI) have significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.dredge) # AICs
anova(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.dredge)  #ANOVA F-Test
#the dredged model is preffered

#checking concurvity
check_concurvity(LM_add.gam_CA.terrain_dist.dredge) #low concurvity

#elevation and twi have high concurvity

# model without TWI
LM_add.gam_CA.terrain_dist.no.twi <- gam(Canopy_area_lg ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness), 
                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.no.twi) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) and somewhat s(northness) have significant non-linear functions

# model without elevation
LM_add.gam_CA.terrain_dist.no.elev <- gam(Canopy_area_lg ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.no.elev) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(northness) has significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist.no.elev, LM_add.gam_CA.terrain_dist.no.twi) # AICs

# The model with elevation instead of TWI is significantly stronger

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.terrain_dist.no.twi) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants all of the variables

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CA.terrain_dist.no.twi 
summary(LM_add.gam_CA.terrain_dist.no.twi)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.terrain_dist.no.twi) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_CA.terrain_dist.no.twi)

#Chosen model: LM_add.gam_CA.terrain_dist.no.twi

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.terrain_dist.no.twi, select=2, 
         all.terms=T, xlab = 'Slope (º)', ylab = expression(f[1]*'Slope (º)'))
plot.gam(LM_add.gam_CA.terrain_dist.no.twi, select=3, 
         all.terms=T, xlab = "LM_Northness", 
         ylab = expression(f[1]*'(LM_Northness)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.terrain_dist.no.twi, select=1, 
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
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#northness and TWI showing signs of significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM after the dredge
LM_add.gam_CS.terrain_dist.dredge <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
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
LM_add.gam_CS.terrain_dist.no.twi <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + s(LM_Northness), 
                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.no.twi) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation) and s(northness) have significant non-linear functions

# model without elevation
LM_add.gam_CS.terrain_dist.no.elev <- gam(Crown_spread ~ s(LM_slope_raster_15_data_pts) + s(LM_Northness) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.no.elev) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#s(northness) and s(TWI) have significant non-linear functions

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CS.terrain_dist.no.twi, LM_add.gam_CS.terrain_dist.no.elev) # AICs
#the model without twi is significantly better

#checking concurvity
check_concurvity(LM_add.gam_CS.terrain_dist.no.twi) #low concurvity

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CS.terrain_dist.no.twi 
summary(LM_add.gam_CS.terrain_dist.no.twi)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.terrain_dist.no.twi) #pretty normal residuals and no Heteroscedasticity

#looking at significance
summary(LM_add.gam_CS.terrain_dist.no.twi)

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
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
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
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only elevation is showing signs of significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(eastness), s(twi)

# Checking a GAM without aspect
LM_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~ s(LM_Eastness) + s(LM_slope_raster_15_data_pts) + s(LM_TWI_values), 
                                             data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness and twi showing somehwat significance

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.dredge) # AICs
anova(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.dredge)  #ANOVA F-Test
#the more complicated model is significantly stronger

# Checking a GAM without aspect
LM_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~  s(LM_slope_raster_15_data_pts) + s(LM_Eastness) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness and twi showing somehwat significance

# without twi
LM_add.gam_DBH.terrain_dist.dredge.no.twi <- gam(DBH_ag_lg ~  s(LM_slope_raster_15_data_pts) + s(LM_Eastness) + s(Elevation..m.FIXED), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.no.twi) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness and twi showing somehwat significance

# without elevation
LM_add.gam_DBH.terrain_dist.dredge.no.elev <- gam(DBH_ag_lg ~  s(LM_slope_raster_15_data_pts) + s(LM_Eastness) + s(LM_TWI_values), 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.no.elev) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#eastness and twi showing somehwat significance

AIC(LM_add.gam_DBH.terrain_dist.dredge.no.twi, LM_add.gam_DBH.terrain_dist.dredge.no.elev)
#the one with TWI is preferred over the one with elevation

#checking concurvity
check_concurvity(LM_add.gam_DBH.terrain_dist.dredge) #low concurvity

# while this model is not significantly stronger than the original model, it has low concurvity so I will choose it

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LM_add.gam_DBH.terrain_dist.dredge)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_DBH.terrain_dist.dredge)

#Chosen model: LM_add.gam_DBH.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = 'Eastness', ylab = expression(f[1]*'(Eastness)'))
plot.gam(LM_add.gam_DBH.terrain_dist.dredge, select=3, 
         all.terms=T, xlab = "TWI", 
         ylab = expression(f[1]*'(TWI)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.terrain_dist.dredge, select=1, 
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
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation) and aspect

# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short_lg ~ s(d) + s(LC_Eastness) + s(X.1, Y.1), 
                                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
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
check_concurvity(LC_add.gam_SCA.terrain_dist.dredge.2) #low concurvity


# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.dredge.2 <- gam(Canopy_short_lg ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + s(LC_TWI_values), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has a significant non-linear fit

#checking the significant differences between models 
AIC(LC_add.gam_SCA.terrain_dist.dredge, LC_add.gam_SCA.terrain_dist.dredge.2)

# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.dredge.4 <- gam(Canopy_short_lg ~ s(d) + s(LC_Eastness) + s(LC_TWI_values) + s(LC_slope_raster_15_data_pts) , 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#location has a significant non-linear fit
check_concurvity(LC_add.gam_SCA.terrain_dist.dredge.4)

#checking the AICs
AIC(LC_add.gam_SCA.terrain_dist.dredge, LC_add.gam_SCA.terrain_dist.dredge.2, LC_add.gam_SCA.terrain_dist.dredge.4) # AICs

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_SCA.terrain_dist.dredge 
summary(LC_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_SCA.terrain_dist.dredge$residuals)

#looking at significance
summary(LC_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LC_add.gam_SCA.terrain_dist.dredge

#plotting all of the variables for observation
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge.3, select=2, 
         all.terms=T, xlab = 'Eastness', ylab = expression(f[1]*'Eastness'))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge.3, select=1, 
         all.terms=T, xlab = "distance to river", 
         ylab = expression(f[1]*'(distance to river)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.terrain_dist.dredge.3, select=3, 
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
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
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
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, elevation, and slope do not have significant non-linear function, aspect is significant linearly for W

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation)

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.dredge <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_Eastness) + s(LC_Northness), 
                                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and eastness is significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge) # AICs
anova(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants all four variables

#check concurvity
check_concurvity(LC_add.gam_LCA.terrain_dist.dredge)
# eachness and northness have high concurvity

# with eastness
LC_add.gam_LCA.terrain_dist.dredge.east <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_Eastness), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.east) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and eastness is significant

# with northness
LC_add.gam_LCA.terrain_dist.dredge.north <- gam(Canopy_long_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_Northness), 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.north) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#elevation and eastness is significant

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.dredge, LC_add.gam_LCA.terrain_dist.dredge.east, LC_add.gam_LCA.terrain_dist.dredge.north)
#LC_add.gam_LCA.terrain_dist.dredge.east has the lowest AIC values and no concurivty

#checking concurvity
check_concurvity(LC_add.gam_LCA.terrain_dist.dredge.east) #low concurvity

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_LCA.terrain_dist.dredge 
summary(LC_add.gam_LCA.terrain_dist.dredge.east)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.terrain_dist.dredge.east) #pretty normal residuals and no Heteroscedasticity 
shapiro.test(LC_add.gam_LCA.terrain_dist.dredge.east$residuals)

#looking at significance
summary(LC_add.gam_LCA.terrain_dist.dredge.east)

#Chosen model: LC_add.gam_LCA.terrain_dist.dredge.east

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.east, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.east, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.terrain_dist.dredge.east, select=3, 
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
AIC(LC_add.gam_LCA.inter, LC_add.gam_LCA.terrain_dist)

#overall best model: LC_add.gam_LCA.terrain_dist

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#location has significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(location), s(eastness)

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge <- gam(Canopy_area_lg ~ s(d) + s(X.1, Y.1) + s(LC_Eastness), 
                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.dredge) # AICs
anova(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(location), s(eastness)

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.2 <- gam(Canopy_area_lg ~ s(X.1, Y.1) + s(Elevation..m.FIXED), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.2)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit


# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge.3 <- gam(Canopy_area_lg ~ s(d) + s(LC_Eastness), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)

# Checking a GAM after the dredge
LC_add.gam_CA.terrain_dist.dredge.4 <- gam(Canopy_area_lg ~ s(d) + s(LC_TWI_values) + s(LC_Eastness), 
                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist.dredge, LC_add.gam_CA.terrain_dist.dredge.2, LC_add.gam_CA.terrain_dist.dredge.3, LC_add.gam_CA.terrain_dist.dredge.4) # AICs
#the dredge model with s(d), s(eastness), s(location) is preferred

#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_CA.terrain_dist.dredge) #storing the residuals of the GAM
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

# the model without s(x,y) does not meet the control of spatial autocorrelation, so we need s(x,y)
# using gam.check, we observed that we do not need to adjust the k value

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CA.terrain_dist.dredge 
summary(LC_add.gam_CA.terrain_dist.dredge)
#but all models seem to do similarly well

#check concurvity
check_concurvity(LC_add.gam_CA.terrain_dist.dredge) #low concurvity

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
LC_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + 
                                    s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#no variables have significant non-linear functions

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge <- gam(Crown_spread ~ s(d) + s(LC_Eastness) + s(X.1, Y.1), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(eastness) is borderline significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.2 <- gam(Crown_spread ~ s(LC_Eastness) + s(X.1, Y.1), 
                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.2) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(eastness) is borderline significant

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.3 <- gam(Crown_spread ~ s(LC_Eastness) + s(Elevation..m.FIXED), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.3) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(eastness) is borderline significant

# Checking a GAM post dredge
LC_add.gam_CS.terrain_dist.dredge.4 <- gam(Crown_spread ~ s(d) + s(LC_Eastness) + s(Elevation..m.FIXED), 
                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.4) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(eastness) is borderline significant

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.dredge, LC_add.gam_CS.terrain_dist.dredge.2, LC_add.gam_CS.terrain_dist.dredge.3, LC_add.gam_CS.terrain_dist.dredge.4) # AICs
anova(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.dredge)  #ANOVA F-Test
#the more complicated model is not significantly better
# the models are similarly strong, and we do not need to control for spatial autocorrelation
# I will use dredge 4

#LC_add.gam_CS.terrain_dist.dredge is currently the best model


#checking if the model without location meets spatial autocorrelation 

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_CS.terrain_dist.dredge.4) #storing the residuals of the GAM
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

#none of the models are having trouble with spatial autocorrelation
# using gam.check, we observed that we do not need to adjust the k value


#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CS.terrain_dist.no.aspect 
summary(LC_add.gam_CS.terrain_dist.dredge.4)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.terrain_dist.dredge.4) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_CS.terrain_dist.dredge.4)

#Chosen model: LC_add.gam_CS.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CS.terrain_dist.dredge.4, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CS.terrain_dist.dredge.4, select=3, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CS.terrain_dist.dredge.4, select=2, 
         all.terms=T, xlab = "Eastness (º)", ylab = expression(f[1]*'Eastness'), 
         se = TRUE , col = "black")


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_add.gam_CS.terrain_dist.dredge.4$Elevation..m.FIXED, 
        y=LC_add.gam_CS.terrain_dist.dredge.4$d, 
        z=LC_add.gam_CS.terrain_dist.dredge.4$LC_Eastness, 
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, d, LC_Eastness), 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_CS.inter, LC_add.gam_CS.terrain_dist.no.aspect)

#overall best model: LC_add.gam_CS.terrain_dist.no.aspect

## DBH ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_DBH.terrain_dist <- gam(DBH_ag_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + 
                                     s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only location has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(location)

# Checking a GAM without aspect
LC_add.gam_DBH.terrain_dist.dredge <- gam(DBH_ag_lg ~ s(d) + s(X.1,Y.1), 
                                                 data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(location) is significant 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_DBH.terrain_dist.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.dredge) # AICs
anova(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.dredge)  #ANOVA F-Test
#the dredge model is preferred

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.2 <- gam(DBH_ag_lg ~ s(d),
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.2)

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.3 <- gam(DBH_ag_lg ~ s(d) + s(X.1,Y.1) + s(LC_Eastness) + s(LC_TWI_values), 
                                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.3)
#s(distance) is marginally not significantly  non-linear, slope and elevation are not significantly useful with the linear fit

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.4 <- gam(DBH_ag_lg ~ s(d) + s(LC_Eastness) + s(LC_TWI_values), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.4)
#s(eastiness) is significantly non-linear 

#comparing the models 
AIC(LC_add.gam_DBH.terrain_dist.dredge, LC_add.gam_DBH.terrain_dist.dredge.2, LC_add.gam_DBH.terrain_dist.dredge.3, LC_add.gam_DBH.terrain_dist.dredge.4)
#LC_add.gam_DBH.terrain_dist.dredge.3 has lowest AIC

# checking to see if we need to control for spatial autocorrelation

#creating the matrix of coordinates to evaluate if spatial autocorrelation as control
res <- residuals(LC_add.gam_DBH.terrain_dist.dredge.4) #storing the residuals of the GAM
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

#LC_add.gam_DBH.terrain_dist.dredge controls for  spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.2 does not control for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.3 controls for spatial autocorrelation
#LC_add.gam_DBH.terrain_dist.dredge.4 does not for spatial autocorrelation


#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_DBH.terrain_dist.dredge.3 
summary(LC_add.gam_DBH.terrain_dist.dredge.3)
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
         all.terms=T, xlab = "Eastness (º)", ylab = expression(f[1]*'Eastness'), 
         se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.terrain_dist.dredge.3, select=4, 
         all.terms=T, xlab = "TWI", ylab = expression(f[1]*'TWI'), 
         se = TRUE , col = "black")

# # 3d plotting in plotly and with gg3D
# plot_ly(x=LC_add.gam_DBH.terrain_dist.dredge.3$Elevation..m.FIXED, 
#         y=LC_add.gam_DBH.terrain_dist.dredge.3$X., 
#         z=LC_add.gam_DBH.terrain_dist.dredge.3$LC_Eastness, 
#         color = LC_add.gam_DBH.terrain_dist.dredge.3$LC_TWI_values,
#         type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_DBH.inter <- gam(DBH_ag ~ ti(LC_Eastness, X.1, Y.1, d, LC_TWI_values), 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_DBH.inter, LC_add.gam_DBH.terrain_dist)
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
SD_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + s(LC_Northness) + 
                                     s(LC_Eastness) + s(LC_TWI_values) + s(X.1, Y.1), 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope has significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM with no aspect
SD_add.gam_SCA.terrain_dist.no.aspect <- gam(Canopy_short ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                             data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#slope of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_SCA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect) # AICs
anova(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model without aspect is better

# setting up the dredge GAM with only slope being smoothed
SD_add.gam_SCA.terrain_dist.slope.smooth <- gam(Canopy_short ~ d + Elevation..m. + s(SD_slope_raster_15_data_pts), 
                                                data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_SCA.terrain_dist.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants elevation, s(slope)

# setting up the dredge GAM with only slope being smoothed, no distance
SD_add.gam_SCA.terrain_dist.slope.smooth.dredge <- gam(Canopy_short ~ Elevation..m. + s(SD_slope_raster_15_data_pts), 
                                                       data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_SCA.terrain_dist.slope.smooth.dredge)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_SCA.terrain_dist.slope.smooth.dredge) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the models 
AIC(SD_add.gam_SCA.terrain_dist, SD_add.gam_SCA.terrain_dist.no.aspect,  SD_add.gam_SCA.terrain_dist.slope.smooth, SD_add.gam_SCA.terrain_dist.slope.smooth.dredge)
summary(SD_add.gam_SCA.terrain_dist.no.aspect)
#SD_add.gam_SCA.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_SCA.terrain_dist.dredge 
summary(SD_add.gam_SCA.terrain_dist.no.aspect)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.terrain_dist) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_SCA.terrain_dist.no.aspect)

#Chosen model: SD_add.gam_SCA.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_SCA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_SCA.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
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
SD_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function and sw was significant aspect

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(slope) and aspect

# Checking a GAM with only slope smooth
SD_add.gam_LCA.terrain_dist.slope.smooth <- gam(Canopy_long ~ d + Elevation..m. + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.slope.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only west is significant
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_LCA.terrain_dist.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#making a dredge model with just smooth elevation
SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist <- gam(Canopy_long ~ s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                       data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)
#s(slope) and NE and SW are significant

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(SD_add.gam_LCA.terrain_dist, SD_add.gam_LCA.terrain_dist.slope.smooth, SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)
#SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist 
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist) #pretty normal residuals and no Heteroscedasticity 

#looking at model significance
summary(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#Chosen model: SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist, select=1, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_LCA.inter, SD_add.gam_LCA.terrain_dist.dredge.no.elev.dist)

#overall best model: SD_add.gam_LCA.inter

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
SD_add.gam_CA.terrain_dist.no.aspect <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                            data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect) # AICs
anova(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and slope and no aspect
SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth <- gam(Canopy_area_lg ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                                data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#making a dredge model with just smooth elevation and slope
SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist <- gam(Canopy_area_lg ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                                        data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(SD_add.gam_CA.terrain_dist, SD_add.gam_CA.terrain_dist.no.aspect, SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth, SD_add.gam_CA.terrain_dist.dredge.just.elev.slope.smooth.no.dist)
#SD_add.gam_CA.terrain_dist.no.aspect  has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CA.terrain_dist.no.aspect 
summary(SD_add.gam_CA.terrain_dist.no.aspect)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.terrain_dist.no.aspect) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_CA.terrain_dist.no.aspect)

#Chosen model: SD_add.gam_CA.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_CA.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_CA.inter, SD_add.gam_CA.terrain_dist.no.aspect)

#overall best model: SD_add.gam_CA.terrain_dist.no.aspect

## CS ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#slope have significant non-linear function and Southwest are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope), and aspect

# Checking a GAM without d smoothed
SD_add.gam_CS.terrain_dist.no.dist.smooth <- gam(Crown_spread ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                 data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.no.dist.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation and slope) are significant
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_CS.terrain_dist.no.dist.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope), aspect

#making a dredge model with just smooth elevation and no aspect
SD_add.gam_CS.terrain_dist.no.dist <- gam(Crown_spread ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                          data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_CS.terrain_dist.no.dist)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models
AIC(SD_add.gam_CS.terrain_dist, SD_add.gam_CS.terrain_dist.no.dist.smooth, SD_add.gam_CS.terrain_dist.no.dist)
#SD_add.gam_CS.terrain_dist.no.dist has lowest AIC, but none are significantlt better than the others

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_CS.terrain_dist.no.aspect 
summary(SD_add.gam_CS.terrain_dist.no.dist)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.terrain_dist.no.dist) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_CS.terrain_dist.no.dist)

#Chosen model: SD_add.gam_CS.terrain_dist.no.dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CS.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_CS.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CS.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_CS.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                           data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_CS.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_CS.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_CS.inter, SD_add.gam_CS.terrain_dist)

#overall best model: SD_add.gam_CS.terrain_dist


## DBH ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
SD_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                   data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only slope have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

# Checking a GAM without smoothing distance 
SD_add.gam_DBH.terrain_dist.no.dist.smooth <- gam(DBH_ag ~ d + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.no.dist.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(slope) is significant
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_DBH.terrain_dist.no.dist.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#making a dredge model with just smooth elevation and no aspect
SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth <- gam(DBH_ag ~ s(Elevation..m.) + s(SD_slope_raster_15_data_pts), 
                                                                 data = SD_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#s(slope) is significantly non-linear but not s(elevation), slope and distance are not significantly useful with the linear fit
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the models 
AIC(SD_add.gam_DBH.terrain_dist, SD_add.gam_DBH.terrain_dist.no.dist.smooth, SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: SD_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(SD_add.gam_DBH.terrain_dist.dredge.no.smooth)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)

#Chosen model: SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth, select=1, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(SD_add.gam_DBH.terrain_dist, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m., 
        y=SD_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=SD_fixed_field_data_processed_terrain_dist_no_NA$SD_slope_raster_15_data_pts, 
        color = SD_fixed_field_data_processed_terrain_dist_no_NA$SD_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
SD_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m., d, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                            data = SD_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(SD_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(SD_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(SD_add.gam_DBH.inter, SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth)
#the interaction is not significant

#overall best model: SD_add.gam_DBH.terrain_dist.dredge.just.elev.slope.smooth 




