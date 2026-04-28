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

#removing the spatial geometry to be able to use the GAM function
LM_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LM_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation has a significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(d), s(elevation), s(slope) as allowing for the best model

# Checking a GAM with just elevation smoothing splines s()
LM_add.gam_SCA.terrain_dist.just.elevation.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model model, rule of thumb is that a difference of 2 is a significant difference
# it wants only s(elevation)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.just.elevation.smooth) # AICs
anova(LM_add.gam_SCA.terrain_dist, LM_add.gam_SCA.terrain_dist.just.elevation.smooth)  #ANOVA F-Test
#the more complex model with smoothing splines on all three quantitative variables is preferable

# setting up the dredge GAM of the larger equation, with just quantitative variables
LM_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #NA fail makes sure the later dredge does not have to worry about NAs

#comparing model with smoothing with and without aspect
AIC(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist)
summary(LM_add.gam_SCA.terrain_dist.dredge)
#aspect does not appear to be necessary

#making a GAM with just elevation smoothed, no aspect
LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model best, rule of thumb is that a difference of 2 is a significant difference
#just s(elevation is preferred)

#making a dredge model with just smooth elevation
LM_add.gam_SCA.terrain_dist.dredge.just.elev <- gam(Canopy_short ~ s(Elevation..m.FIXED), 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth)
#the linear fits for distance and slope were not significant

#comparing the models with smoothed distance, elevation, and slope to one with just smoothed elevation and the other variables, and one with only s(elevation)
AIC(LM_add.gam_SCA.terrain_dist.dredge, LM_add.gam_SCA.terrain_dist.dredge.just.elev.smooth, LM_add.gam_SCA.terrain_dist.dredge.just.elev)
#LM_add.gam_SCA.terrain_dist.dredge has lowest AIC, the one with smoothing splines on all quantiative variables and no aspect

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_SCA.terrain_dist.dredge 
summary(LM_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler, so it could be an alternative choice

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LM_add.gam_SCA.terrain_dist.dredge

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.terrain_dist.dredge, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_SCA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_SCA.inter, LM_add.gam_SCA.terrain_dist.dredge)

#overall best model:LM_add.gam_SCA.terrain_dist.dredge


## LCA ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, elevation, and slope have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(d), s(elevation), s(slope) as allowing for the best model

# Checking a GAM without aspect
LM_add.gam_LCA.terrain_dist.no.aspect <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(d), s(elevation), s(slope) as allowing for the best model

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation
LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_long ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                    data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth)

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_LCA.terrain_dist, LM_add.gam_LCA.terrain_dist.no.aspect, LM_add.gam_LCA.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_LCA.terrain_dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_LCA.terrain_dist.dredge 
summary(LM_add.gam_LCA.terrain_dist)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.terrain_dist) #pretty normal residuals and no heteroscedasticity 

#looking at significance
summary(LM_add.gam_LCA.terrain_dist)

#Chosen model: LM_add.gam_LCA.terrain_dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_LCA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_LCA.inter <- gam(Canopy_long ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                            data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_LCA.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_LCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LM_add.gam_LCA.inter, LM_add.gam_LCA.terrain_dist)

#overall best model: LM_add.gam_LCA.terrain_dist

## CA ##

# I am using the logged transformation of canopy area to get more normal residuals and less heteroscedasticity to better meet the conditions for the GAMs

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_CA.terrain_dist.no.aspect <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                             data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth <- gam(Canopy_area_lg ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_CA.terrain_dist, LM_add.gam_CA.terrain_dist.no.aspect, LM_add.gam_CA.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_CA.terrain_dist has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CA.terrain_dist.no.aspect 
summary(LM_add.gam_CA.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.terrain_dist.no.aspect) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_CA.terrain_dist.no.aspect)

#Chosen model: LM_add.gam_CA.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_CA.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
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
LM_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, slope, elevation have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_CS.terrain_dist.no.aspect <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CS.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth <- gam(Crown_spread ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts, 
                                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_CS.terrain_dist, LM_add.gam_CS.terrain_dist.no.aspect, LM_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#LM_add.gam_CS.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_CS.terrain_dist.no.aspect 
summary(LM_add.gam_CS.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.terrain_dist.no.aspect) #pretty normal residuals and no Heteroscedasticity

#looking at significance
summary(LM_add.gam_CS.terrain_dist.no.aspect)

#Chosen model: LM_add.gam_CS.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CS.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_CS.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
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

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                  data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#only elevation have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LM_add.gam_DBH.terrain_dist.no.aspect <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                            data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_DBH.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect) # AICs
anova(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is significantly better

#making a dredge model with just smooth elevation and no aspect
LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth <- gam(DBH_ag ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                          data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit

#making a model with no smoothing splines
LM_add.gam_DBH.terrain_dist.dredge.no.smooth <- gam(DBH_ag ~ d + Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LM_add.gam_DBH.terrain_dist, LM_add.gam_DBH.terrain_dist.no.aspect, LM_add.gam_DBH.terrain_dist.dredge.just.elev.smooth, LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#LM_add.gam_DBH.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LM_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)
#the linear model seems to do the best, the GAM model with all smoothing is close behind

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.terrain_dist.dredge.no.smooth) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LM_add.gam_DBH.terrain_dist.dredge.no.smooth)

#Chosen model: LM_add.gam_DBH.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LM_add.gam_DBH.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LM_add.gam_DBH.terrain_dist, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LM_fixed_field_data_processed_terrain_dist_no_NA$LM_slope_raster_15_data_pts, 
        color = LM_fixed_field_data_processed_terrain_dist_no_NA$LM_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")


#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LM_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m.FIXED, d, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                           data = LM_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LM_add.gam_DBH.inter)
#there was no significant interaction term after checking combinations

#interaction plots
plot.gam(LM_add.gam_DBH.inter, select=1, 
         all.terms=T, main = "s(Elevation:Slope:clay content)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º)):clay content (‰)'), se = TRUE,
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
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #Crown Spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed

## SCA ##

#removing the spatial geometry to be able to use the GAM function
LC_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LC_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#none of the variables has significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation) and aspect

# Checking a GAM with just elevation smoothing splines s()
LC_add.gam_SCA.terrain_dist.just.elevation.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                         data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_SCA.terrain_dist.just.elevation.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.terrain_dist.just.elevation.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants only s(elevation) and aspect

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth) # AICs
anova(LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth)  #ANOVA F-Test
#the more complex model with smoothing splines on all three quantitative variables is preferable

# setting up the dredge model aspect and s(elevation)
LC_add.gam_SCA.terrain_dist.dredge <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs

#comparing model with smoothing with and without aspect
AIC(LC_add.gam_SCA.terrain_dist.dredge, LC_add.gam_SCA.terrain_dist, LC_add.gam_SCA.terrain_dist.just.elevation.smooth)
#LC_add.gam_SCA.terrain_dist.dredge has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_SCA.terrain_dist.dredge 
summary(LC_add.gam_SCA.terrain_dist.dredge)
#but a model with just s(Elevation) seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_SCA.terrain_dist.dredge)

#Chosen model: LC_add.gam_SCA.terrain_dist.dredge

#plotting the chosen function, with no interaction  LC_add.gam_SCA.terrain_dist.dredge
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist.dredge, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
visreg(LC_add.gam_SCA.terrain_dist.dredge, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA")

#plotting all of the variables for observation
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_SCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_SCA.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on SCA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_SCA.inter <- gam(Canopy_short ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                            data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_SCA.inter)
#there was a significant interaction term 

#interaction plots
plot.gam(LC_add.gam_SCA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_SCA.inter, LC_add.gam_SCA.terrain_dist.dredge)

#overall best model: LC_add.gam_SCA.inter

## LCA ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_LCA.terrain_dist <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, elevation, and slope do not have significant non-linear function, aspect is significant linearly for W

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation)

# Checking a GAM without aspect
LC_add.gam_LCA.terrain_dist.no.aspect.slope <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED), 
                                             data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.no.aspect.slope) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.terrain_dist.no.aspect.slope) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope) # AICs
anova(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope)  #ANOVA F-Test
#the model with aspect and slope is significantly better

#making a dredge model with just smooth elevation
LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth <- gam(Canopy_long ~ s(d) + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth)

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_LCA.terrain_dist, LC_add.gam_LCA.terrain_dist.no.aspect.slope, LC_add.gam_LCA.terrain_dist.dredge.just.elev.d.smooth)
#none of the models have much better AIC values 

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_LCA.terrain_dist.dredge 
summary(LC_add.gam_LCA.terrain_dist)
#but a model with aspect seems like it can do similarly as well and is simpler

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.terrain_dist) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_LCA.terrain_dist)

#Chosen model: LC_add.gam_LCA.terrain_dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_LCA.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_LCA.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on LCA") 

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
LC_add.gam_CA.terrain_dist <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), aspect

# Checking a GAM without s(slope) but slope instead
LC_add.gam_CA.terrain_dist.less.s <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.less.s) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.terrain_dist.less.s) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), aspect

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s) # AICs
anova(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s)  #ANOVA F-Test
#the model with s(slope) is significantly better

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CA.terrain_dist.dredge.no.slope <- gam(Canopy_area_lg ~ s(d) + s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_CA.terrain_dist, LC_add.gam_CA.terrain_dist.less.s, LC_add.gam_CA.terrain_dist.dredge.no.slope)
#LC_add.gam_CA.terrain_dist.dredge.no.slope has lowest AIC, but really none are much better

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CA.terrain_dist.dredge.no.slope 
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.terrain_dist.dredge.no.slope) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_CA.terrain_dist.dredge.no.slope)

#Chosen model: LC_add.gam_CA.terrain_dist.dredge.no.slope

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.terrain_dist.dredge.no.slope, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CA.terrain_dist.dredge.no.slope, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_CA.terrain_dist.dredge.no.slope, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_CA.inter <- gam(Canopy_area ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                           data = LC_fixed_field_data_processed_terrain_dist_no_NA,  na.action = na.fail)
summary(LC_add.gam_CA.inter)
#there was is a significant interaction term after checking combinations

#interaction plots
plot.gam(LC_add.gam_CA.inter, select=1, 
         all.terms=T, main = "s(Elevation:Distance:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Distance (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)

#checking to see whether interaction model outperforms our previously selected model
AIC(LC_add.gam_CA.inter, LC_add.gam_CA.terrain_dist.no.aspect)

#overall best model: LC_add.gam_CA.terrain_dist.no.aspect

## CS ##

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LC_add.gam_CS.terrain_dist <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                  data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#distance, slope, elevation have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

# Checking a GAM without aspect
LC_add.gam_CS.terrain_dist.no.aspect <- gam(Crown_spread ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts), 
                                            data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.no.aspect) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.terrain_dist.no.aspect) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d), s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect) # AICs
anova(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect)  #ANOVA F-Test
#the model with aspect is not significantly better

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth <- gam(Crown_spread ~ d + s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts, 
                                                          data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#s(elevation) is significantly non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models with smoothed distance, elevation, and slope with and without aspect, and one with only s(elevation) and no aspect
AIC(LC_add.gam_CS.terrain_dist, LC_add.gam_CS.terrain_dist.no.aspect, LC_add.gam_CS.terrain_dist.dredge.just.elev.smooth)
#LC_add.gam_CS.terrain_dist.no.aspect has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_CS.terrain_dist.no.aspect 
summary(LC_add.gam_CS.terrain_dist.no.aspect)
#but all models seem to do similarly well

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.terrain_dist.no.aspect) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_CS.terrain_dist.no.aspect)

#Chosen model: LC_add.gam_CS.terrain_dist.no.aspect

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CS.terrain_dist.no.aspect, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_CS.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 


# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_CS.inter <- gam(Crown_spread ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
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
LC_add.gam_DBH.terrain_dist <- gam(DBH_ag ~ s(d) + s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                   data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#none have significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_DBH.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(d)

# Checking a GAM without aspect
LC_add.gam_DBH.terrain_dist.only.distance <- gam(DBH_ag ~ s(d), 
                                             data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.only.distance) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#only s(elevation is significant)
#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_DBH.terrain_dist.only.distance) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# it wants s(elevation), s(slope)

#comparing the AIC of the model by comparing their AICs and using an ANOVA F-Test
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance) # AICs
anova(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance)  #ANOVA F-Test
#the model with only distance is sufficient

#making a dredge model with just smooth elevation and no aspect
LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth <- gam(DBH_ag ~ s(d) + Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                           data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth)
#s(distance) is marginally not significantly  non-linear, slope and elevation are not significantly useful with the linear fit

#making a model with no smoothing splines
LC_add.gam_DBH.terrain_dist.dredge.no.smooth <- gam(DBH_ag ~ d + Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                    data = LC_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LC_add.gam_DBH.terrain_dist.dredge.no.smooth)
#s(elevation) is marginally not significantly  non-linear, slope and distance are not significantly useful with the linear fit

#comparing the models 
AIC(LC_add.gam_DBH.terrain_dist, LC_add.gam_DBH.terrain_dist.only.distance, LC_add.gam_DBH.terrain_dist.dredge.just.dist.smooth, LC_add.gam_DBH.terrain_dist.dredge.no.smooth)
#LC_add.gam_DBH.terrain_dist.dredge.no.smooth has lowest AIC

#Based on the comparisons (AIC/Anova) of these models, the best model seems to be: LC_add.gam_DBH.terrain_dist.dredge.no.smooth 
summary(LC_add.gam_DBH.terrain_dist)
#the model with just distance smooth is marginally better, so I will use the one with smoothing splines on all of the quantitative variables

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.terrain_dist) #pretty normal residuals and no Heteroscedasticity 

#looking at significance
summary(LC_add.gam_DBH.terrain_dist)

#Chosen model: LC_add.gam_DBH.terrain_dist

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_DBH.terrain_dist, select=1, 
         all.terms=T, xlab = 'Distance (m)', ylab = expression(f[1]*'(Distance)'))
plot.gam(LC_add.gam_DBH.terrain_dist, select=2, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.terrain_dist, select=3, 
         all.terms=T, xlab = "Slope (º)", ylab = expression(f[1]*'Slope'), 
         se = TRUE , col = "black")
visreg(LC_add.gam_DBH.terrain_dist, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on CA") 

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_dist_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_dist_no_NA$d, 
        z=LC_fixed_field_data_processed_terrain_dist_no_NA$LC_slope_raster_15_data_pts, 
        color = LC_fixed_field_data_processed_terrain_dist_no_NA$LC_aspect_raster_15_data_pts_8_categorical,
        type="scatter3d", mode="markers")

#checking for significant interaction terms  

#creating an interaction model using tensor interaction to get interaction smooths
LC_add.gam_DBH.inter <- gam(DBH_ag ~ ti(Elevation..m.FIXED, d, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
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
SD_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
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
gam.check(SD_add.gam_SCA.terrain_dist.no.aspect) #pretty normal residuals and no Heteroscedasticity 

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


#### Session Info ####
# 
# R version 4.4.3 (2025-02-28)
# Platform: aarch64-apple-darwin20
# Running under: macOS Sequoia 15.2
# 
# Matrix products: default
# BLAS:   /System/Library/Frameworks/Accelerate.framework/Versions/A/Frameworks/vecLib.framework/Versions/A/libBLAS.dylib 
# LAPACK: /Library/Frameworks/R.framework/Versions/4.4-arm64/Resources/lib/libRlapack.dylib;  LAPACK version 3.12.0
# 
# locale:
#   [1] en_US.UTF-8/en_US.UTF-8/en_US.UTF-8/C/en_US.UTF-8/en_US.UTF-8
# 
# time zone: America/New_York
# tzcode source: internal
# 
# attached base packages:
#   [1] stats     graphics  grDevices utils     datasets  methods   base     
# 
# other attached packages:
#   [1] lmtest_0.9-40          zoo_1.8-15             gdalUtilities_1.2.5   
# [4] car_3.1-5              carData_3.0-6          googledrive_2.1.1     
# [7] spdep_1.3-10           spData_2.3.4           performance_0.16.0    
# [10] ggsoiltexture_1.1.1    whitebox_2.4.3         spatialEco_2.0-3      
# [13] visreg_2.7.0           MuMIn_1.48.11          plotly_4.10.4         
# [16] mgcv_1.9-1             tmaptools_3.2          geostatsp_2.0.8       
# [19] terra_1.8-29           Matrix_1.7-2           starsExtra_0.2.8      
# [22] stars_0.6-8            abind_1.4-8            ggnewscale_0.5.1      
# [25] rstatix_0.7.2          raster_3.6-32          sp_2.2-1              
# [28] spatstat_3.3-1         spatstat.linnet_3.2-5  spatstat.model_3.3-4  
# [31] rpart_4.1.24           spatstat.explore_3.5-2 nlme_3.1-167          
# [34] spatstat.random_3.4-1  spatstat.geom_3.5-0    spatstat.univar_3.1-4 
# [37] spatstat.data_3.1-6    geomtextpath_0.1.5     PMCMRplus_1.9.12      
# [40] ggpmisc_0.6.1          ggpp_0.5.8-1           smatr_3.4-8           
# [43] sf_1.0-24              moments_0.14.1         lubridate_1.9.4       
# [46] forcats_1.0.0          stringr_1.6.0          dplyr_1.2.0           
# [49] purrr_1.2.1            readr_2.2.0            tidyr_1.3.2           
# [52] tibble_3.3.1           ggplot2_4.0.2          tidyverse_2.0.0       
# [55] emmeans_2.0.2         
# 
# loaded via a namespace (and not attached):
#   [1] RColorBrewer_1.1-3    wk_0.9.4              rstudioapi_0.18.0    
# [4] jsonlite_2.0.0        magrittr_2.0.4        TH.data_1.1-5        
# [7] estimability_1.5.1    spatstat.utils_3.1-5  SuppDists_1.1-9.8    
# [10] farver_2.1.2          fs_1.6.7              vctrs_0.7.2          
# [13] memoise_2.0.1         usethis_3.2.1         htmltools_0.5.9      
# [16] polynom_1.4-1         curl_6.2.1            broom_1.0.12         
# [19] s2_1.1.7              BWStest_0.2.3         Formula_1.2-5        
# [22] KernSmooth_2.23-26    htmlwidgets_1.6.4     sandwich_3.1-1       
# [25] cachem_1.1.0          igraph_2.2.2          lifecycle_1.0.5      
# [28] pkgconfig_2.0.3       R6_2.6.1              fastmap_1.2.0        
# [31] digest_0.6.39         numDeriv_2016.8-1.1   colorspace_2.1-2     
# [34] tensor_1.5            pkgload_1.4.1         nngeo_0.4.8          
# [37] textshaping_1.0.0     labeling_0.4.3        lwgeom_0.2-14        
# [40] spatstat.sparse_3.1-0 timechange_0.3.0      httr_1.4.7           
# [43] polyclip_1.10-7       compiler_4.4.3        gargle_1.5.2         
# [46] remotes_2.5.0         proxy_0.4-27          withr_3.0.2          
# [49] backports_1.5.0       DBI_1.2.3             pkgbuild_1.4.8       
# [52] MASS_7.3-64           quantreg_6.1          sessioninfo_1.2.3    
# [55] classInt_0.4-11       tools_4.4.3           units_0.8-6          
# [58] goftest_1.2-3         glue_1.8.0            grid_4.4.3           
# [61] generics_0.1.4        gtable_0.3.6          tzdb_0.5.0           
# [64] class_7.3-23          data.table_1.17.0     hms_1.1.4            
# [67] pillar_1.11.1         splines_4.4.3         lattice_0.22-9       
# [70] survival_3.8-3        gmp_0.7-5             deldir_2.0-4         
# [73] SparseM_1.84-2        tidyselect_1.2.1      stats4_4.4.3         
# [76] devtools_2.4.6        stringi_1.8.7         boot_1.3-31          
# [79] lazyeval_0.2.2        codetools_0.2-20      kSamples_1.2-10      
# [82] multcompView_0.1-10   cli_3.6.5             xtable_1.8-4         
# [85] systemfonts_1.2.3     munsell_0.5.1         dichromat_2.0-0.1    
# [88] Rcpp_1.1.1            coda_0.19-4.1         XML_3.99-0.20        
# [91] parallel_4.4.3        ellipsis_0.3.2        MatrixModels_0.5-4   
# [94] Rmpfr_1.0-0           viridisLite_0.4.3     mvtnorm_1.3-6        
# [97] scales_1.4.0          e1071_1.7-16          insight_1.4.6        
# [100] crayon_1.5.3          rlang_1.1.7           multcomp_1.4-30  
# 


