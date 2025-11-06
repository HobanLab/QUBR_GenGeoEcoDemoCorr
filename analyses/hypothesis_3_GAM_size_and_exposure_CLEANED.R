# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei size/shape relate to their exposure%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size and shape of Quercus brandegeei 
# individuals across all sites is impacted by their slope and/or elevation and/or aspect 
# using generalized additive models (GAMs)

# We used Generalized Additive Models (GAMs) after having trouble with multiple linear 
# regressions because of issues with the normality conditions and some nervousness about linearity. 

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data includes:
      # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
      # Extracting and processing slope, elevation, and aspect (4 and 8 cardinal directions) data using 15 m res rasters,
# 2) Creating Generalized Additive Models to look at the relationships between the shape/size of trees and their aspect, elevation, and slope

# NOTE: Uncomment and run line 49, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) #ggplot extension
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(car) #to create added variable plots and to run levene's test for checking ANOVA conditions
library(stars) # to convert raster into stars
library(gdalUtilities) #to be able to use gdalwarp
library(mgcv) #needed for the gam function for generalized additive models
library(starsExtra) #to use dist_to_nearest
library(MuMIn) #to use the dredge function
library(rpart) #to use the function rpart to check recurissive binary
library(visreg) #package to be able to plot effects of categorical variables
library(gratia) #using the function smooth estimates
library(gridExtra) #way to arrange ggplots in one plot
library(plotly) #3d plotting
devtools::install_github("AckerDWM/gg3D") #3d plotting
library("gg3D") #3d plotting
library(mgcViz) #3d plotting
library(rgl) #3d plotting

# loading in the processed tree data 
# NOTE: Uncomment and run line 49, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")

#### Generalized Additive Models ####

# For each population/all populations and each size/shape metric we created
     #generalized additive models whereby...

    # a) We removed NAs from the explanatory and response variables
    # b) We then created four base GAM models: 
              #one with no smoothing (basically a linear regression) with all explanatory variables (elevation, slope, and aspect)
              #one with smoothing on all quantitative explanatory variables (elevation, slope),
              #one with smoothing on the first quantitative explanatory variable (elevation),
              #one with smoothing on the second quantitative explanatory variable (slope)
    # c) We compared the four models AIC values to see which one had the lowest, supporting which model fits the data the best
    # d) We then checked if the conditions of the GAM were met well (normal distribution and equal variance 
        #of the residuals since we are using a GAM with a Gaussian distribution)
    # e) Using the dredge function to determine which explanatory variables allow for the best fitting model
    # f) Comparing the previously selected function and the dredge function to see which one is a better fit 
        #using AIC and an ANOVA F-test (only works if one models variables are entirely present in the other model)
    # g) We then ran a K check on the chosen model (if it used smoothing) to see if the K dimension choices for the model are adequate
        #p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
    # h) We then plotted the GAMs to be able to describe the relationship between the explanatory variables and the size/shape characteristics 
    # i) We then looked for interactions and compared the interaction model to the previous best model and if we have a significant interaction term we have to consider
    # j) We then plotted the interactions both in 2D and 3D

#We used only the 8 categories to get a more specific idea of how direction may be influencing size 

# all points 

#had to remove points 174 and 175 because they had NAs in the slope data and there was a NA in elevation we needed to remove to continue the analysis
all_points_fixed_field_data_processed_terrain_no_NA <- all_points_fixed_field_data_processed_terrain %>%
  filter(is.na(all_points_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(all_points_aspect_raster_15_data_pts_8_categorical) == F)

#Cook's D

# plot(all_points_multiple_lm_SCA)
# all_points_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain_no_NA)
# all_points_mlm_SCA_cooks <- cooks.distance(all_points_mlm_SCA) #calculating the cook.s D for each point
# plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
# influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (2 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
# influential

# removing the NAs from the explanatory and response variables to avoid issues while making the GAMs
LM_fixed_field_data_processed_terrain_dist_no_NA <- LM_fixed_field_data_processed_terrain_dist %>% 
  filter(!is.na(d)) %>% #distance NAs removed
  filter(!is.na(Elevation..m.FIXED)) %>% #Elevation NAs removed
  filter(!is.na(LM_slope_raster_15_data_pts)) %>%  #slope NAs removed
  filter(!is.na(LM_aspect_raster_15_data_pts_8_categorical)) %>% #aspect NAs removed
  filter(!is.na(Canopy_short)) %>% #short canopy axis NAs removed
  filter(!is.na(Canopy_long)) %>% #long canopy axis NAs removed
  filter(!is.na(Canopy_area)) %>% #canopy area NAs removed
  filter(!is.na(Crown_spread)) %>% #crown spread NAs removed
  filter(!is.na(DBH_ag)) #DBH NAs removed


## SCA ##

#removing the spatial geometry to be able to use the GAM function
LM_fixed_field_data_processed_terrain_dist_no_NA <- st_drop_geometry(LM_fixed_field_data_processed_terrain_dist_no_NA)

# Checking a GAM with smoothing splines s(), note we cannot put splines on a categorical variable
LM_add.gam_SCA.terrain_dist <- gam(Canopy_short ~ s(d) + s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                   data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist) #looking at which variables are significant in the linear vs. non-linear model based on the p-values

#elevation has significant non-linear function 

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
# the dredge selected s(d), s(elevation), s(slope) as allowing for the best model

# Checking a GAM with just elevation smoothing splines s()
LM_add.gam_SCA.terrain_dist.just.elevation.smooth <- gam(Canopy_short ~ d + s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                         data = LM_fixed_field_data_processed_terrain_dist_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
summary(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #looking at which variables are significant in the linear vs. non-linear model based on the p-values
#none of the linear fits are significant

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.terrain_dist.just.elevation.smooth) #using the dredge model to narrow the models down to the best choice
dredge[1:5,] #looking at the top five best models, lowest AIC is the best model, rule of thumb is that a difference of 2 is a significant difference
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
dredge[1:5,] #looking at the top five best models, lowest AIC is the best best, rule of thumb is that a difference of 2 is a significant difference
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

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.terrain_dist.dredge) #pretty normal residuals and no heterodescadisticty 

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



#SCA

#creating the basic GAM models

#regular linear regression
all_points_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na.fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
all_points_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
all_points_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(all_points_add.gam_SCA, all_points_add.gam_SCA.smoothed, all_points_add.gam_SCA.smoothed_first_term, 
    all_points_add.gam_SCA.smoothed_second_term)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(all_points_add.gam_SCA)
summary(all_points_add.gam_SCA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(all_points_add.gam_SCA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(all_points_add.gam_SCA.smoothed.dredge, all_points_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(all_points_add.gam_SCA.smoothed.dredge, all_points_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_SCA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#ow p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(all_points_add.gam_SCA.smoothed.dredge)
k.check(all_points_add.gam_SCA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
 plot.gam(all_points_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_SCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
             gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
all_points_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
all_points_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_SCA.smoothed.inter)
summary(all_points_add.gam_SCA.smoothed)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_SCA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)


# LCA

#creating the basic GAM models

#regular linear regression
all_points_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
all_points_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
all_points_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
all_points_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(all_points_add.gam_LCA, all_points_add.gam_LCA.smoothed, all_points_add.gam_LCA.smoothed_first_term, 
    all_points_add.gam_LCA.smoothed_second_term)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(all_points_add.gam_LCA)
summary(all_points_add.gam_LCA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(all_points_add.gam_LCA.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 
#the full model is the dredge output

#Chosen model: all_points_add.gam_LCA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(all_points_add.gam_LCA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_LCA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
all_points_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
all_points_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_LCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_SCA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

# CA

#we choose to log Canopy Area because the residuals were highly non-normal and it chances elevations significance

#creating the basic GAM models

#regular linear regression
all_points_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                              data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
all_points_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
all_points_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
all_points_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                   data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lowered the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(all_points_add.gam_CA, all_points_add.gam_CA.smoothed, all_points_add.gam_CA.smoothed_first_term, 
    all_points_add.gam_CA.smoothed_second_term)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(all_points_add.gam_CA)
summary(all_points_add.gam_CA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(all_points_add.gam_CA.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                               data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(all_points_add.gam_CA.smoothed.dredge, all_points_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(all_points_add.gam_CA.smoothed.dredge, all_points_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_CA.smoothed.dredge
summary(all_points_add.gam_CA.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(all_points_add.gam_CA.smoothed.dredge)
k.check(all_points_add.gam_CA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_CA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
all_points_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                       data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
all_points_add.gam_CA.smoothed.inter <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                             data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_CA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_CA.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

# CS

#creating the basic GAM models

#regular linear regression
all_points_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                             data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
all_points_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
all_points_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                 data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
all_points_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(all_points_add.gam_CS, all_points_add.gam_CS.smoothed, all_points_add.gam_CS.smoothed_first_term, 
    all_points_add.gam_CS.smoothed_second_term)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(all_points_add.gam_CS)
summary(all_points_add.gam_CS.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(all_points_add.gam_CS.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                              data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(all_points_add.gam_CS.smoothed.dredge, all_points_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(all_points_add.gam_CS.smoothed.dredge, all_points_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_CA.smoothed.dredge
summary(all_points_add.gam_CS.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(all_points_add.gam_CS.smoothed.dredge)
k.check(all_points_add.gam_CS.smoothed)


#plotting the chosen function, with no interaction
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_CA.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
all_points_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
all_points_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                            data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_CS.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_CS.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

# DBH_ag

#creating the basic GAM models

#regular linear regression
all_points_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                             data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
all_points_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
all_points_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + all_points_slope_raster_15_data_pts + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                 data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
all_points_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                                  data = all_points_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(all_points_add.gam_DBH, all_points_add.gam_DBH.smoothed, all_points_add.gam_DBH.smoothed_first_term, 
    all_points_add.gam_DBH.smoothed_second_term)
anova(all_points_add.gam_DBH, all_points_add.gam_DBH.smoothed_first_term, 
    all_points_add.gam_DBH.smoothed_second_term, all_points_add.gam_DBH.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(all_points_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(all_points_add.gam_DBH)
summary(all_points_add.gam_DBH.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(all_points_add.gam_DBH.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
all_points_add.gam_DBH.smoothed.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts), 
                                              data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(all_points_add.gam_DBH.smoothed.dredge, all_points_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(all_points_add.gam_DBH.smoothed.dredge, all_points_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: all_points_add.gam_DBH.smoothd.dredge
summary(all_points_add.gam_DBH.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(all_points_add.gam_DBH.smoothed.dredge)
k.check(all_points_add.gam_DBH.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(all_points_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(all_points_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(all_points_add.gam_DBH.smoothed, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
all_points_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                      data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
all_points_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, all_points_slope_raster_15_data_pts) + all_points_aspect_raster_15_data_pts_8_categorical, 
                                            data = all_points_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(all_points_add.gam_DBH.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(all_points_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(all_points_add.gam_DBH.smoothed.inter, "all_points_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=all_points_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=all_points_fixed_field_data_processed_terrain_no_NA$all_points_slope_raster_15_data_pts, 
        z=all_points_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=all_points_fixed_field_data_processed_terrain_no_NA$all_points_aspect_raster_15_data_pts_8_categorical)

# LM

#removing the NAs from slopes and elevation
LM_fixed_field_data_processed_terrain_no_NA <- LM_fixed_field_data_processed_terrain %>%
  filter(is.na(LM_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F)

#running a Cook's D to determine which points are influential/outliers
LM_mlm_SCA <- lm(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain_no_NA)
LM_mlm_SCA_cooks <- cooks.distance(LM_mlm_SCA) #calculating the cook.s D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which cook's D are unsually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with cooks D that are bigger than 3 times the mean cook's D
influential

#removing outliers based on which points were deemed influential using the Cook's D
#in this case, we decided to keep the outliers to get a more accurate view of what is occuring in the population
#LM_fixed_field_data_processed_terrain_no_NA_No_outliers <- LM_fixed_field_data_processed_terrain_no_NA[-c(24,26,27),]


# SCA

#creating the basic GAM models

#regular linear regression
LM_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LM_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LM_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LM_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)


#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LM_add.gam_SCA, LM_add.gam_SCA.smoothed_first_term, 
    LM_add.gam_SCA.smoothed_second_term, LM_add.gam_SCA.smoothed)
anova(LM_add.gam_SCA, LM_add.gam_SCA.smoothed_first_term, 
      LM_add.gam_SCA.smoothed_second_term, LM_add.gam_SCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal variance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LM_add.gam_SCA)
summary(LM_add.gam_SCA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_SCA.smoothed) #using the dredge model to narrowww the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED),
                                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_SCA.smoothed.dredge

summary(LM_add.gam_SCA.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LM_add.gam_SCA.smoothed.dredge)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_SCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LM_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
LM_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                             data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_SCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_SCA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


# LCA

#creating the basic GAM models

#regular linear regression
LM_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LM_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LM_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LM_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LM_add.gam_LCA, LM_add.gam_LCA.smoothed, LM_add.gam_LCA.smoothed_first_term, 
    LM_add.gam_LCA.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LM_add.gam_LCA, LM_add.gam_LCA.smoothed_second_term, 
      LM_add.gam_LCA.smoothed_first_term, LM_add.gam_LCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LM_add.gam_LCA)
summary(LM_add.gam_LCA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_LCA.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 
#the full model is the dredge output

#fitting the dredged model
LM_add.gam_SCA.smoothed.dredge <-  gam(Canopy_long ~ s(Elevation..m.FIXED),
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_LCA.smoothed)
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LM_add.gam_SCA.smoothed.dredge, LM_add.gam_LCA.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_SCA.smoothed.dredge

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LM_add.gam_SCA.smoothed.dredge)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_LCA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LM_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
LM_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_LCA.smoothed.inter)
#there is a not a significant interaction term

#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_LCA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


# CA

#creating the basic GAM models

#regular linear regression
LM_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                             data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LM_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                      data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LM_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LM_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LM_add.gam_CA, LM_add.gam_CA.smoothed, LM_add.gam_CA.smoothed_first_term, 
    LM_add.gam_CA.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LM_add.gam_CA, LM_add.gam_CA.smoothed_first_term, 
    LM_add.gam_CA.smoothed_second_term, LM_add.gam_CA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LM_add.gam_CA)
summary(LM_add.gam_CA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CA.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED), 
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LM_add.gam_CA.smoothed.dredge, LM_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LM_add.gam_CA.smoothed.dredge, LM_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_CA.smoothed.dredge
summary(LM_add.gam_CA.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LM_add.gam_CA.smoothed.dredge)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_CA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LM_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                               data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
LM_add.gam_CA.smoothed.inter <- gam(Canopy_area ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                     data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_CA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_CA.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


# CS

#creating the basic GAM models

#regular linear regression
LM_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                             data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LM_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                      data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LM_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LM_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LM_add.gam_CS, LM_add.gam_CS.smoothed, LM_add.gam_CS.smoothed_first_term, 
    LM_add.gam_CS.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LM_add.gam_CS, LM_add.gam_CS.smoothed_first_term, 
    LM_add.gam_CS.smoothed_second_term, LM_add.gam_CS.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LM_add.gam_CS)
summary(LM_add.gam_CS.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_CS.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED),
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LM_add.gam_CS.smoothed.dredge, LM_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LM_add.gam_CS.smoothed.dredge, LM_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_CA.smoothed.dredge
summary(LM_add.gam_CS.smoothed.dredge)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LM_add.gam_CS.smoothed.dredge)
k.check(LM_add.gam_CS.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_CA.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LM_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
LM_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                    data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_CS.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_CS.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


# DBH_ag

#creating the basic GAM models

#regular linear regression
LM_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LM_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                       data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LM_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LM_slope_raster_15_data_pts + LM_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LM_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LM_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LM_add.gam_DBH, LM_add.gam_DBH.smoothed, LM_add.gam_DBH.smoothed_first_term, 
    LM_add.gam_DBH.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LM_add.gam_DBH, LM_add.gam_DBH.smoothed_first_term, 
    LM_add.gam_DBH.smoothed_second_term, LM_add.gam_DBH.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LM_add.gam_DBH.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LM_add.gam_DBH)
summary(LM_add.gam_DBH.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LM_add.gam_DBH.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
LM_add.gam_DBH.smoothd.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts), 
                                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LM_add.gam_DBH.smoothd.dredge, LM_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LM_add.gam_DBH.smoothd.dredge, LM_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: LM_add.gam_DBH.smoothed
summary(LM_add.gam_DBH.smoothed)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LM_add.gam_DBH.smoothed.dredge)
k.check(LM_add.gam_DBH.smoothed)

#plotting the chosen function, with no interaction
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LM_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LM_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LM_add.gam_DBH.smoothed, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LM_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                              data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#interaction model
LM_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, LM_slope_raster_15_data_pts) + LM_aspect_raster_15_data_pts_8_categorical, 
                                    data = LM_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
summary(LM_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LM_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LM_add.gam_DBH.smoothed.inter, "LM_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LM_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LM_fixed_field_data_processed_terrain_no_NA$LM_slope_raster_15_data_pts, 
        z=LM_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=LM_fixed_field_data_processed_terrain_no_NA$LM_aspect_raster_15_data_pts_8_categorical)


# LC

#removing the NAs from slope, elevation, aspect, and short canopy aspect
LC_fixed_field_data_processed_terrain_no_NA <- LC_fixed_field_data_processed_terrain %>%
  filter(is.na(LC_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(LC_aspect_raster_15_data_pts_8_categorical) == F) %>%
  filter(is.na(Canopy_short) == F)

# SCA


#creating the basic GAM models
  
#regular linear regression  
LC_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                      data = LC_fixed_field_data_processed_terrain_no_NA, na.action = "na.omit") #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LC_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                          data = LC_fixed_field_data_processed_terrain_no_NA)
#smoothing the second explanatory variable (slope)
LC_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                           data = LC_fixed_field_data_processed_terrain_no_NA)

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LC_add.gam_SCA, LC_add.gam_SCA.smoothed_first_term, 
    LC_add.gam_SCA.smoothed_second_term, LC_add.gam_SCA.smoothed)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LC_add.gam_SCA, LC_add.gam_SCA.smoothed_first_term, 
      LC_add.gam_SCA.smoothed_second_term, LC_add.gam_SCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LC_add.gam_SCA)
summary(LC_add.gam_SCA.smoothed)
summary(LC_add.gam_SCA.smoothed.dredge)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_SCA.smoothed) #using the dredge model to narroww the models down to the best choice
dredge[1,] 

#fitting the dredged model
LC_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                               data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LC_add.gam_SCA.smoothed.dredge, LC_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LC_add.gam_SCA.smoothed.dredge, LC_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: LC_add.gam_SCA.smoothed.dredge

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LC_add.gam_SCA.smoothed)

#plotting the chosen function, with no interaction
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_SCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

#looking for interaction
LC_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_SCA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_SCA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


# LCA

#creating the basic GAM models

#regular linear regression
LC_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LC_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
#smoothing the second explanatory variable (slope)
LC_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)


#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LC_add.gam_LCA, LC_add.gam_LCA.smoothed, LC_add.gam_LCA.smoothed_first_term, 
    LC_add.gam_LCA.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LC_add.gam_LCA, LC_add.gam_LCA.smoothed_first_term, 
    LC_add.gam_LCA.smoothed_second_term, LC_add.gam_LCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LC_add.gam_LCA)
summary(LC_add.gam_LCA.smoothed)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_LCA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
LC_add.gam_LCA.smoothed.dredge <-  gam(Canopy_long ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LC_add.gam_LCA.smoothed.dredge, LC_add.gam_LCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LC_add.gam_LCA.smoothed.dredge, LC_add.gam_LCA.smoothed) 
#results show marginal differences

#Chosen model: LC_add.gam_LCA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LC_add.gam_LCA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_LCA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LC_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
#interaction model
LC_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_LCA.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_LCA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)

# CA

#creating the basic GAM models

#regular linear regression
LC_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                             data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LC_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LC_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LC_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LC_add.gam_CA, LC_add.gam_CA.smoothed, LC_add.gam_CA.smoothed_first_term, 
    LC_add.gam_CA.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LC_add.gam_CA, LC_add.gam_CA.smoothed_first_term, 
    LC_add.gam_CA.smoothed_second_term, LC_add.gam_CA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LC_add.gam_CA)
summary(LC_add.gam_CA.smoothed)
summary(LC_add.gam_CA.smoothed.slope.less)

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
LC_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LC_add.gam_CA.smoothed.dredge, LC_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LC_add.gam_CA.smoothed.dredge, LC_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: LC_add.gam_CA.smoothed.dredge

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LC_add.gam_CA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_CA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LC_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
LC_add.gam_CA.smoothed.inter <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                    data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_CA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_CA.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)


# CS

#creating the basic GAM models

#regular linear regression
LC_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                             data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the first explanatory variable (elevation)
LC_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                 data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
LC_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LC_add.gam_CS, LC_add.gam_CS.smoothed, LC_add.gam_CS.smoothed_first_term, 
    LC_add.gam_CS.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LC_add.gam_CS, LC_add.gam_CS.smoothed_first_term, 
    LC_add.gam_CS.smoothed_second_term, LC_add.gam_CS.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(LC_add.gam_CS.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(LC_add.gam_CS)
summary(LC_add.gam_CS.smoothed)
summary(LC_add.gam_CS.smoothed_first_term)
#we can see that elevation and slope does not appear to be as useful

#using the dredge function to determine which explanatory variables allow for the best fitting model
dredge <- dredge(LC_add.gam_CS.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
LC_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LC_add.gam_CS.smoothed.dredge, LC_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LC_add.gam_CS.smoothed.dredge, LC_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: LC_add.gam_CS.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LC_add.gam_CS.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_CA.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LC_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA)
#interaction model
LC_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                    data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_CS.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_CS.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)

# DBH_ag

#creating the basic GAM models

#regular linear regression
LC_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                              data = LC_fixed_field_data_processed_terrain_no_NA) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                       data = LC_fixed_field_data_processed_terrain_no_NA)
#smoothing the first explanatory variable (elevation)
LC_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + LC_slope_raster_15_data_pts + LC_aspect_raster_15_data_pts_8_categorical, 
                                                  data = LC_fixed_field_data_processed_terrain_no_NA)
#smoothing the second explanatory variable (slope)
LC_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                                   data = LC_fixed_field_data_processed_terrain_no_NA)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(LC_add.gam_DBH, LC_add.gam_DBH.smoothed, LC_add.gam_DBH.smoothed_first_term, 
    LC_add.gam_DBH.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(LC_add.gam_DBH, LC_add.gam_DBH.smoothed_first_term, 
    LC_add.gam_DBH.smoothed_second_term, LC_add.gam_DBH.smoothed)

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(LC_add.gam_DBH.smoothed)
#comparing the models with an ANOVA test to further provide evidence for the chosen model

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
summary(LC_add.gam_DBH)
summary(LC_add.gam_DBH.smoothed)
summary(LC_add.gam_DBH.smoothed_first_term)
#based on these results we can see that the normality condition is not well met, so we can try

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(LC_add.gam_DBH.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
LC_add.gam_DBH.smoothed.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + LC_aspect_raster_15_data_pts_8_categorical, 
                                      data = LC_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(LC_add.gam_DBH.smoothed.dredge, LC_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(LC_add.gam_DBH.smoothed.dredge, LC_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: LC_add.gam_DBH.smoothed.dredge
summary(LC_add.gam_DBH.smoothed)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(LC_add.gam_DBH.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(LC_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(LC_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(LC_add.gam_DBH.smoothed, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
LC_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                               data = LC_fixed_field_data_processed_terrain_no_NA)
#interaction model
LC_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, LC_slope_raster_15_data_pts) + LC_aspect_raster_15_data_pts_8_categorical, 
                                     data = LC_fixed_field_data_processed_terrain_no_NA)
summary(LC_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(LC_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(LC_add.gam_DBH.smoothed.inter, "LC_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=LC_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=LC_fixed_field_data_processed_terrain_no_NA$LC_slope_raster_15_data_pts, 
        z=LC_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=LC_fixed_field_data_processed_terrain_no_NA$LC_aspect_raster_15_data_pts_8_categorical)

# SD

#removing the NAs from slope, elevation, and aspect
SD_fixed_field_data_processed_terrain_no_NA <- SD_fixed_field_data_processed_terrain %>%
  filter(is.na(SD_slope_raster_15_data_pts) == F) %>%
  filter(is.na(Elevation..m.FIXED) == F) %>%
  filter(is.na(SD_aspect_raster_15_data_pts_8_categorical) == F)

# SCA

#creating the basic GAM models

#regular linear regression
SD_add.gam_SCA <- gam(Canopy_short ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                      data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the first explanatory variable (elevation)
SD_add.gam_SCA.smoothed_first_term <- gam(Canopy_short ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                          data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
SD_add.gam_SCA.smoothed_second_term <- gam(Canopy_short ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                           data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(SD_add.gam_SCA, SD_add.gam_SCA.smoothed_first_term, 
    SD_add.gam_SCA.smoothed_second_term, SD_add.gam_SCA.smoothed)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(SD_add.gam_SCA, SD_add.gam_SCA.smoothed_first_term, 
      SD_add.gam_SCA.smoothed_second_term, SD_add.gam_SCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_SCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(SD_add.gam_SCA)
summary(SD_add.gam_SCA.smoothed)

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(SD_add.gam_SCA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
SD_add.gam_SCA.smoothed.dredge <-  gam(Canopy_short ~ s(Elevation..m.FIXED) + s(Elevation..m.FIXED), 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(SD_add.gam_SCA.smoothed.dredge, SD_add.gam_SCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(SD_add.gam_SCA.smoothed.dredge, SD_add.gam_SCA.smoothed) 
#results show marginal differences

#Chosen model: SD_add.gam_SCA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(SD_add.gam_SCA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_SCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_SCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_SCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

#looking for interaction
SD_add.gam_SCA.smoothed <- gam(Canopy_short ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
SD_add.gam_SCA.smoothed.inter <- gam(Canopy_short ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_SCA.smoothed.inter)
#there is a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_SCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_SCA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Short Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_short, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

# LCA

#creating the basic GAM models

#regular linear regression
SD_add.gam_LCA <- gam(Canopy_long ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
SD_add.gam_LCA.smoothed_first_term <- gam(Canopy_long ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
SD_add.gam_LCA.smoothed_second_term <- gam(Canopy_long ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(SD_add.gam_LCA, SD_add.gam_LCA.smoothed, SD_add.gam_LCA.smoothed_first_term, 
    SD_add.gam_LCA.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(SD_add.gam_LCA, SD_add.gam_LCA.smoothed_first_term, 
    SD_add.gam_LCA.smoothed_second_term, SD_add.gam_LCA.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_LCA.smoothed)
#based on these results we can see that the normality condition is not well met, so we can try

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(SD_add.gam_LCA)
summary(SD_add.gam_LCA.smoothed)

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(SD_add.gam_LCA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
SD_add.gam_LCA.smoothed.dredge <-  gam(Canopy_long ~ s(Elevation..m.FIXED) + s(Elevation..m.FIXED), 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(SD_add.gam_LCA.smoothed.dredge, SD_add.gam_LCA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(SD_add.gam_LCA.smoothed.dredge, SD_add.gam_LCA.smoothed) 
#results show marginal differences

#Chosen model: SD_add.gam_LCA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(SD_add.gam_LCA.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_LCA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_LCA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_LCA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
SD_add.gam_LCA.smoothed <- gam(Canopy_long ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
#interaction model
SD_add.gam_LCA.smoothed.inter <- gam(Canopy_long ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_LCA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_LCA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_LCA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Long Canopy Axis")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_long, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

# CA

#creating the basic GAM models

#regular linear regression
SD_add.gam_CA <- gam(log(Canopy_area) ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                             data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
SD_add.gam_CA.smoothed <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
SD_add.gam_CA.smoothed_first_term <- gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                 data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
SD_add.gam_CA.smoothed_second_term <- gam(log(Canopy_area) ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(SD_add.gam_CA, SD_add.gam_CA.smoothed, SD_add.gam_CA.smoothed_first_term, 
    SD_add.gam_CA.smoothed_second_term)

summary(SD_add.gam_CA.smoothed)
#we can see that elevation does not appear to be as useful

#checking overall fit and potential issues
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CA.smoothed)

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(SD_add.gam_CA)
summary(SD_add.gam_CA.smoothed)

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(SD_add.gam_CA.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
SD_add.gam_CA.smoothed.dredge <-  gam(log(Canopy_area) ~ s(Elevation..m.FIXED) + s(Elevation..m.FIXED), 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(SD_add.gam_CA.smoothed.dredge, SD_add.gam_CA.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(SD_add.gam_CA.smoothed.dredge, SD_add.gam_CA.smoothed) 
#results show marginal differences

#Chosen model: SD_add.gam_CA.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(SD_add.gam_CA.smoothed)

#plotting the chosen function, with no interaction
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_CA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
SD_add.gam_CA.smoothed <- gam(Canopy_area ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
#interaction model
SD_add.gam_CA.smoothed.inter <- gam(Canopy_area ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                    data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_CA.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_CA.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_CA.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Canopy_area, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


# CS


#creating the basic GAM models

#regular linear regression
SD_add.gam_CS <- gam(Crown_spread ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                             data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
SD_add.gam_CS.smoothed_first_term <- gam(Crown_spread ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                 data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
SD_add.gam_CS.smoothed_second_term <- gam(Crown_spread ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(SD_add.gam_CS, SD_add.gam_CS.smoothed, SD_add.gam_CS.smoothed_first_term, 
    SD_add.gam_CS.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(SD_add.gam_CS, SD_add.gam_CS.smoothed_first_term, 
    SD_add.gam_CS.smoothed_second_term, SD_add.gam_CS.smoothed)

summary(SD_add.gam_CS.smoothed)
#we can see that elevation does not appear to be as useful

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_CS.smoothed)

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(SD_add.gam_CS)
summary(SD_add.gam_CS.smoothed)

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(SD_add.gam_CS.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
SD_add.gam_CS.smoothed.dredge <-  gam(Crown_spread ~ s(Elevation..m.FIXED) + s(Elevation..m.FIXED), 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(SD_add.gam_CS.smoothed.dredge, SD_add.gam_CS.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(SD_add.gam_CS.smoothed.dredge, SD_add.gam_CS.smoothed) 
#results show marginal differences

#Chosen model: SD_add.gam_CS.smoothed

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(SD_add.gam_CS.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_CA.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_CA.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_CA.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Canopy Area")  # Uses ggplot2 for a cleaner plot


#checking for significant interaction terms

#chosen function
SD_add.gam_CS.smoothed <- gam(Crown_spread ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA)
#interaction model
SD_add.gam_CS.smoothed.inter <- gam(Crown_spread ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                    data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_CS.smoothed.inter)
#there is not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_CS.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_CS.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on Crown Spread")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$Crown_spread, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)


# DBH_ag

#creating the basic GAM models

#regular linear regression
SD_add.gam_DBH <- gam(DBH_ag ~ Elevation..m.FIXED + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                              data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail) #na fail makes sure the later dredge does not have to worry about NAs
#smoothing both quantitative variables
SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                       data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)
#smoothing the first explanatory variable (elevation)
SD_add.gam_DBH.smoothed_first_term <- gam(DBH_ag ~ s(Elevation..m.FIXED) + SD_slope_raster_15_data_pts + SD_aspect_raster_15_data_pts_8_categorical, 
                                                  data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#smoothing the second explanatory variable (slope)
SD_add.gam_DBH.smoothed_second_term <- gam(DBH_ag ~ Elevation..m.FIXED + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                                   data = SD_fixed_field_data_processed_terrain_no_NA, na.action = na.fail)
#logging canopy area lower the AIC significantly

#comparing the models' AIC values, the model with the lowest value shows the best fit
AIC(SD_add.gam_DBH, SD_add.gam_DBH.smoothed, SD_add.gam_DBH.smoothed_first_term, 
    SD_add.gam_DBH.smoothed_second_term)
#comparing the models with an ANOVA test to further provide evidence for the chosen model
anova(SD_add.gam_DBH, SD_add.gam_DBH.smoothed_first_term, 
    SD_add.gam_DBH.smoothed_second_term, SD_add.gam_DBH.smoothed)

#checking conditions for our GAM which assumes a Gaussian distributed (normal distribution and equal vairance of residuals assumption)
par(mfrow = c(2, 2))
gam.check(SD_add.gam_DBH.smoothed)
gam.check(SD_add.gam_DBH.smoothed.dredge)

#comparing the model's the models GCV summary values to see which is lowest as another method of comparing the fit of the models
summary(SD_add.gam_DBH)
summary(SD_add.gam_DBH.smoothed)
summary(SD_add.gam_DBH.smoothed.dredge)

#using the dredge function to determine which explanatory variables allows for the best fitting model
dredge <- dredge(SD_add.gam_DBH.smoothed) #using the dredge model to narrow the models down to the best choice
dredge[1,] 

#fitting the dredged model
SD_add.gam_DBH.smoothed.dredge <-  gam(DBH_ag ~ s(Elevation..m.FIXED) + s(Elevation..m.FIXED), 
                                      data = SD_fixed_field_data_processed_terrain_no_NA,  na.action = na.fail)

#Comparing the dredge model and the previously-considered best model

#Anova F-test comparing strength of dredge vs. full model demonstrates dredge performs just as well 
#significance means the model with more variables explains significantly more of the response
anova(SD_add.gam_DBH.smoothed.dredge, SD_add.gam_DBH.smoothed, test = "F")
#AIC comparing the dredge and full model to see which one is a better fit to the data
AIC(SD_add.gam_DBH.smoothed.dredge, SD_add.gam_DBH.smoothed) 
#results show marginal differences

#Chosen model: SD_add.gam_DBH.smoothed
summary(SD_add.gam_DBH.smoothed)

#checking K to see if our smoothing terms are K dimension choices for the model are adequate
#p-values may indicate that the basis dimension, k, has been set too low, especially if the reported edf is close to k
k.check(SD_add.gam_DBH.smoothed)

#plotting the chosen function, with no interaction 
par(mfrow = c(3,2), mar = c(4.5, 4.5, 1, 1))
plot.gam(SD_add.gam_DBH.smoothed, select=1, 
         all.terms=T, xlab = "Elevation (m)", 
         ylab = expression(f[1]*'(Elevation)'), se = TRUE , col = "black")
plot.gam(SD_add.gam_DBH.smoothed, select=2, 
         all.terms=T, xlab = "Slope (º)", ylab = "f_1 (Slope)")
visreg(SD_add.gam_DBH.smoothed, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

#checking for significant interaction terms

#chosen function
SD_add.gam_DBH.smoothed <- gam(DBH_ag ~ s(Elevation..m.FIXED) + s(SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                               data = SD_fixed_field_data_processed_terrain_no_NA)
#interaction model
SD_add.gam_DBH.smoothed.inter <- gam(DBH_ag ~ s(Elevation..m.FIXED, SD_slope_raster_15_data_pts) + SD_aspect_raster_15_data_pts_8_categorical, 
                                     data = SD_fixed_field_data_processed_terrain_no_NA)
summary(SD_add.gam_DBH.smoothed.inter)
#there is a not a significant interaction term

#interaction plots
par(mfrow = c(2,2), mar = c(4.5, 4.5, 2, 2))
plot.gam(SD_add.gam_DBH.smoothed.inter, select=1, 
         all.terms=T, xlab = "s(Elevation (m):Slope (º))", main = "s(Elevation:Slope)", 
         ylab = expression(f[1]*'(Elevation (m):Slope (º))'), se = TRUE,
         cex.axis = 1, cex.main = 1, cex.lab = 1)
legend("topright", col = c("lightgreen", "black", "#F08080"), lty = c(3, 1, 2), legend = c("+1 SE", "Fit", "-1 SE"))
visreg(SD_add.gam_DBH.smoothed.inter, "SD_aspect_raster_15_data_pts_8_categorical",
       gg = F, xlab = "Aspect", ylab = "Effect on DBH")  # Uses ggplot2 for a cleaner plot

# 3d plotting in plotly and with gg3D
plot_ly(x=SD_fixed_field_data_processed_terrain_no_NA$Elevation..m.FIXED, 
        y=SD_fixed_field_data_processed_terrain_no_NA$SD_slope_raster_15_data_pts, 
        z=SD_fixed_field_data_processed_terrain_no_NA$DBH_ag, type="scatter3d", mode="markers", 
        color=SD_fixed_field_data_processed_terrain_no_NA$SD_aspect_raster_15_data_pts_8_categorical)

