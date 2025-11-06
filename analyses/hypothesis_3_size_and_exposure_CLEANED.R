# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei's size/shape is influenced by either elevation, slope, and/or aspect%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#the purpose of the script is to determine how the elevation, slope, and aspect each 
#relate to size/shape of the trees using single linear regressions. 

#The script is broken up into these sections:

  # 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data includes:
          # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
          # Extracting and processing slope, elevation, and aspect (4 and 8 cardinal directions) data using 15 m res rasters,
          # Extracting the distance to the river of each tree for each population,
    # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
  # 2) Comparing size of the trees to their elevation to look for a relationship using single linear regression
  # 3) Comparing size of the trees to their slope to look for a relationship using single linear regression
  # 4) Comparing size of the trees to their aspect to look for a relationship using ANOVA / Kruskal-Wallis Tests

# NOTE: Uncomment and run line 40, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

#### Loading libraries and relevant data ####

library(googledrive) #to download files from google drive
library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(raster) #for working with the rast files
library(terra) # for extracting the slope and aspect from the DEM elevation files
library(perm.t.test) #permutation t test 
library(car) #to create added variable plots and to run levene's test for checking ANOVA conditions
library(stars) # to convert raster into stars
library(gdalUtilities) #to be able to use gdalwarp

# loading in the processed tree data 
# NOTE: Uncomment and run line 40, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")

#### Sizes vs. Elevation ####

# For all populations/each population and size/shape metric we created single variable regressions by...
  #a) creating the single variable linear regressions with the un-transformed response variable, logged variable, and 
          # square root of the variable, respectively, either with or without outliers
  #b) testing which model best satisfies the conditions for the analysis: LINES
          # Linearity, Independence, Normality of residuals, Equal variance of residuals, and simple random sample 
                # we tested Linearity by looking at the scatterplots,
                # we tested Independence by thinking about the explanatory and response variables across the points,
                # we tested Normality of Residuals using histograms, qq norm plots, and the Shapiro-Wilk's test,
                # we tested Equal Variance of Residuals using a fitted vs. residuals plot,
                # we tested Simple Random Sample by thinkg about the data collection method.
  #c) We then looked for significant associations (slopes/correlations)
        # 1) If the LINES conditions are met...
               # we ran a slope test and a Pearson's correlation test to see if there is a significant association
        # 2) If the LINES conditions are not met...
               # we ran a Mann-Kendall test (non-parametric test) to look for a significant correlation/tau 

# For all trees

# removing NAs
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_trans_coordinates %>%
  drop_na(Elevation..m.FIXED) #removing NAs in elevation 

#SCA

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_SCA <- lm(Canopy_short ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates) #creating a linear regression to use to calculate the Cook's D
all_points_slr_SCA_cooks <- cooks.distance(all_points_slr_SCA) #calculating the Cook's D for each point
plot(all_points_slr_SCA_cooks, type = 'h') #checking to see which cook's D are unusually high
influential <- all_points_slr_SCA_cooks[(all_points_slr_SCA_cooks > (3 * mean(all_points_slr_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(as.numeric(names(influential))),]


#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#creating the linear regression without transformations
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_short_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Canopy_long_lg ~ fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_sca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Canopy_long_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_sca_no_outliers$Elevation..m.FIXED)


# square root transformation does the best job of meeting the conditions

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_sca_elev, aes(x= all_points_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_sca_elev, aes(sample = all_points_lm_sca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_sca_elev$residuals) #only not significant for a square root transformation, we could use mann-kendall test for non-parametric data or the square root transformation

#checking equal variance with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = all_points_lm_sca_elev, aes(x = all_points_lm_sca_elev$fitted.values, y = all_points_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test, obtained from the summary of the linear regression reslts
summary(all_points_lm_sca_elev)

#correlation test
cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates) #creating a linear regression to use to calculate the Cook's D
all_points_slr_LCA_cooks <- cooks.distance(all_points_slr_LCA) #calculating the Cook's D for each point
plot(all_points_slr_LCA_cooks, type = 'h') #checking to see which cook's D are unusually high
influential <- all_points_slr_LCA_cooks[(all_points_slr_LCA_cooks > (3 * mean(all_points_slr_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 49, 88, 89, 105, 126, 169, 177, 189, 
                                                                                                                      190, 194, 208, 210, 212, 27, 218, 219, 242, 250,
                                                                                                                      252, 254, 258, 270, 290, 291, 295, 304, 305, 306, 307,                                                                                                               
                                                                                                                      480, 494, 514, 572, 643, 648),]
#creating the linear regressions

#linear regression without transformations
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_long_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long_lg ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_lca_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Canopy_long_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_lca_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(all_points_lm_lca_elev, aes(x= all_points_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot plot
ggplot(all_points_lm_lca_elev, aes(sample = all_points_lm_lca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(all_points_lm_lca_elev$residuals) #it is sig despite transformations and removal of outliers, so we have to use Mann-Kendall non-parametric test

#checking equal variance with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = all_points_lm_lca_elev, aes(x = all_points_lm_lca_elev$fitted.values, y = all_points_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_lca_elev)

#non-parametric Mann-Kendall Test, non-parametric test
all_points_tau_result_LCA <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_long,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_LCA)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Canopy_long ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = (fixed_field_data_processed_sf_trans_coordinates$Canopy_long), color = "blue")) +
  geom_line(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = all_points_trend_line_LCA), color = "red") +
  labs(x = "Elevation", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#Canopy Area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates) #creating a linear regression to use to calculate the Cook's D
all_points_slr_CA_cooks <- cooks.distance(all_points_slr_CA) #calculating the Cook's D for each point
plot(all_points_slr_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_CA_cooks[(all_points_slr_CA_cooks > (3 * mean(all_points_slr_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 49, 89, 163, 169, 177,
                                                                                                                     190, 208, 210, 212, 217, 219, 249, 251,
                                                                                                                     252, 253, 254, 270, 304, 305, 306,
                                                                                                                     309, 318, 321, 338, 343, 375, 376, 377, 378, 379, 384, 
                                                                                                                     398, 433, 449, 472, 
                                                                                                                     480, 494, 514, 572, 648, 660),]
                                                                                                                     
#creating the linear regressions

#linear regression without transformations
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Canopy_area_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area_lg ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_CA_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Canopy_area_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_ca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(all_points_lm_CA_elev, aes(x= all_points_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot  
ggplot(all_points_lm_CA_elev, aes(sample = all_points_lm_CA_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(all_points_lm_CA_elev$residuals) #shapiro-wilk test is sig depsite transformations and removal of outliers, so we need to use a non-parametric mann-kendall test

#checking equal variance with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = all_points_lm_CA_elev, aes(x = all_points_lm_CA_elev$fitted.values, y = all_points_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_CA_elev)

#non-parametric Mann-Kendall Test, non-parametric test
all_points_tau_result_CA <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Canopy_area,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_CA)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Canopy_area ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = (fixed_field_data_processed_sf_trans_coordinates$Canopy_area), color = "blue")) +
  geom_line(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = all_points_trend_line_LCA), color = "red") +
  labs(x = "Elevation (m)", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#Crown Spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_CS <- lm(Crown_spread ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates) #creating a linear regression to use to calculate the Cook's D
all_points_slr_CS_cooks <- cooks.distance(all_points_slr_CS) #calculating the Cook's D for each point
plot(all_points_slr_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_CS_cooks[(all_points_slr_CS_cooks > (3 * mean(all_points_slr_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(42, 44, 45, 88, 89, 163, 169, 177,
                                                                                                                    189, 190, 194, 208, 210, 212, 217, 219, 242, 249, 251,
                                                                                                                     252, 253, 254, 270, 290, 291, 295, 304, 305, 306,
                                                                                                                     309, 318, 321, 338, 343, 375, 376, 377, 378, 379, 381, 384, 
                                                                                                                     398, 433, 449, 
                                                                                                                     480, 494, 514, 572, 648, 660),]

#creating the linear regressions

#linear regression without transformations
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$Crown_spread_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread_lg ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_CS_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Crown_spread_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_cs_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(all_points_lm_CS_elev, aes(x= all_points_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_CS_elev, aes(sample = all_points_lm_CS_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_CS_elev$residuals) #all were sig, despite transformations and removal of outliers, so need to use non-parametric mann-kendall test

#checking equal variance with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = all_points_lm_CS_elev, aes(x = all_points_lm_CS_elev$fitted.values, y = all_points_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_CS_elev)

#non-parametric Mann-Kendall Test, non-parametric test
all_points_tau_result_CS <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$Crown_spread,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_CS)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$Crown_spread ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#removing outliers

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_DBH <- lm(DBH_ag ~ Elevation..m.FIXED, data = fixed_field_data_processed_sf_trans_coordinates) #creating a linear regression to use to calculate the Cook's D
all_points_slr_DBH_cooks <- cooks.distance(all_points_slr_DBH) #calculating the Cook's D for each point
plot(all_points_slr_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_DBH_cooks[(all_points_slr_DBH_cooks > (3 * mean(all_points_slr_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers <- fixed_field_data_processed_sf_trans_coordinates[-c(16,49, 53, 65, 73, 79, 89, 94, 96, 125, 152, 160,
                                                                                                                      163, 165, 169, 177, 208, 210, 212, 214, 219, 244, 
                                                                                                                      245, 246, 247, 249, 250, 251,
                                                                                                                      252, 253, 254, 255, 257, 258, 270, 280, 283, 285,
                                                                                                                      286, 290, 304, 305, 306,
                                                                                                                      314, 318, 321, 338, 343, 353, 369, 373, 374, 
                                                                                                                      376, 378, 379, 384, 
                                                                                                                      398, 449, 613, 648, 660),]

#creating the linear regressions

#linear regression without transformations
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_lg ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates$DBH_ag_sqrt ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag_lg ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
all_points_lm_DBH_elev  <- lm(fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$DBH_ag_sqrt ~ fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(all_points_lm_DBH_elev, aes(x= all_points_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_DBH_elev, aes(sample = all_points_lm_DBH_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_DBH_elev$residuals) #all were sig, despite transformations and removal of outliers, so need to use non-parametric mann-kendall test

#plotting the scatterplot and linear model in ggplot
ggplot(data = fixed_field_data_processed_sf_trans_coordinates_dbh_no_outliers, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#checking equal variance with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = all_points_lm_DBH_elev, aes(x = all_points_lm_DBH_elev$fitted.values, y = all_points_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_DBH_elev)

#non-parametric Mann-Kendall Test, non-parametric test
all_points_tau_result_DBH <- cor.test(fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, fixed_field_data_processed_sf_trans_coordinates$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_DBH)

# Calculating the trend line for plotting
all_points_trend_line_DBH <- predict(loess(fixed_field_data_processed_sf_trans_coordinates$DBH_ag ~ fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = (fixed_field_data_processed_sf_trans_coordinates$DBH_ag), color = "blue")) +
  geom_line(aes(x = fixed_field_data_processed_sf_trans_coordinates$Elevation..m.FIXED, y = all_points_trend_line_DBH), color = "red") +
  labs(x = "Elevation", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

# LM 

#removing NAs in elevation
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  drop_na(Elevation..m.FIXED)

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis (m)")+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12))

#using Cook's D to check for highly influential points that may skew the linear model results
LM_slr_SCA <- lm(Canopy_short ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LM_slr_SCA_cooks <- cooks.distance(LM_slr_SCA) #calculating the Cook's D for each point
plot(LM_slr_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_slr_SCA_cooks[(LM_slr_SCA_cooks > (3 * mean(LM_slr_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_sca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 44, 45, 113, 125, 126, 147, 
                                                                                  169, 177, 189, 194, 208, 210, 219),]

#creating the linear regressions

#linear regression without transformations
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed$Canopy_short_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed_sca_no_outliers$Canopy_short ~ LM_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed_sca_no_outliers$Canopy_short_lg ~ LM_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_sca_elev  <- lm(LM_fixed_field_data_processed_sca_no_outliers$Canopy_short_sqrt ~ LM_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LM_lm_sca_elev, aes(x= LM_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_elev, aes(sample = LM_lm_sca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_sca_elev$residuals) #normality met when using outliers and square root of SCA

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_sca_elev, aes(x = LM_lm_sca_elev$fitted.values, y = LM_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_sca_elev)

#correlation test
cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed) #calculating the Cook's D for each point
LM_slr_LCA_cooks <- cooks.distance(LM_slr_LCA) #calculating the Cook's D for each point
plot(LM_slr_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_slr_LCA_cooks[(LM_slr_LCA_cooks > (3 * mean(LM_slr_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_lca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 42, 49, 88, 113, 126, 138, 141, 147, 169, 
                                                                                  177, 189, 185, 189, 190, 194, 208, 210, 212, 217, 218, 219),]


#creating the linear regression

#linear regression without transformations
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed$Canopy_long_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long_lg ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_lca_elev  <- lm(LM_fixed_field_data_processed_lca_no_outliers$Canopy_long_sqrt ~ LM_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LM_lm_lca_elev, aes(x= LM_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_lca_elev, aes(sample = LM_lm_lca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(LM_lm_lca_elev$residuals) #sign not normally distributed with just transformations,but non sign with square root trans and removal of outliers

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_lca_elev, aes(x = LM_lm_lca_elev$fitted.values, y = LM_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_lca_elev)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_LCA <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_long,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(LM_tau_result_LCA)

# Calculating the trend line for plotting
LM_trend_line_LCA <- predict(loess(LM_fixed_field_data_processed$Canopy_long ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LM_slr_CA_cooks <- cooks.distance(LM_slr_CA) #calculating the Cook's D for each point
plot(LM_slr_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_slr_CA_cooks[(LM_slr_CA_cooks > (3 * mean(LM_slr_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_ca_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 44, 89, 113, 125, 126, 147, 169, 177, 
                                                                                 189, 208, 210, 212, 217, 219),]

#creating the linear regression

#linear regression without transformations
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed$Canopy_area_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area_lg ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_CA_elev  <- lm(LM_fixed_field_data_processed_ca_no_outliers$Canopy_area_sqrt ~ LM_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(Canopy Area_")

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LM_lm_CA_elev, aes(x= LM_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_CA_elev, aes(sample = LM_lm_CA_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_CA_elev$residuals) #with just transformations sig, but with removal of outliers and square root trans it is normal

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_CA_elev)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_CA <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_area,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(LM_tau_result_CA)

# Calculating the trend line for plotting
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed$Canopy_area ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_slr_CS <- lm(Crown_spread ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LM_slr_CS_cooks <- cooks.distance(LM_slr_CS) #calculating the Cook's D for each point
plot(LM_slr_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_slr_CS_cooks[(LM_slr_CS_cooks > (3 * mean(LM_slr_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_cs_no_outliers <- LM_fixed_field_data_processed[-c(21, 29, 44, 113, 125, 126, 147, 169, 177, 
                                                                                 189, 194, 208, 210, 212, 217, 219),]



#creating the linear regression

#linear regression without transformations
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed$Crown_spread_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)


#creating the linear regressions without any outliers 

#linear regression without transformations
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed_cs_no_outliers$Crown_spread ~ LM_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed_cs_no_outliers$Crown_spread_lg ~ LM_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_CS_elev  <- lm(LM_fixed_field_data_processed_cs_no_outliers$Crown_spread_sqrt ~ LM_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LM_lm_CS_elev, aes(x= LM_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_CS_elev, aes(sample = LM_lm_CS_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_CS_elev$residuals) #not significant with a square root transformation 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(Crown Spread)")

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_CS_elev, aes(x = LM_lm_CS_elev$fitted.values, y = LM_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_CS_elev)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_CS <- cor.test(LM_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED, LM_fixed_field_data_processed_cs_no_outliers$Crown_spread,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(LM_tau_result_CS)

# Calculating the trend line for plotting
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed_cs_no_outliers$Crown_spread ~ LM_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()


#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_slr_DBH <- lm(DBH_ag ~ Elevation..m.FIXED, data = LM_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LM_slr_DBH_cooks <- cooks.distance(LM_slr_DBH) #calculating the Cook's D for each point
plot(LM_slr_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_slr_DBH_cooks[(LM_slr_DBH_cooks > (3 * mean(LM_slr_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_dbh_no_outliers <- LM_fixed_field_data_processed[-c(16, 21, 29, 53, 73, 79, 94, 96, 114, 125, 126, 152, 160, 162, 165,
                                                                                 169, 177, 208, 210, 212, 214, 219),]





#creating the linear regressions

#linear regression without transformations
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_lg ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed$DBH_ag_sqrt ~ LM_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed_dbh_no_outliers$DBH_ag ~ LM_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed_dbh_no_outliers$DBH_ag_lg ~ LM_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
LM_lm_DBH_elev  <- lm(LM_fixed_field_data_processed_dbh_no_outliers$DBH_ag_sqrt ~ LM_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("log(DBH)")+
  abline(LM_lm_DBH_elev)+
  theme(axis.title.x = element_text(size=15),
        axis.title.y = element_text(size=15),
        axis.text.x = element_text(size=12), 
        axis.text.y = element_text(size=12))

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LM_lm_DBH_elev, aes(x= LM_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_DBH_elev, aes(sample = LM_lm_DBH_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_DBH_elev$residuals) #significant for all transformations, have to ue non-parametric method

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_DBH_elev, aes(x = LM_lm_DBH_elev$fitted.values, y = LM_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_DBH_elev)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_DBH <- cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(LM_tau_result_DBH)

# Calculating the trend line for plotting
LM_trend_line_DBH <- predict(loess(LM_fixed_field_data_processed$DBH_ag ~ LM_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_processed$Elevation..m.FIXED, y = (LM_fixed_field_data_processed$DBH_ag), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_processed$Elevation..m.FIXED, y = LM_trend_line_DBH), color = "red") +
  labs(x = "Elevation", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

## LC linear models ##

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_slr_SCA <- lm(Canopy_short ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LC_slr_SCA_cooks <- cooks.distance(LC_slr_SCA) #calculating the Cook's D for each point
plot(LC_slr_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_slr_SCA_cooks[(LC_slr_SCA_cooks > (3 * mean(LC_slr_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_sca_no_outliers <- LC_fixed_field_data_processed[-c(34, 38, 54, 74, 77, 85, 99, 119,
                                                                                  177, 195, 203),]

#creating the linear regressions

#linear regression without transformations
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed_sca_no_outliers$Canopy_short ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_sca_elev  <- lm(LC_fixed_field_data_processed$Canopy_short_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#plotting the scatterplot and linear model in ggplot with a logged SCA value
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_short_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("log(Short Canopy Axis)")

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LC_lm_sca_elev, aes(x= LC_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_elev, aes(sample = LC_lm_sca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_sca_elev$residuals) #not significant with logged transformation of sca

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_sca_elev, aes(x = LC_lm_sca_elev$fitted.values, y = LC_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_sca_elev)

#correlation test
cor.test(LM_fixed_field_data_processed$Elevation..m.FIXED, LM_fixed_field_data_processed$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Canopy_long_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LC_slr_LCA_cooks <- cooks.distance(LC_slr_LCA) #calculating the Cook's D for each point
plot(LC_slr_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_slr_LCA_cooks[(LC_slr_LCA_cooks > (3 * mean(LC_slr_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_lca_no_outliers <- LC_fixed_field_data_processed[-c(3, 13, 19, 32, 34, 38, 54, 55, 74, 77,
                                                                                  85, 99, 119, 126, 159, 177),]

#creating the linear regressions

#linear regression without transformations
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed$Canopy_long_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed_lca_no_outliers$Canopy_long ~ LC_fixed_field_data_processed_lca_no_outliers$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed_lca_no_outliers$Canopy_long_lg ~ LC_fixed_field_data_processed_lca_no_outliers$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_lca_elev  <- lm(LC_fixed_field_data_processed_lca_no_outliers$Canopy_long_sqrt ~ LC_fixed_field_data_processed_lca_no_outliers$Elevation..m.)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LC_lm_lca_elev, aes(x= LC_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_lca_elev, aes(sample = LC_lm_lca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_lca_elev$residuals) #not significant with logged transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_lca_elev, aes(x = LC_lm_lca_elev$fitted.values, y = LC_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LC_slr_CA_cooks <- cooks.distance(LC_slr_CA) #calculating the Cook's D for each point
plot(LC_slr_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_slr_CA_cooks[(LC_slr_CA_cooks > (3 * mean(LC_slr_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_ca_no_outliers <- LC_fixed_field_data_processed[-c(3, 34, 38, 54, 74, 77, 85, 99, 119, 177, 203),]


#creating the linear regressions

#linear regression without transformations
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed$Canopy_area_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#creating the linear regressions without any outliers

LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed_ca_no_outliers$Canopy_area ~ LC_fixed_field_data_processed_ca_no_outliers$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed_ca_no_outliers$Canopy_area_lg ~ LC_fixed_field_data_processed_ca_no_outliers$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_CA_elev  <- lm(LC_fixed_field_data_processed_ca_no_outliers$Canopy_area_sqrt ~ LC_fixed_field_data_processed_ca_no_outliers$Elevation..m.)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LC_lm_CA_elev, aes(x= LC_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_CA_elev, aes(sample = LC_lm_CA_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_CA_elev$residuals) #not significant with log transformation of CA

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_CA_elev, aes(x = LM_lm_CA_elev$fitted.values, y = LM_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_CA_elev)

#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_slr_CS <- lm(Crown_spread ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LC_slr_CS_cooks <- cooks.distance(LC_slr_CS) #calculating the Cook's D for each point
plot(LC_slr_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_slr_CS_cooks[(LC_slr_CS_cooks > (3 * mean(LC_slr_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_cs_no_outliers <- LC_fixed_field_data_processed[-c(3, 13, 34, 38, 54, 55, 74, 77, 85, 99, 119, 177, 203),]

#creating the linear regressions

#linear regression without transformations
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed$Crown_spread_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed_cs_no_outliers$Crown_spread ~ LC_fixed_field_data_processed_cs_no_outliers$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed_cs_no_outliers$Crown_spread_lg ~ LC_fixed_field_data_processed_cs_no_outliers$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_CS_elev  <- lm(LC_fixed_field_data_processed_cs_no_outliers$Crown_spread_sqrt ~ LC_fixed_field_data_processed_cs_no_outliers$Elevation..m.)

#plotting the scatterplot and linear model in ggplot with logged crown spread
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=Crown_spread_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("log(Crown Spread)")

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LC_lm_CS_elev, aes(x= LC_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_CS_elev, aes(sample = LC_lm_CS_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_CS_elev$residuals) #not significant with logged transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_CS_elev, aes(x = LC_lm_CS_elev$fitted.values, y = LC_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_CS_elev)

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_slr_DBH <- lm(DBH_ag ~ Elevation..m.FIXED, data = LC_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
LC_slr_DBH_cooks <- cooks.distance(LC_slr_DBH) #calculating the Cook's D for each point
plot(LC_slr_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_slr_DBH_cooks[(LC_slr_DBH_cooks > (3 * mean(LC_slr_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_dbh_no_outliers <- LC_fixed_field_data_processed[-c(3, 33, 38, 48, 51, 54, 56, 58, 60, 61,
                                                                                  64, 67, 69, 72, 73, 75, 77, 84, 85, 99, 119, 126, 161),]

#creating the linear regressions

#linear regression without transformations
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag_lg ~ LC_fixed_field_data_processed$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed$DBH_ag_sqrt ~ LC_fixed_field_data_processed$Elevation..m.)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed_dbh_no_outliers$DBH_ag ~ LC_fixed_field_data_processed_dbh_no_outliers$Elevation..m.)

#linear regression with log transformation of response variable
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed_dbh_no_outliers$DBH_ag_lg ~ LC_fixed_field_data_processed_dbh_no_outliers$Elevation..m.)

#linear regression with square root transformation of response variable
LC_lm_DBH_elev  <- lm(LC_fixed_field_data_processed_dbh_no_outliers$DBH_ag_sqrt ~ LC_fixed_field_data_processed_dbh_no_outliers$Elevation..m.)

#plotting the scatterplot and linear model in ggplot with a square root transformation of DBH
ggplot(data = LC_fixed_field_data_processed, (aes(x=Elevation..m., y=DBH_ag_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("sqrt(DBH)")

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(LC_lm_DBH_elev, aes(x= LC_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_DBH_elev, aes(sample = LC_lm_DBH_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_DBH_elev$residuals) #not significant with square root transformation and ithout outliers

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_DBH_elev, aes(x = LC_lm_DBH_elev$fitted.values, y = LC_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_DBH_elev)

#non-parametric Mann-Kendall Test, non-parametric test
LC_tau_result_DBH <- cor.test(LC_fixed_field_data_processed$Elevation..m.FIXED, LC_fixed_field_data_processed$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(LC_tau_result_DBH)

# Calculating the trend line for plotting
LC_trend_line_DBH <- predict(loess(LC_fixed_field_data_processed$DBH_ag ~ LC_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_processed$Elevation..m.FIXED, y = (LC_fixed_field_data_processed$DBH_ag), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_processed$Elevation..m.FIXED, y = LC_trend_line_DBH), color = "red") +
  labs(x = "Elevation", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()


## SD linear models

#dropping the NAs from SD elevation values
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  drop_na(Elevation..m.FIXED)

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_slr_SCA <- lm(Canopy_short ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
SD_slr_SCA_cooks <- cooks.distance(SD_slr_SCA) #calculating the Cook's D for each point
plot(SD_slr_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_slr_SCA_cooks[(SD_slr_SCA_cooks > (3 * mean(SD_slr_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_sca_no_outliers <- SD_fixed_field_data_processed[-c(1, 49, 69, 70, 74, 83, 84, 86,
                                                                                  97, 117, 121, 122, 155, 156, 158, 210),]

#creating the linear regressions

#linear regression without transformations
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed$Canopy_short_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed_sca_no_outliers$Canopy_short ~ SD_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed_sca_no_outliers$Canopy_short_lg ~ SD_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_sca_elev  <- lm(SD_fixed_field_data_processed_sca_no_outliers$Canopy_short_sqrt ~ SD_fixed_field_data_processed_sca_no_outliers$Elevation..m.FIXED)

#plotting the scatterplot and linear model in ggplot with the logged transformation of SCA
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_short_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("log(Short Canopy Axis)")

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  #histogram of the residuals
ggplot(SD_lm_sca_elev, aes(x= SD_lm_sca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_elev, aes(sample = SD_lm_sca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_sca_elev$residuals) #significant with all transformations and removing of outliers

#will need to use a non-parametric test of association because I cannot get it to pass the normality condition

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_sca_elev, aes(x = SD_lm_sca_elev$fitted.values, y = SD_lm_sca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_sca_elev)

#correlation test
cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$Canopy_short)

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_SCA <- cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$Canopy_short,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(SD_tau_result_SCA)

# Calculating the trend line for plotting
SD_trend_line_SCA <- predict(loess(SD_fixed_field_data_processed$Canopy_short ~ SD_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = (SD_fixed_field_data_processed$Canopy_short), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = SD_trend_line_SCA), color = "red") +
  labs(x = "Elevation", y = "Short Canopy Area", title = "Trend Line Plot") +
  theme_minimal()


#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_slr_LCA <- lm(Canopy_long ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
SD_slr_LCA_cooks <- cooks.distance(SD_slr_LCA) #calculating the Cook's D for each point
plot(SD_slr_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_slr_LCA_cooks[(SD_slr_LCA_cooks > (3 * mean(SD_slr_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_lca_no_outliers <- SD_fixed_field_data_processed[-c(2, 3, 70, 74, 83, 84, 85, 86, 88, 97,
                                                                                  117, 121, 122, 157, 191),]

#creating the linear regressions

#linear regression without transformations
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed$Canopy_long_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers 

#linear regression without transformations
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed_lca_no_outliers$Canopy_long ~ SD_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed_lca_no_outliers$Canopy_long_lg ~ SD_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_lca_elev  <- lm(SD_fixed_field_data_processed_lca_no_outliers$Canopy_long_sqrt ~ SD_fixed_field_data_processed_lca_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_lca_elev, aes(x= SD_lm_lca_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_lca_elev, aes(sample = SD_lm_lca_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_lca_elev$residuals) # not sign when using a logged transformation or a square root one, but the normality condition was met better with square rooting

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_lca_elev, aes(x = SD_lm_lca_elev$fitted.values, y = SD_lm_lca_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_lca_elev)

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_slr_CA <- lm(Canopy_area ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
SD_slr_CA_cooks <- cooks.distance(SD_slr_CA) #calculating the Cook's D for each point
plot(SD_slr_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_slr_CA_cooks[(SD_slr_CA_cooks > (3 * mean(SD_slr_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_ca_no_outliers <- SD_fixed_field_data_processed[-c(1, 3, 49, 74, 83, 84, 85, 86, 97,
                                                                                  117, 121, 122, 155, 156, 157, 158, 210),]


#creating the linear regressions

#linear regression without transformations
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed$Canopy_area_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed_ca_no_outliers$Canopy_area ~ SD_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed_ca_no_outliers$Canopy_area_lg ~ SD_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_CA_elev  <- lm(SD_fixed_field_data_processed_ca_no_outliers$Canopy_area_sqrt ~ SD_fixed_field_data_processed_ca_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_CA_elev, aes(x= SD_lm_CA_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_CA_elev, aes(sample = SD_lm_CA_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_CA_elev$residuals) #not sign with any transformations or removal of outliers

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_CA_elev, aes(x = SD_lm_CA_elev$fitted.values, y = SD_lm_CA_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_CA_elev)

#because the conditions for the slope test were not met, we can use the non-parametric Mann-Kendall's test to check for significant associations

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_CA <- cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$Canopy_area,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(SD_tau_result_CA)

# Calculating the trend line for plotting
SD_trend_line_CA <- predict(loess(SD_fixed_field_data_processed$Canopy_area ~ SD_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = (SD_fixed_field_data_processed$Canopy_area), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = SD_trend_line_CA), color = "red") +
  labs(x = "Elevation", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()


#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_slr_CS <- lm(Crown_spread ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
SD_slr_CS_cooks <- cooks.distance(SD_slr_CS) #calculating the Cook's D for each point
plot(SD_slr_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_slr_CS_cooks[(SD_slr_CS_cooks > (3 * mean(SD_slr_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_cs_no_outliers <- SD_fixed_field_data_processed[-c(1, 2, 3, 69, 70, 74, 83, 84, 85, 86, 97,
                                                                                 117, 121, 122, 158, 210),]

#creating the linear regressions

#linear regression without transformations
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed$Crown_spread_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed_cs_no_outliers$Crown_spread ~ SD_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed_cs_no_outliers$Crown_spread_lg ~ SD_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_CS_elev  <- lm(SD_fixed_field_data_processed_cs_no_outliers$Crown_spread_sqrt ~ SD_fixed_field_data_processed_cs_no_outliers$Elevation..m.FIXED)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_CS_elev, aes(x= SD_lm_CS_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_CS_elev, aes(sample = SD_lm_CS_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_CS_elev$residuals) #significant without outlies and all transformations, need to use non-parametric Mann-Kendall's test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_CS_elev, aes(x = SD_lm_CS_elev$fitted.values, y = SD_lm_CS_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_CS_elev)

#because the conditions for the slope test were not met, we can use the non-parametric Mann-Kendall's test to check for significant associations

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_CS <- cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$Crown_spread,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(SD_tau_result_CS)

# Calculating the trend line for plotting
SD_trend_line_CS <- predict(loess(SD_fixed_field_data_processed$Crown_spread ~ SD_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = (SD_fixed_field_data_processed$Crown_spread), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = SD_trend_line_CS), color = "red") +
  labs(x = "Elevation", y = "Crown Spread", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_slr_DBH <- lm(DBH_ag ~ Elevation..m.FIXED, data = SD_fixed_field_data_processed) #creating a linear regression to use to calculate the Cook's D
SD_slr_DBH_cooks <- cooks.distance(SD_slr_DBH) #calculating the Cook's D for each point
plot(SD_slr_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_slr_DBH_cooks[(SD_slr_DBH_cooks > (3 * mean(SD_slr_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_dbh_no_outliers <- SD_fixed_field_data_processed[-c(as.numeric(names(influential))),]

#creating the linear regressions

#linear regression without transformations
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag_lg ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed$DBH_ag_sqrt ~ SD_fixed_field_data_processed$Elevation..m.FIXED)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed_dbh_no_outliers$DBH_ag ~ SD_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with log transformation of response variable
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed_dbh_no_outliers$DBH_ag_lg ~ SD_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#linear regression with square root transformation of response variable
SD_lm_DBH_elev  <- lm(SD_fixed_field_data_processed_dbh_no_outliers$DBH_ag_sqrt ~ SD_fixed_field_data_processed_dbh_no_outliers$Elevation..m.FIXED)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_DBH_elev, aes(x= SD_lm_DBH_elev$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_DBH_elev, aes(sample = SD_lm_DBH_elev$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_DBH_elev$residuals) #significant without outliers and all transformations

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_DBH_elev, aes(x = SD_lm_DBH_elev$fitted.values, y = SD_lm_DBH_elev$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_DBH_elev)

#because the conditions for the slope test were not met, we can use the non-parametric Mann-Kendall's test to check for significant associations

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_DBH <- cor.test(SD_fixed_field_data_processed$Elevation..m.FIXED, SD_fixed_field_data_processed$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(SD_tau_result_DBH)

# Calculating the trend line for plotting
SD_trend_line_DBH <- predict(loess(SD_fixed_field_data_processed$DBH_ag ~ SD_fixed_field_data_processed$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = (SD_fixed_field_data_processed$DBH_ag), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed$Elevation..m.FIXED, y = SD_trend_line_DBH), color = "red") +
  labs(x = "Elevation", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()


#### Sizes vs. Slope ####

# linear models comparing slope to size/shape

# For all populations/each population and size/shape metric we created single variable regressions by...
      #a) creating the single variable linear regressions with the un-transformed response variable, logged variable, and 
              # square root of the variable, respectively, either with or without outliers
      #b) testing which model best satisfies the conditions for the analysis: LINES
             # Linearity, Independence, Normality of residuals, Equal variance of residuals, and simple random sample 
                    # we tested Linearity by looking at the scatterplots,
                    # we tested Independence by thinking about the explanatory and response variables across the points,
                    # we tested Normality of Residuals using histograms, qq norm plots, and the Shapiro-Wilk's test,
                    # we tested Equal Variance of Residuals using a fitted vs. residuals plot,
                    # we tested Simple Random Sample by thinkg about the data collection method.
      #c) We then looked for significant associations (slopes/correlations)
            # 1) If the LINES conditions are met...
                    # we ran a slope test and a Pearson's correlation test to see if there is a significant association
            # 2) If the LINES conditions are not met...
                    # we ran a Mann-Kendall test (non-parametric test) to look for a significant correlation/tau 


#all points 

#removing NAs in SCA and slope from all points to run tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_slope_raster_15_data_pts)

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x= all_points_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_SCA <- lm(Canopy_short ~ all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
all_points_slr_SCA_cooks <- cooks.distance(all_points_slr_SCA) #calculating the Cook's D for each point
plot(all_points_slr_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_SCA_cooks[(all_points_slr_SCA_cooks > (3 * mean(all_points_slr_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
all_points_fixed_field_data_processed_terrain_sca_no_outliers <- all_points_fixed_field_data_processed_terrain[-c(45, 164, 170, 178, 209, 211, 213, 220, 243, 253, 271, 305,
                                                                                                                  306, 307, 309, 319, 334, 336, 339, 344, 358, 359, 364, 378, 379,
                                                                                                                  380, 381, 385, 431, 435, 451, 474, 482, 496, 516, 574, 615, 629, 630, 634, 635,
                                                                                                                  636, 637, 638, 641, 645, 647, 650, 662),]


#creating the linear regressions

#linear regression without transformations
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain_sca_no_outliers$Canopy_short ~ all_points_fixed_field_data_processed_terrain_sca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain_sca_no_outliers$Canopy_short_lg ~ all_points_fixed_field_data_processed_terrain_sca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_sca_slope  <- lm(all_points_fixed_field_data_processed_terrain_sca_no_outliers$Canopy_short_sqrt ~ all_points_fixed_field_data_processed_terrain_sca_no_outliers$all_points_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_sca_slope, aes(x= all_points_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_sca_slope, aes(sample = all_points_lm_sca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_sca_slope$residuals) #never not significant without outliers or with transformations, use the mann-kendall test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = all_points_lm_sca_slope, aes(x = all_points_lm_sca_slope$fitted.values, y = all_points_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_sca_slope)

#correlation test
cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_short)

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_SCA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_short,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_SCA)

# Calculating the trend line for plotting
LC_trend_line_LCA <- predict(loess(LC_fixed_field_data_all_focal_trees$Canopy_long ~ LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_LCA <- lm(Canopy_long ~ all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
all_points_slr_LCA_cooks <- cooks.distance(all_points_slr_LCA) #calculating the Cook's D for each point
plot(all_points_slr_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_LCA_cooks[(all_points_slr_LCA_cooks > (3 * mean(all_points_slr_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
all_points_fixed_field_data_processed_terrain_lca_no_outliers <- all_points_fixed_field_data_processed_terrain[-c(127, 170, 209, 211, 213, 218, 219, 220, 243, 296, 305,
                                                                                                                  306, 307, 308, 309, 310, 319, 334, 339, 344, 355, 356, 358, 359, 360, 362, 364, 378, 379,
                                                                                                                  380, 381, 383, 385, 400, 435, 451, 482, 496, 516, 574, 615, 627, 629, 631, 632, 634, 635,
                                                                                                                  641, 642, 645, 646, 647, 650),]

#creating the linear regressions

#linear regression without transformations
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain_lca_no_outliers$Canopy_long ~ all_points_fixed_field_data_processed_terrain_lca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain_lca_no_outliers$Canopy_long_lg ~ all_points_fixed_field_data_processed_terrain_lca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_lca_slope  <- lm(all_points_fixed_field_data_processed_terrain_lca_no_outliers$Canopy_long_sqrt ~ all_points_fixed_field_data_processed_terrain_lca_no_outliers$all_points_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_lca_slope, aes(x= all_points_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_lca_slope, aes(sample = all_points_lm_lca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_lca_slope$residuals) #not significant without outliers and square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = all_points_lm_lca_slope, aes(x = all_points_lm_lca_slope$fitted.values, y = all_points_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_lca_slope)

#does not meet the condition of homoscedasticity, so we need to use a non-parametric Mann-Kendall's test

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_LCA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_long,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_LCA)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = (all_points_fixed_field_data_processed_terrain$Canopy_long), color = "blue")) +
  geom_line(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = all_points_trend_line_LCA), color = "red") +
  labs(x = "Slope", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results 
all_points_slr_CA <- lm(Canopy_area ~ all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
all_points_slr_CA_cooks <- cooks.distance(all_points_slr_CA) #calculating the Cook's D for each point
plot(all_points_slr_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_CA_cooks[(all_points_slr_CA_cooks > (3 * mean(all_points_slr_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
all_points_fixed_field_data_processed_terrain_ca_no_outliers <- all_points_fixed_field_data_processed_terrain[-c(45, 90, 164, 170, 178, 209, 211, 213, 218, 220, 253, 255, 271, 305,
                                                                                                                  306, 307, 309, 339, 344, 364, 377, 378, 379,
                                                                                                                  380, 381, 383, 385, 386, 400, 431, 435, 451, 474, 482, 496, 
                                                                                                                  516, 574, 615, 629, 635,
                                                                                                                  641, 645, 646, 647, 650, 662),]

#creating the linear regressions

#linear regression without transformations
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain_ca_no_outliers$Canopy_area ~ all_points_fixed_field_data_processed_terrain_ca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain_ca_no_outliers$Canopy_area_lg ~ all_points_fixed_field_data_processed_terrain_ca_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_CA_slope  <- lm(all_points_fixed_field_data_processed_terrain_ca_no_outliers$Canopy_area_sqrt ~ all_points_fixed_field_data_processed_terrain_ca_no_outliers$all_points_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_CA_slope, aes(x= all_points_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_CA_slope, aes(sample = all_points_lm_CA_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_CA_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = all_points_lm_CA_slope, aes(x = all_points_lm_CA_slope$fitted.values, y = all_points_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_CA_slope)

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_CA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_area,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_CA)

# Calculating the trend line for plotting
all_points_trend_line_CA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = (all_points_fixed_field_data_processed_terrain$Canopy_area), color = "blue")) +
  geom_line(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = all_points_trend_line_CA), color = "red") +
  labs(x = "Slope", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()


#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_CS <- lm(Crown_spread ~ all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
all_points_slr_CS_cooks <- cooks.distance(all_points_slr_CS) #calculating the Cook's D for each point
plot(all_points_slr_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_CS_cooks[(all_points_slr_CS_cooks > (3 * mean(all_points_slr_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
all_points_fixed_field_data_processed_terrain_cs_no_outliers <- all_points_fixed_field_data_processed_terrain[-c(45, 164, 170, 178, 209, 211, 213, 218, 220, 243, 253, 271, 296, 305,
                                                                                                                 306, 307, 308, 309, 319, 334, 336, 339, 344, 355, 358, 359,
                                                                                                                 360, 362, 364, 366, 378, 379,
                                                                                                                 380, 381, 383, 385, 400, 435, 451, 482, 496, 
                                                                                                                 516, 574, 615, 627, 629, 631, 632, 635,
                                                                                                                 641, 645, 646, 647, 650, 662),]
#creating the linear regressions

#linear regression without transformations
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain_cs_no_outliers$Crown_spread ~ all_points_fixed_field_data_processed_terrain_cs_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain_cs_no_outliers$Crown_spread_lg ~ all_points_fixed_field_data_processed_terrain_cs_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_CS_slope  <- lm(all_points_fixed_field_data_processed_terrain_cs_no_outliers$Crown_spread_sqrt ~ all_points_fixed_field_data_processed_terrain_cs_no_outliers$all_points_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_CS_slope, aes(x= all_points_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_CS_slope, aes(sample = all_points_lm_CS_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(all_points_lm_CS_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = all_points_lm_CS_slope, aes(x = all_points_lm_CS_slope$fitted.values, y = all_points_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_CS_slope)

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_CS <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Crown_spread,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_CS)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = (LC_fixed_field_data_all_focal_trees$Canopy_long), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_LCA_over_distance, y = LC_trend_line_LCA), color = "red") +
  labs(x = "LCA over Distance", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
all_points_slr_DBH <- lm(DBH_ag ~ all_points_slope_raster_15_data_pts, data = all_points_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
all_points_slr_DBH_cooks <- cooks.distance(all_points_slr_DBH) #calculating the Cook's D for each point
plot(all_points_slr_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- all_points_slr_DBH_cooks[(all_points_slr_DBH_cooks > (3 * mean(all_points_slr_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
all_points_fixed_field_data_processed_terrain_dbh_no_outliers <- all_points_fixed_field_data_processed_terrain[-c(16, 50, 54, 80, 90, 97, 126, 153, 161, 166, 
                                                                                                                  170, 178, 211, 213, 220, 245, 250, 251, 253, 254, 258,
                                                                                                                  259, 284, 286, 287, 305,
                                                                                                                 306, 307, 309, 330, 339, 344, 355, 358, 359,
                                                                                                                 360, 362, 364, 376, 
                                                                                                                 380, 381, 400, 451, 466, 496, 
                                                                                                                 470, 496, 615, 628, 629, 630, 631, 635, 636, 641, 645, 646, 
                                                                                                                 647, 640, 662),]
#creating the linear regressions

#linear regression without transformations
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag_lg ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain_dbh_no_outliers$DBH_ag ~ all_points_fixed_field_data_processed_terrain_dbh_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain_dbh_no_outliers$DBH_ag_lg ~ all_points_fixed_field_data_processed_terrain_dbh_no_outliers$all_points_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
all_points_lm_DBH_slope  <- lm(all_points_fixed_field_data_processed_terrain_dbh_no_outliers$DBH_ag_sqrt ~ all_points_fixed_field_data_processed_terrain_dbh_no_outliers$all_points_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(all_points_lm_DBH_slope, aes(x= all_points_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(all_points_lm_DBH_slope, aes(sample = all_points_lm_DBH_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(all_points_lm_DBH_slope$residuals) #shapiro-wilk test, all versions of the models are sig meaning we must use a Mann-Kendall test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = all_points_lm_DBH_slope, aes(x = all_points_lm_DBH_slope$fitted.values, y = all_points_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(all_points_lm_DBH_slope)

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_DBH <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_DBH)

# Calculating the trend line for plotting
all_points_trend_line_DBH <- predict(loess(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = (all_points_fixed_field_data_processed_terrain$DBH_ag), color = "blue")) +
  geom_line(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = all_points_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

# LM 

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x= LM_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_lm_focal_SCA <- lm(Canopy_short ~ LM_slope_raster_15_data_pts, data = LM_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LM_lm_focal_SCA_cooks <- cooks.distance(LM_lm_focal_SCA) #calculating the Cook's D for each point
plot(LM_lm_focal_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_lm_focal_SCA_cooks[(LM_lm_focal_SCA_cooks > (3 * mean(LM_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_terrain_no_sca_outliers <- LM_fixed_field_data_processed_terrain[-c(45, 116, 118, 119, 126,
                                                                                              151, 152, 164, 170, 178, 209, 211, 213, 220),]

#creating the linear regressions

#linear regression without transformations
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short ~ LM_fixed_field_data_processed_terrain_no_sca_outliers$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_lg ~ LM_fixed_field_data_processed_terrain_no_sca_outliers$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_sca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_sqrt ~ LM_fixed_field_data_processed_terrain_no_sca_outliers$LM_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LM_lm_sca_slope, aes(x= LM_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_sca_slope, aes(sample = LM_lm_sca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(LM_lm_sca_slope$residuals) #shapiro wilk test not significant when using square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_sca_slope, aes(x = LM_lm_sca_slope$fitted.values, y = LM_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_sca_slope)

#correlation test
cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, LM_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_lm_focal_LCA <- lm(Canopy_long ~ LM_slope_raster_15_data_pts, data = LM_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LM_lm_focal_LCA_cooks <- cooks.distance(LM_lm_focal_LCA) #calculating the Cook's D for each point
plot(LM_lm_focal_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_lm_focal_LCA_cooks[(LM_lm_focal_LCA_cooks > (3 * mean(LM_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_terrain_no_lca_outliers <- LM_fixed_field_data_processed_terrain[-c(43, 50, 106, 116, 127,
                                                                                                  164, 170, 178, 209, 211, 213, 218, 219, 220),]

#creating the linear regressions

#linear regression without transformations
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long ~ LM_fixed_field_data_processed_terrain_no_lca_outliers$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_lg ~ LM_fixed_field_data_processed_terrain_no_lca_outliers$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_lca_slope  <- lm(LM_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_sqrt ~ LM_fixed_field_data_processed_terrain_no_lca_outliers$LM_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LM_lm_lca_slope, aes(x= LM_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_lca_slope, aes(sample = LM_lm_lca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_lca_slope$residuals) #none of the models without outliers and with transformations have non-significant shapiro wilks test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_lca_slope, aes(x = LM_lm_lca_slope$fitted.values, y = LM_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_lca_slope)

#non-parametric Mann-Kendall Test, non-parametric test for the version without outliers
all_points_tau_result_LCA <- cor.test(all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, all_points_fixed_field_data_processed_terrain$Canopy_long,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value
print(all_points_tau_result_LCA)

# Calculating the trend line for plotting
all_points_trend_line_LCA <- predict(loess(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = (all_points_fixed_field_data_processed_terrain$Canopy_long), color = "blue")) +
  geom_line(aes(x = all_points_fixed_field_data_processed_terrain$all_points_slope_raster_15_data_pts, y = all_points_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "Long Canopy Axis", title = "Trend Line Plot") +
  theme_minimal()

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y = Canopy_area_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_lm_focal_CA <- lm(Canopy_area ~ LM_slope_raster_15_data_pts, data = LM_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LM_lm_focal_CA_cooks <- cooks.distance(LM_lm_focal_CA) #calculating the Cook's D for each point
plot(LM_lm_focal_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_lm_focal_CA_cooks[(LM_lm_focal_CA_cooks > (3 * mean(LM_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_terrain_no_ca_outliers <- LM_fixed_field_data_processed_terrain[-c(45, 90, 116, 126, 127,
                                                                                                  164, 170, 178, 209, 211, 213, 218, 220),]


#creating the linear regressions

#linear regression without transformations
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area ~ LM_fixed_field_data_processed_terrain_no_ca_outliers$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain_no_ca_outliers$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_CA_slope  <- lm(LM_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_sqrt ~ LM_fixed_field_data_processed_terrain_no_ca_outliers$LM_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LM_lm_CA_slope, aes(x= LM_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_CA_slope, aes(sample = LM_lm_CA_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_CA_slope$residuals) #non-significant when using without outliers and using square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_CA_slope, aes(x = LM_lm_CA_slope$fitted.values, y = LM_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_CA_slope)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_CA <- cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, LM_fixed_field_data_processed_terrain$Canopy_area_lg,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value 
print(LM_tau_result_CA)

# Calculating the trend line for plotting
LM_trend_line_CA <- predict(loess(LM_fixed_field_data_processed_terrain$Canopy_area_lg ~ LM_fixed_field_data_processed_terrain$sum_CA_over_distance))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = (LC_fixed_field_data_all_focal_trees$Crown_spread), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_all_focal_trees$sum_CS_over_distance, y = LC_trend_line_CS), color = "red") +
  labs(x = "CS over Distance", y = "Crown Spread ", title = "Trend Line Plot") +
  theme_minimal()

#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_lm_focal_CS <- lm(Crown_spread ~ LM_slope_raster_15_data_pts, data = LM_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LM_lm_focal_CS_cooks <- cooks.distance(LM_lm_focal_CS) #calculating the Cook's D for each point
plot(LM_lm_focal_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_lm_focal_CS_cooks[(LM_lm_focal_CS_cooks > (3 * mean(LM_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_terrain_no_cs_outliers <- LM_fixed_field_data_processed_terrain[-c(45, 116, 118, 126, 127, 152,
                                                                                                 164, 170, 178, 209, 211, 213, 218, 220),]


#creating the linear regressions

#linear regression without transformations
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation 
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation 
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread ~ LM_fixed_field_data_processed_terrain_no_cs_outliers$LM_slope_raster_15_data_pts)

#linear regression with log transformation 
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_lg ~ LM_fixed_field_data_processed_terrain_no_cs_outliers$LM_slope_raster_15_data_pts)

#linear regression with square root transformation 
LM_lm_CS_slope  <- lm(LM_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_sqrt ~ LM_fixed_field_data_processed_terrain_no_cs_outliers$LM_slope_raster_15_data_pts)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LM_lm_CS_slope, aes(x= LM_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_CS_slope, aes(sample = LM_lm_CS_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_CS_slope$residuals) #not significant with a square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_CS_slope, aes(x = LM_lm_CS_slope$fitted.values, y = LM_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_CS_slope)

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain, (aes(x=LM_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
LM_lm_focal_DBH <- lm(DBH_ag ~ LM_slope_raster_15_data_pts, data = LM_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LM_lm_focal_DBH_cooks <- cooks.distance(LM_lm_focal_DBH) #calculating the Cook's D for each point
plot(LM_lm_focal_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LM_lm_focal_DBH_cooks[(LM_lm_focal_DBH_cooks > (3 * mean(LM_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LM_fixed_field_data_processed_terrain_no_dbh_outliers <- LM_fixed_field_data_processed_terrain[-c(16, 54, 97, 116, 126, 127, 153,
                                                                                                 161, 164, 165, 166, 170, 178, 211, 213, 215, 220),]

#creating the linear regressions

#linear regression without transformations
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_lg ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag ~ LM_fixed_field_data_processed_terrain_no_dbh_outliers$LM_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_lg ~ LM_fixed_field_data_processed_terrain_no_dbh_outliers$LM_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LM_lm_DBH_slope  <- lm(LM_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_sqrt ~ LM_fixed_field_data_processed_terrain_no_dbh_outliers$LM_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LM_lm_DBH_slope, aes(x= LM_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LM_lm_DBH_slope, aes(sample = LM_lm_DBH_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LM_lm_DBH_slope$residuals) #creating the linear regressions without any outliers and with transformations they are all significant

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LM_lm_DBH_slope, aes(x = LM_lm_DBH_slope$fitted.values, y = LM_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LM_lm_DBH_slope)

#non-parametric Mann-Kendall Test, non-parametric test
LM_tau_result_DBH <- cor.test(LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, LM_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value 
print(LM_tau_result_DBH)

# Calculating the trend line for plotting
LM_trend_line_DBH <- predict(loess(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, y = (LM_fixed_field_data_processed_terrain$DBH_ag), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, y = LM_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

#LC linear models

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_lm_focal_SCA <- lm(Canopy_short ~ LC_slope_raster_15_data_pts, data = LC_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LC_lm_focal_SCA_cooks <- cooks.distance(LC_lm_focal_SCA) #calculating the Cook's D for each point
plot(LC_lm_focal_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_lm_focal_SCA_cooks[(LC_lm_focal_SCA_cooks > (3 * mean(LC_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_terrain_no_sca_outliers <- LC_fixed_field_data_processed_terrain[-c(34, 48, 54, 77, 85, 90, 91, 94, 99, 105, 108, 119, 117),]

#creating the linear regressions

#linear regression without transformations
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with logged transformation of short canopy axis
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation 
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short ~ LC_fixed_field_data_processed_terrain_no_sca_outliers$LC_slope_raster_15_data_pts)

#linear regression with logged transformation of short canopy axis
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_lg ~ LC_fixed_field_data_processed_terrain_no_sca_outliers$LC_slope_raster_15_data_pts)

#linear regression with square root transformation 
LC_lm_sca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_sqrt ~ LC_fixed_field_data_processed_terrain_no_sca_outliers$LC_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LC_lm_sca_slope, aes(x= LC_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_sca_slope, aes(sample = LC_lm_sca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_sca_slope$residuals) #not significant with a log transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_sca_slope, aes(x = LC_lm_sca_slope$fitted.values, y = LC_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_sca_slope)

#correlation test
cor.test(LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts, LC_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Canopy_long_lg)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_lm_focal_LCA <- lm(Canopy_long ~ LC_slope_raster_15_data_pts, data = LC_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LC_lm_focal_LCA_cooks <- cooks.distance(LC_lm_focal_LCA) #calculating the Cook's D for each point
plot(LC_lm_focal_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_lm_focal_LCA_cooks[(LC_lm_focal_LCA_cooks > (3 * mean(LC_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_terrain_no_lca_outliers <- LC_fixed_field_data_processed_terrain[-c(3, 38, 54, 77, 85, 91, 92, 94, 99, 105, 108, 114, 119, 119, 126, 159, 177),]


#creating the linear regressions

#linear regression without transformations
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear transformation with logged long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear transformation with square root transformation long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long ~ LC_fixed_field_data_processed_terrain_no_lca_outliers$LC_slope_raster_15_data_pts)

#linear transformation with logged long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_lg ~ LC_fixed_field_data_processed_terrain_no_lca_outliers$LC_slope_raster_15_data_pts)

#linear transformation with square root transformation long canopy axis
LC_lm_lca_slope  <- lm(LC_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_sqrt ~ LC_fixed_field_data_processed_terrain_no_lca_outliers$LC_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LC_lm_lca_slope, aes(x= LC_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_lca_slope, aes(sample = LC_lm_lca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_lca_slope$residuals) #non significant with a log transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_lca_slope, aes(x = LC_lm_lca_slope$fitted.values, y = LC_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_lca_slope)

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_lm_focal_CA <- lm(Canopy_area ~ LC_slope_raster_15_data_pts, data = LC_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LC_lm_focal_CA_cooks <- cooks.distance(LC_lm_focal_CA) #calculating the Cook's D for each point
plot(LC_lm_focal_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_lm_focal_CA_cooks[(LC_lm_focal_CA_cooks > (3 * mean(LC_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_terrain_no_ca_outliers <- LC_fixed_field_data_processed_terrain[-c(3, 34, 38, 54, 77, 85, 91, 94, 99, 105, 108, 114, 119, 177),]

#creating the linear regressions

#linear regression without transformations
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area ~ LC_fixed_field_data_processed_terrain_no_ca_outliers$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_lg ~ LC_fixed_field_data_processed_terrain_no_ca_outliers$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_CA_slope  <- lm(LC_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_sqrt ~ LC_fixed_field_data_processed_terrain_no_ca_outliers$LC_slope_raster_15_data_pts)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LC_lm_CA_slope, aes(x= LC_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_CA_slope, aes(sample = LC_lm_CA_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test 
shapiro.test(LC_lm_CA_slope$residuals) #non significant with a log transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_CA_slope, aes(x = LC_lm_CA_slope$fitted.values, y = LC_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_CA_slope)

#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Crown Spread")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_lm_focal_CS <- lm(Crown_spread ~ LC_slope_raster_15_data_pts, data = LC_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LC_lm_focal_CS_cooks <- cooks.distance(LC_lm_focal_CS) #calculating the Cook's D for each point
plot(LC_lm_focal_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_lm_focal_CS_cooks[(LC_lm_focal_CS_cooks > (3 * mean(LC_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_terrain_no_cs_outliers <- LC_fixed_field_data_processed_terrain[-c(3, 34, 38, 54, 77, 85, 91, 94, 99, 105, 108, 114, 119, 177),]

#creating the linear regressions

#linear regression without transformations
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread ~ LC_fixed_field_data_processed_terrain_no_cs_outliers$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_lg ~ LC_fixed_field_data_processed_terrain_no_cs_outliers$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_CS_slope  <- lm(LC_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_sqrt ~ LC_fixed_field_data_processed_terrain_no_cs_outliers$LC_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LC_lm_CS_slope, aes(x= LC_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_CS_slope, aes(sample = LC_lm_CS_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_CS_slope$residuals) #not significant with a log transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_CS_slope, aes(x = LC_lm_CS_slope$fitted.values, y = LC_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_CS_slope)

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain, (aes(x=LC_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
LC_lm_focal_DBH <- lm(DBH_ag ~ LC_slope_raster_15_data_pts, data = LC_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
LC_lm_focal_DBH_cooks <- cooks.distance(LC_lm_focal_DBH) #calculating the Cook's D for each point
plot(LC_lm_focal_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- LC_lm_focal_DBH_cooks[(LC_lm_focal_DBH_cooks > (3 * mean(LC_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
LC_fixed_field_data_processed_terrain_no_dbh_outliers <- LC_fixed_field_data_processed_terrain[-c(3, 54, 61, 69, 72, 73, 75, 77, 84, 85, 94, 99, 105, 108, 119, 
                                                                                                 126, 159, 161),]

#creating the linear regressions

#linear regression without transformations
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_lg ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag ~ LC_fixed_field_data_processed_terrain_no_dbh_outliers$LC_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_lg ~ LC_fixed_field_data_processed_terrain_no_dbh_outliers$LC_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
LC_lm_DBH_slope  <- lm(LC_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_sqrt ~ LC_fixed_field_data_processed_terrain_no_dbh_outliers$LC_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(LC_lm_DBH_slope, aes(x= LC_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(LC_lm_DBH_slope, aes(sample = LC_lm_DBH_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_lm_DBH_slope$residuals) #non significant without outliers

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = LC_lm_DBH_slope, aes(x = LC_lm_DBH_slope$fitted.values, y = LC_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(LC_lm_DBH_slope)

#non-parametric Mann-Kendall Test, non-parametric test
LC_tau_result_DBH <- cor.test(LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts, LC_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value 
print(LC_tau_result_DBH)

# Calculating the trend line for plotting
LM_trend_line_DBH <- predict(loess(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts, y = (LC_fixed_field_data_processed_terrain$DBH_ag), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_processed_terrain$LC_slope_raster_15_data_pts, y = LC_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

#SD linear models

#short canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Short Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_lm_focal_SCA <- lm(Canopy_short ~ SD_slope_raster_15_data_pts, data = SD_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
SD_lm_focal_SCA_cooks <- cooks.distance(SD_lm_focal_SCA) #calculating the Cook's D for each point
plot(SD_lm_focal_SCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_lm_focal_SCA_cooks[(SD_lm_focal_SCA_cooks > (3 * mean(SD_lm_focal_SCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_terrain_no_sca_outliers <- SD_fixed_field_data_processed_terrain[-c(21, 49, 83, 85, 87, 97, 117, 122, 156, 157, 
                                                                                                  159, 163, 202, 208, 211),]
#creating the linear regressions

#linear regression without transformations
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_short_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short ~ SD_fixed_field_data_processed_terrain_no_sca_outliers$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_lg ~ SD_fixed_field_data_processed_terrain_no_sca_outliers$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_sca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_sca_outliers$Canopy_short_sqrt ~ SD_fixed_field_data_processed_terrain_no_sca_outliers$SD_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_sca_slope, aes(x= SD_lm_sca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_sca_slope, aes(sample = SD_lm_sca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_sca_slope$residuals) #non significant with square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_sca_slope, aes(x = SD_lm_sca_slope$fitted.values, y = SD_lm_sca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_sca_slope)

#correlation test
cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, SD_fixed_field_data_processed_terrain$Canopy_short)

#long canopy axis

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Long Canopy Axis")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_lm_focal_LCA <- lm(Canopy_long ~ SD_slope_raster_15_data_pts, data = SD_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
SD_lm_focal_LCA_cooks <- cooks.distance(SD_lm_focal_LCA) #calculating the Cook's D for each point
plot(SD_lm_focal_LCA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_lm_focal_LCA_cooks[(SD_lm_focal_LCA_cooks > (3 * mean(SD_lm_focal_LCA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_terrain_no_lca_outliers <- SD_fixed_field_data_processed_terrain[-c(21, 74, 83, 85, 97, 117, 122, 156, 158, 
                                                                                                  159, 161, 163, 188, 192, 202, 206, 208, 211),]

#creating the linear regressions

#linear regression without transformations
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_long_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long ~ SD_fixed_field_data_processed_terrain_no_lca_outliers$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_lg ~ SD_fixed_field_data_processed_terrain_no_lca_outliers$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_lca_slope  <- lm(SD_fixed_field_data_processed_terrain_no_lca_outliers$Canopy_long_sqrt ~ SD_fixed_field_data_processed_terrain_no_lca_outliers$SD_slope_raster_15_data_pts)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_lca_slope, aes(x= SD_lm_lca_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_lca_slope, aes(sample = SD_lm_lca_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_lca_slope$residuals) #non significant with a log transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_lca_slope, aes(x = SD_lm_lca_slope$fitted.values, y = SD_lm_lca_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_lca_slope)

#canopy area

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y = Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("Canopy Area")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_lm_focal_CA <- lm(Canopy_area ~ SD_slope_raster_15_data_pts, data = SD_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
SD_lm_focal_CA_cooks <- cooks.distance(SD_lm_focal_CA) #calculating the Cook's D for each point
plot(SD_lm_focal_CA_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_lm_focal_CA_cooks[(SD_lm_focal_CA_cooks > (3 * mean(SD_lm_focal_CA_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_terrain_no_ca_outliers <- SD_fixed_field_data_processed_terrain[-c(31, 49, 83, 85, 97, 117, 122, 155, 156, 157, 158, 
                                                                                                  159, 163, 202, 208, 211, 223),]

#creating the linear regressions

#linear regression without transformations
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain$Canopy_area_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area ~ SD_fixed_field_data_processed_terrain_no_ca_outliers$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_lg ~ SD_fixed_field_data_processed_terrain_no_ca_outliers$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_CA_slope  <- lm(SD_fixed_field_data_processed_terrain_no_ca_outliers$Canopy_area_sqrt ~ SD_fixed_field_data_processed_terrain_no_ca_outliers$SD_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_CA_slope, aes(x= SD_lm_CA_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_CA_slope, aes(sample = SD_lm_CA_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_CA_slope$residuals) #significant without outliers and with transformations need to use non-parametric test

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_CA_slope, aes(x = SD_lm_CA_slope$fitted.values, y = SD_lm_CA_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_CA_slope)

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_CA <- cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, SD_fixed_field_data_processed_terrain$Canopy_area,  method = "kendall")

# Print Kendall's tau (a correlation metric) and its associated p-value 
print(SD_tau_result_CA)

# Calculating the trend line for plotting
LM_trend_line_DBH <- predict(loess(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, y = (LM_fixed_field_data_processed_terrain$DBH_ag), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_processed_terrain$LM_slope_raster_15_data_pts, y = LM_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

#crown spread

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=Crown_spread_sqrt)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("sqrt(Crown Spread)")+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

#using Cook's D to check for highly influential points that may skew the linear model results
SD_lm_focal_CS <- lm(Crown_spread ~ SD_slope_raster_15_data_pts, data = SD_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
SD_lm_focal_CS_cooks <- cooks.distance(SD_lm_focal_CS) #calculating the Cook's D for each point
plot(SD_lm_focal_CS_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_lm_focal_CS_cooks[(SD_lm_focal_CS_cooks > (3 * mean(SD_lm_focal_CS_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_terrain_no_cs_outliers <- SD_fixed_field_data_processed_terrain[-c(21, 49, 74, 83, 85, 97, 117, 122, 156, 158, 
                                                                                                 159, 163, 202, 206, 208, 211),]

#creating the linear regressions

#linear regression without transformations
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with logged transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain$Crown_spread_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread ~ SD_fixed_field_data_processed_terrain_no_cs_outliers$SD_slope_raster_15_data_pts)

#linear regression with logged transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_lg ~ SD_fixed_field_data_processed_terrain_no_cs_outliers$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of crown spread
SD_lm_CS_slope  <- lm(SD_fixed_field_data_processed_terrain_no_cs_outliers$Crown_spread_sqrt ~ SD_fixed_field_data_processed_terrain_no_cs_outliers$SD_slope_raster_15_data_pts)


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_CS_slope, aes(x= SD_lm_CS_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_CS_slope, aes(sample = SD_lm_CS_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_CS_slope$residuals) #not significant with square root transformation

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_CS_slope, aes(x = SD_lm_CS_slope$fitted.values, y = SD_lm_CS_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_CS_slope)

#DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain, (aes(x=SD_slope_raster_15_data_pts, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (degrees)")+
  ylab("DBH")

#using Cook's D to check for highly influential points that may skew the linear model results
SD_lm_focal_DBH <- lm(DBH_ag ~ SD_slope_raster_15_data_pts, data = SD_fixed_field_data_processed_terrain) #creating a linear regression to use to calculate the Cook's D
SD_lm_focal_DBH_cooks <- cooks.distance(SD_lm_focal_DBH) #calculating the Cook's D for each point
plot(SD_lm_focal_DBH_cooks, type = 'h') #checking to see which Cook's D are unusually high
influential <- SD_lm_focal_DBH_cooks[(SD_lm_focal_DBH_cooks > (3 * mean(SD_lm_focal_DBH_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
influential

#removing points that were deemed too influential on the linear model fit
SD_fixed_field_data_processed_terrain_no_dbh_outliers <- SD_fixed_field_data_processed_terrain[-c(31, 37, 85, 117, 122, 154, 158, 
                                                                                                 159, 176, 208, 211, 223),]

#creating the linear regressions

#linear regression without transformations
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_lg ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain$DBH_ag_sqrt ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts)

#creating the linear regressions without any outliers

#linear regression without transformations
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag ~ SD_fixed_field_data_processed_terrain_no_dbh_outliers$SD_slope_raster_15_data_pts)

#linear regression with log transformation of response variable
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_lg ~ SD_fixed_field_data_processed_terrain_no_dbh_outliers$SD_slope_raster_15_data_pts)

#linear regression with square root transformation of response variable
SD_lm_DBH_slope  <- lm(SD_fixed_field_data_processed_terrain_no_dbh_outliers$DBH_ag_sqrt ~ SD_fixed_field_data_processed_terrain_no_dbh_outliers$SD_slope_raster_15_data_pts)

#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test  

#histogram of the residuals
ggplot(SD_lm_DBH_slope, aes(x= SD_lm_DBH_slope$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope (degrees)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(SD_lm_DBH_slope, aes(sample = SD_lm_DBH_slope$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(SD_lm_DBH_slope$residuals) #significant with all transformations and without outliers, need non-parametric tests

#checking equal variance with a residuals vs. fitted values plot
ggplot(data = SD_lm_DBH_slope, aes(x = SD_lm_DBH_slope$fitted.values, y = SD_lm_DBH_slope$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope (degrees)")

#Slope Test, obtained from the summary of the linear regression results
summary(SD_lm_DBH_slope)

#non-parametric Mann-Kendall Test, non-parametric test
SD_tau_result_DBH <- cor.test(SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, SD_fixed_field_data_processed_terrain$DBH_ag,  method = "kendall") #DBH_ag_sqrt

# Print Kendall's tau (a correlation metric) and its associated p-value
print(SD_tau_result_DBH)

# Calculating the trend line for plotting
SD_trend_line_DBH <- predict(loess(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, y = (SD_fixed_field_data_processed_terrain$DBH_ag), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed_terrain$SD_slope_raster_15_data_pts, y = SD_trend_line_DBH), color = "red") +
  labs(x = "Slope", y = "DBH", title = "Trend Line Plot") +
  theme_minimal()

#### Sizes vs. Aspect ####

# we ran ANOVAs to test difference in size means between cardinal directions

#ANOVAs assume residuals are normally distributed, homoscedasticity, and independence

# For all populations/each population and size/shape metric we tested for significant relationships with 
   #the 8 or 4 categories of aspect by...
       
   #a) creating box plots of tree size/shape 
   #b) fitting an ANOVA (analysis of variance) model and seeing if there are any pairwise significant comparisons
   #c) checking to see if the ANOVA model meets the conditions (normal distribution of the residuals, homoscedasticity, independence)
            #normal residuals is checked with histograms, qq norm plot, and Shapiro-Wilks test
            #equal variance is tested with Fligner-Killeen test which is more useful when dealing with non-normal residuals and when outliers are present
                     #Levene's Test is also used to check equal variance, but is not super robust to strong differences to normality
                     #Finally, a rule of thumb test is also used to look for equal variance.
   #d) checking for significant difference in means between size/shape between aspect categories
            # 1) if conditions are met...
                      #we used the summary results from the ANOVA model to check for a significant difference in means
            # 2) if conditions are not met...
                      #we used the Kruskal-Wallis test (non-parametric test) to look for a general significant difference in means
                           #and if it is significant, we used the post-hoc Wilcoxon rank sum test to look at pairwise comparisons of
                           #significant differences in mean size/shape between aspect categories

#8 categories for direction


#all points 

#removing NAs in SCA and aspect to allow us to run tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_aspect_raster_15_data_pts_8_categorical)

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
all_points_aov_SCA_aspect_8 <- aov(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_SCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_SCA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


#checking normality of residuals with a histogram, qqnorm plot, and Shapiro-Wilk Test
hist(all_points_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals are not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test

#Fligner-Killeen, more useful than the other equal variance tests when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_SCA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                    all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_SCA, na.rm = T) / min(all_points_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two, the test did not pass

#variances are equally distributed

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
all_points_aov_LCA_aspect_8 <- aov(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_LCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_LCA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_long, 
                                          all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test, if significant, have to run a non-parametric test

#not normally distributed

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_LCA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_LCA, na.rm = T) / min(all_points_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances condition is met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
all_points_aov_CA_aspect_8 <- aov(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_CA_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_area, 
                                         all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_CA <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_CA, na.rm = T) / min(all_points_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances condition is met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
all_points_aov_CS_aspect_8 <- aov(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CS_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_CS_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Crown_spread, 
                                         all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_CS <- tapply(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_CS, na.rm = T) / min(all_points_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
all_points_aov_DBH_aspect_8 <- aov(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_DBH_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_DBH_aspect_8 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$DBH_ag, 
                                          all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(all_points_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

#residuals are not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
all_points_thumb_test_DBH <- tapply(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_DBH, na.rm = T) / min(all_points_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_8_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "none") #version with no p-value adjustment

pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method


# LM

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LM_aov_SCA_aspect_8 <- aov(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect_8) #ANOVA summary

# checking to see if residuals are normal
hist(LM_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LM_thumb_test_SCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_SCA, na.rm = T) / min(LM_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LM_aov_LCA_aspect_8 <- aov(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_LCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_LCA_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_long, 
                                          LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LM_thumb_test_LCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_LCA, na.rm = T) / min(LM_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
LM_aov_CA_aspect_8 <- aov(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_CA_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_area, 
                                         LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LM_thumb_test_CA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_CA, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
LM_aov_CS_aspect_8 <- aov(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CS_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_CS_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Crown_spread, 
                                         LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LM_thumb_test_CS <- tapply(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_CS, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
LM_aov_DBH_aspect_8 <- aov(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_DBH_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_DBH_aspect_8 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$DBH_ag, 
                                          LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(LM_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$DBH_ag ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LM_thumb_test_DBH <- tapply(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_DBH, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ LM_aspect_raster_15_data_pts_8_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# LC

#short canopy axis
 
#had to remove the NAs to be able to run the function
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(LC_aspect_raster_15_data_pts_8_categorical)


#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LC_aov_SCA_aspect_8 <- aov(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
esummary(LC_aov_SCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_SCA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LC_aov_LCA_aspect_8 <- aov(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_LCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_LCA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_long, 
                                          LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LC_thumb_test_LCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_LCA, na.rm = T) / min(LC_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
LC_aov_CA_aspect_8 <- aov(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_CA_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_area, 
                                         LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_t_test_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(LC_t_test_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LC_thumb_test_CA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_CA, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
LC_aov_CS_aspect_8 <- aov(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CS_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_CS_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Crown_spread, 
                                         LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_t_test_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(LC_t_test_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LC_thumb_test_CS <- tapply(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_CS, na.rm = T) / min(LC_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method


#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
LC_aov_DBH_aspect_8 <- aov(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_DBH_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_DBH_aspect_8 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$DBH_ag, 
                                          LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_t_test_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(LC_t_test_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(LC_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ LC_aspect_raster_15_data_pts_8_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# SD

#short canopy axis

#had to remove the NAs to be able to run the function
SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(SD_aspect_raster_15_data_pts_8_categorical)

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
SD_aov_SCA_aspect_8 <- aov(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_SCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_aov_SCA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_SCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_SCA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_SCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
SD_thumb_test_SCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_SCA, na.rm = T) / min(SD_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
SD_aov_LCA_aspect_8 <- aov(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_LCA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_LCA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_long, 
                                          SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_LCA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect")

qqnorm(SD_aov_LCA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_LCA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
SD_thumb_test_LCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_LCA, na.rm = T) / min(SD_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
SD_aov_CA_aspect_8 <- aov(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CA_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_CA_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_area, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_CA_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_CA_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_CA_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
SD_thumb_test_CA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_CA, na.rm = T) / min(SD_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
SD_aov_CS_aspect_8 <- aov(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CS_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_CS_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Crown_spread, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_CS_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_CS_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_CS_aspect_8$residuals) #Shapiro-Wilk test

#residual not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
SD_thumb_test_CS <- tapply(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_CS, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method

#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")+
  ylim(c(0,1.5))+
  theme(
    axis.title.x = element_text(size = 16),
    axis.text.x = element_text(size = 14),
    axis.title.y = element_text(size = 16),
    axis.text.y = element_text(size = 14))

#generating the ANOVA and ANOVA summary
SD_aov_DBH_aspect_8 <- aov(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_DBH_aspect_8) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_DBH_aspect_8 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$DBH_ag, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_DBH_aspect_8$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_DBH_aspect_8$residuals) #qqnorm plot

shapiro.test(SD_aov_DBH_aspect_8$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical)

#Rule of Thumb Test
SD_thumb_test_DBH <- tapply(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_DBH, na.rm = T) / min(SD_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ SD_aspect_raster_15_data_pts_8_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_8_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#4 categories for direction

#all points 

#removing NAs preventing us from running tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_aspect_raster_15_data_pts_4_categorical)

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
all_points_aov_SCA_aspect_4 <- aov(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_SCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_SCA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_short, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categoricalv, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_short ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
all_points_thumb_test_SCA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_SCA_4, na.rm = T) / min(all_points_thumb_test_SCA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variances not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_short, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
all_points_aov_LCA_aspect_4 <- aov(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_LCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_LCA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_long, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_long ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
all_points_thumb_test_LCA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_LCA_4, na.rm = T) / min(all_points_thumb_test_LCA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_long, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
all_points_aov_CA_aspect_4 <- aov(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_CA_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Canopy_area, 
                                                 all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)


#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Canopy_area ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
all_points_thumb_test_CA_4 <- tapply(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_CA_4, na.rm = T) / min(all_points_thumb_test_CA_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Canopy_area, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
all_points_aov_CS_aspect_4 <- aov(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_CS_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_CS_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$Crown_spread, 
                                                 all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(all_points_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal 

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$Crown_spread ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
all_points_thumb_test_CS_4 <- tapply(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_CS_4, na.rm = T) / min(all_points_thumb_test_CS_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$Crown_spread, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = all_points_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = all_points_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
all_points_aov_DBH_aspect_4 <- aov(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)
summary(all_points_aov_DBH_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
all_points_t_test_DBH_aspect_4 <- pairwise.t.test(all_points_fixed_field_data_processed_terrain$DBH_ag, 
                                                  all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(all_points_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(all_points_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(all_points_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(all_points_fixed_field_data_processed_terrain$DBH_ag ~ all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
all_points_thumb_test_DBH_4 <- tapply(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(all_points_thumb_test_DBH_4, na.rm = T) / min(all_points_thumb_test_DBH_4, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ all_points_aspect_raster_15_data_pts_4_categorical, data = all_points_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(all_points_fixed_field_data_processed_terrain$DBH_ag, all_points_fixed_field_data_processed_terrain$all_points_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# LM

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LM_aov_SCA_aspect_4 <- aov(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_SCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_SCA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_short, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_short ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LM_thumb_test_SCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_SCA, na.rm = T) / min(LM_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition is met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_short, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LM_aov_LCA_aspect_4 <- aov(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LCMaov_LCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_LCA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_long, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf") 

# checking to see if residuals are normal
hist(LM_t_test_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LM_t_test_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_long ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LM_thumb_test_LCA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_LCA, na.rm = T) / min(LM_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_long, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
LM_aov_CA_aspect_4 <- aov(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_CA_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Canopy_area, 
                                       LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Canopy_area ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LM_thumb_test_CA <- tapply(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_CA, na.rm = T) / min(LM_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Canopy_area, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
LM_aov_CS_aspect_4 <- aov(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_CS_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_CS_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$Crown_spread, 
                                       LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LM_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$Crown_spread ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LM_thumb_test_CS <- tapply(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_CS, na.rm = T) / min(LM_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$Crown_spread, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = LM_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LM_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
LM_aov_DBH_aspect_4 <- aov(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)
summary(LM_aov_DBH_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LM_t_test_DBH_aspect_4 <- pairwise.t.test(LM_fixed_field_data_processed_terrain$DBH_ag, 
                                        LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LM_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(LM_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(LM_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LM_fixed_field_data_processed_terrain$LM_t_test_DBH_aspect_4 ~ LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LM_thumb_test_DBH <- tapply(LM_fixed_field_data_processed_terrain$LM_t_test_DBH_aspect_4, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LM_thumb_test_DBH, na.rm = T) / min(LM_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ LM_aspect_raster_15_data_pts_4_categorical, data = LM_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LM_fixed_field_data_processed_terrain$DBH_ag, LM_fixed_field_data_processed_terrain$LM_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# LC

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LC_aov_SCA_aspect_4 <- aov(Canopy_short ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_SCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_SCA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_short, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_short ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LC_thumb_test_SCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_SCA, na.rm = T) / min(LC_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_short, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method


#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
LC_aov_LCA_aspect_4 <- aov(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_LCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_LCA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_long, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test

#residuals were not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_long ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LC_thumb_test_LCA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_LCA, na.rm = T) / min(LC_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_long, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method


# canopy area

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
LC_aov_CA_aspect_4 <- aov(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_CA_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Canopy_area, 
                                       LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(LC_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Canopy_area ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LC_thumb_test_CA <- tapply(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_CA, na.rm = T) / min(LC_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Canopy_area, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
LC_aov_CS_aspect_4 <- aov(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_CS_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_CS_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$Crown_spread, 
                                       LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

#residuals not normal

# checking equal variances with Fligner-Killeen Test, Levene's Test, and Rule of Thumb Test 

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$Crown_spread ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LC_thumb_test_CS <- tapply(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_CS, na.rm = T) / min(LC_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$Crown_spread, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method
#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = LC_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
LC_aov_DBH_aspect_4 <- aov(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)
summary(LC_aov_DBH_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
LC_t_test_DBH_aspect_4 <- pairwise.t.test(LC_fixed_field_data_processed_terrain$DBH_ag, 
                                        LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(LC_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(LC_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(LC_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test, sig so we need to use a non-parametric test

#residuals not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(LC_fixed_field_data_processed_terrain$DBH_ag ~ LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
LC_thumb_test_DBH <- tapply(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(LC_thumb_test_DBH, na.rm = T) / min(LC_thumb_test_DBH, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance is met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ LC_aspect_raster_15_data_pts_4_categorical, data = LC_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(LC_fixed_field_data_processed_terrain$DBH_ag, LC_fixed_field_data_processed_terrain$LC_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method

# SD

#short canopy axis

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_short))+
  xlab("Directions")+
  ylab("Short Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
SD_aov_SCA_aspect_4 <- aov(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_SCA_aspect_4) #ANOVA summary
 
#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_SCA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_short, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_SCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Short Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_SCA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_SCA_aspect_4$residuals) #Shapiro-Wilk test, sig so need to use a non-parametric test

#residuals are not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_short ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
SD_thumb_test_SCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_SCA, na.rm = T) / min(SD_thumb_test_SCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_short ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_short, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method

#long canopy axis

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_long))+
  xlab("Directions")+
  ylab("Long Canopy Axis (m)")

#generating the ANOVA and ANOVA summary
SD_aov_LCA_aspect_4 <- aov(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_LCA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_LCA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_long, 
                                          SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_LCA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Long Canopy Axis vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_LCA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_LCA_aspect_4$residuals) #Shapiro-Wilk test, sig so need to use a non-parametric test

#residuals not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_long ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
SD_thumb_test_LCA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_LCA, na.rm = T) / min(SD_thumb_test_LCA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_long ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_long, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method

# canopy area

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Canopy_area))+
  xlab("Directions")+
  ylab("Canopy Area (m2)")

#generating the ANOVA and ANOVA summary
SD_aov_CA_aspect_4 <- aov(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CA_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_CA_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Canopy_area, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_CA_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Canopy Area vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_CA_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_CA_aspect_4$residuals) #Shapiro-Wilk test

#residuals are not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Canopy_area ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
SD_thumb_test_CA <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_CA, na.rm = T) / min(SD_thumb_test_CA, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition not met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Canopy_area ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method

#crown spread

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = Crown_spread))+
  xlab("Directions")+
  ylab("Crown Spread (m2)")

#generating the ANOVA and ANOVA summary
SD_aov_CS_aspect_4 <- aov(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_CS_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_CS_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$Crown_spread, 
                                       SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")

# checking to see if residuals are normal
hist(SD_aov_CS_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for Crown Spread vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_CS_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_CS_aspect_4$residuals) #Shapiro-Wilk test

#residuals are not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$Crown_spread ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
SD_thumb_test_CS <- tapply(SD_fixed_field_data_processed_terrain$Canopy_area, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_CS, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(Crown_spread ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$Crown_spread, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p value adjusted Crown_spread false discovery rate method


#DBH ag

#boxplot of sizes by the directional categories
ggplot(data = SD_fixed_field_data_processed_terrain)+
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_4_categorical, y = DBH_ag))+
  xlab("Directions")+
  ylab("DBH")

#generating the ANOVA and ANOVA summary
SD_aov_DBH_aspect_4 <- aov(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)
summary(SD_aov_DBH_aspect_4) #ANOVA summary

#pairwise t-test to see significant differences between categories, using a bonferonni adjustment to control for multiple testing
SD_t_test_DBH_aspect_4 <- pairwise.t.test(SD_fixed_field_data_processed_terrain$DBH_ag, 
                                        SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, p.adj = "bonf")


# checking to see if residuals are normal
hist(SD_aov_DBH_aspect_4$residuals, xlab = "Residuals", main = "Distribution of Residuals for DBH vs. Aspect") #histogram of the residuals

qqnorm(SD_aov_DBH_aspect_4$residuals) #qqnorm plot

shapiro.test(SD_aov_DBH_aspect_4$residuals) #Shapiro-Wilk test, sig so have to use a non-parametric test

#residuals are not normal

# checking equal variances with levene's test, Fligner-Killeen, and rule of thumb

#Fligner-Killeen, more useful when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#Levene's Test, not super robust to strong differences to normality
leveneTest(SD_fixed_field_data_processed_terrain$DBH_ag ~ SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical)

#Rule of Thumb Test
SD_thumb_test_DBH <- tapply(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical, sd) #calculating the standard deviation for the response variable across each cardinal direction
max(SD_thumb_test_DBH, na.rm = T) / min(SD_thumb_test_CS, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass

#equal variance condition met

#non-parametric tests

#Kruskal-Wallis test
kruskal.test(DBH_ag ~ SD_aspect_raster_15_data_pts_4_categorical, data = SD_fixed_field_data_processed_terrain)

#post-hoc Wilcoxon rank sum tests to check difference in means/medians
pairwise.wilcox.test(SD_fixed_field_data_processed_terrain$DBH_ag, SD_fixed_field_data_processed_terrain$SD_aspect_raster_15_data_pts_4_categorical,
                     p.adjust.method = "fdr") #p-value adjusted using false discovery rate method


#Boxplots to look at LC NW, SW, and N tended to have larger means

#grouping NW, N, SW, and the other directions
LC_fixed_field_data_processed_terrain <- LC_fixed_field_data_processed_terrain %>%
  add_column(LC_aspect_raster_15_data_pts_8_regrouped = NA) %>% #creating a new column to regroup the directions
  mutate(LC_aspect_raster_15_data_pts_8_regrouped = case_when(LC_aspect_raster_15_data_pts_8_categorical == "NE" ~ "Other", 
                                                              LC_aspect_raster_15_data_pts_8_categorical == "E" ~  "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "NW" ~ "NW",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "N" ~ "N",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "W" ~ "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "SW" ~ "SW",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "S" ~ "Other",
                                                              LC_aspect_raster_15_data_pts_8_categorical == "SE" ~ "Other",
                                                              ))

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_short))+
  xlab("Direction")+
  ylab("Short Canopy Axis")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_long))+
  xlab("Direction")+
  ylab("Long Canopy Axis")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Canopy_area))+
  xlab("Direction")+
  ylab("Canopy Area")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = Crown_spread))+
  xlab("Direction")+
  ylab("Crown Spread")

ggplot(LC_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = LC_aspect_raster_15_data_pts_8_regrouped, y = DBH_ag))+
  xlab("Direction")+
  ylab("DBH")


#Boxplot of SD to look at NE, N, and W tended to have larger means

SD_fixed_field_data_processed_terrain <- SD_fixed_field_data_processed_terrain %>%
  add_column(SD_aspect_raster_15_data_pts_8_regrouped = NA) %>% #creating a new column to regroup the directions
  mutate(SD_aspect_raster_15_data_pts_8_regrouped = case_when(SD_aspect_raster_15_data_pts_8_categorical == "NE" ~ "NE", 
                                                              SD_aspect_raster_15_data_pts_8_categorical == "E" ~  "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "NW" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "N" ~ "N",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "W" ~ "W",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "SW" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "S" ~ "Other",
                                                              SD_aspect_raster_15_data_pts_8_categorical == "SE" ~ "Other",
  ))



ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_short))+
  xlab("Direction")+
  ylab("Short Canopy Axis")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_long))+
  xlab("Direction")+
  ylab("Long Canopy Axis")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Canopy_area))+
  xlab("Direction")+
  ylab("Canopy Area")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = Crown_spread))+
  xlab("Direction")+
  ylab("Crown Spread")

ggplot(SD_fixed_field_data_processed_terrain) + #generate the base plot
  geom_boxplot(aes(x = SD_aspect_raster_15_data_pts_8_regrouped, y = DBH_ag))+
  xlab("Direction")+
  ylab("DBH")


# Plotting size of shape/size values for points in population with coloring based on aspect

# LC

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_short, color = LC_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_long, color = LC_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Canopy_area, color = LC_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = Crown_spread, color = LC_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_LC_trans)+
  geom_sf(data = LC_fixed_field_data_processed_terrain, aes(size = DBH_ag, color = LC_aspect_raster_15_data_pts_8_regrouped))


#SD

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_short, color = SD_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_long, color = SD_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Canopy_area, color = SD_aspect_raster_15_data_pts_8_regrouped))

ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = Crown_spread, color = SD_aspect_raster_15_data_pts_8_regrouped))


ggplot()+
  geom_sf(data = river_SD_trans)+
  geom_sf(data = SD_fixed_field_data_processed_terrain, aes(size = DBH_ag, color = SD_aspect_raster_15_data_pts_8_regrouped))

