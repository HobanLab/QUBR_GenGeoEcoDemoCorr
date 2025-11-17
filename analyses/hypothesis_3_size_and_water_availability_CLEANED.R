# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei's size/shape is influenced by either elevation, slope, aspect, and/or distance to rivers%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#the purpose of the script is to determine how the elevation, slope, and distance to river each 
#relate to size/shape of the trees using single linear regressions and how aspect relates to the size/shape
#of the trees using difference in means tests/ANOVAs/Kruskal-Wallis Tests

# A significant relationship could indicate how elevation, slope, aspect, and distance to river may influence the size and shape of the 
# trees and potentially help explain their distribution. 

#The script is broken up into these sections:

  # 1) Loading and processing the packages and processed data for the trees, topography, soil metrics, and distance to river in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
          # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
              #generating river and point buffers and bounding boxes,
          # Extracting and processing slope, elevation, and aspect (4 and 8 cardinal directions) data using 15 m res rasters,
          # Extracting the distance to the river of each tree for each population,
    # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
  # 2) Comparing size of the trees to their elevation to look for a relationship using single linear regression
  # 3) Comparing size of the trees to their slope to look for a relationship using single linear regression
  # 4) Comparing size of the trees to their aspect to look for a relationship using ANOVA / Kruskal-Wallis Tests
  # 5) Making a function for creating the distance to river single linear regressions, finding the best models (with transformations and removal of outliers or not) 
          # that meet the conditions for SLRs
  # 6) Running the distance to river SLR function and storing the outputs

# NOTE: Uncomment and run line 47, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

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
# NOTE: Uncomment and run line 47, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")


#### Creating a Function for Simple Linear Regressions ####

#comparing distance to river's edge to size values

# For all populations/each population and size/shape metric we created single variable regressions by...
#a) creating the single variable linear regressions with the un-transformed response variable, logged variable, and 
# square root of the variable, respectively, either with or without outliers
#b) testing which model best satisfies the conditions for the analysis: LINES
    # Linearity, Independence, Normality of residuals, Equal variance of residuals, and simple random sample 
      # we tested Linearity by looking at the scatterplots,
      # we tested Independence by thinking about the explanatory and response variables across the points,
      # we tested Normality of Residuals using histograms, qq norm plots, and the Shapiro-Wilk's test,
      # we tested Equal Variance of Residuals using a fitted vs. residuals plot and the Breusch-Pagan Test,
      # we tested Simple Random Sample by thinkg about the data collection method.
#c) We then looked for significant associations (slopes/correlations)
    # 1) If the LINES conditions are met...
            # we ran a slope test and a Pearson's correlation test to see if there is a significant association
    # 2) If the LINES conditions are not met...
            # we ran a Mann-Kendall test (non-parametric test) to look for a significant correlation/tau 



#Function for selecting the Simple Linear Regressions based on meeting the equal variance and normal residuals 
#conditions and then returns the model and slope test results

#The tests for Linearity, Independence, and Simple Random Sampling are not included in the function. Linearity is tested in the 
#"Running the Simple Linear Regressions" section and the other two conditions should be tested by the analyst.

simple_linear_regressions <- function(population, size_variable, explanatory_var){ #input the population name and the size variable/response variable
  
  #assigning the population based on the inputted population
  if (population == "LM"){  #LM trees
    dataframe_metric = LM_fixed_field_data_processed_terrain_dist #assigning the LM tree/topography/distance dataframe
    dataframe_metric$Slope <- dataframe_metric$LM_slope_raster_15_data_pts #adding a slope column with the generic slope name to be able to call the same column across dataframes for different populations
  } else if (population == "LC"){ #LC trees
    dataframe_metric = LC_fixed_field_data_processed_terrain_dist #assigning the LC tree/topography/distance dataframe
    dataframe_metric$Slope <- dataframe_metric$LC_slope_raster_15_data_pts #adding a slope column with the generic slope name to be able to call the same column across dataframes for different populations
  } else if (population == "SD"){ #SD trees
    dataframe_metric = SD_fixed_field_data_processed_terrain_dist #assigning the SD tree/topography/distance dataframe
    dataframe_metric$Slope <- dataframe_metric$SD_slope_raster_15_data_pts #adding a slope column with the generic slope name to be able to call the same column across dataframes for different populations
  } else if (population == "All Points"){ #All trees across all populations
    dataframe_metric = all_points_fixed_field_data_processed_terrain #assigning the all points/population tree/topography/distance dataframe
    dataframe_metric$Slope <- dataframe_metric$all_points_slope_raster_15_data_pts #adding a slope column with the generic slope name to be able to call the same column across dataframes for different populations
  }
  
  # removing NAs
  dataframe_metric <- dataframe_metric %>%
    drop_na(Elevation..m.FIXED) #removing NAs in elevation 
  
  #assigning the size/response variable based on the user input
  if (size_variable == "SCA"){  #Short Canopy Axis
    size_variable_name = "Canopy_short" #storing the name of the size variable we are using
    size_metric = dataframe_metric$Canopy_short #assigning the short canopy axis variable to the size metric variable
  } else if (size_variable == "LCA"){ #Long Canopy Axis
    size_variable_name = "Canopy_long" #storing the name of the size variable we are using
    size_metric = dataframe_metric$Canopy_long #assigning the long canopy axis variable to the size metric variable
  } else if (size_variable == "CA"){ #Canopy Area
    size_variable_name = "Canopy_area" #storing the name of the size variable we are using
    size_metric = dataframe_metric$Canopy_area #assigning the canopy area variable to the size metric variable
  } else if (size_variable == "CS"){ #Crown Spread
    size_variable_name = "Crown_spread" #storing the name of the size variable we are using
    size_metric = dataframe_metric$Crown_spread #assigning the crown spread variable to the size metric variable
  } else if (size_variable == "DBH"){ #DBH
    size_variable_name = "DBH_ag" #storing the name of the size variable we are using
    size_metric = dataframe_metric$DBH_ag #assigning the DBH variable to the size metric variable
  }
  
  #assigning the explanatory variable based on the user input
  if (explanatory_var == "Elevation"){  #elevation
    explanatory_var_name = "Elevation..m.FIXED" #storing the name of the explanatory variable we are using
    explanatory_var_metric = dataframe_metric$Elevation..m.FIXED #assigning the elevation variable to the explanatory variable
  } else if (explanatory_var == "Slope"){ #Slope
    explanatory_var_name = "Slope" #storing the name of the explanatory variable we are using
    explanatory_var_metric = dataframe_metric$Slope #assigning the slope variable to the explanatory variable
  } else if (explanatory_var == "Distance to River"){ #Distance to River
    explanatory_var_name = "d" #storing the name of the explanatory variable we are using
    explanatory_var_metric = dataframe_metric$d #assigning the Distance to River variable to the explanatory variable
  } 
  
  #creating a dataframe with influential/outlier points removed 
  
  #using Cook's D to check for highly influential points that may skew the linear model results
  slr <- lm(size_metric ~ explanatory_var_metric) #creating a linear regression to use to calculate the Cook's D
  slr_cooks <- cooks.distance(slr) #creating a linear regression to use to calculate the Cook's D
  plot(slr_cooks, type = 'h') #checking to see which Cook's D are unusually high
  influential <- slr_cooks[(slr_cooks > (3 * mean(slr_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
  influential
  
  #removing points that were deemed too influential on the linear model fit
  dataframe_metric_no_outliers <- dataframe_metric[-c(as.numeric(names(influential))),]
  
  #creating the linear regressions
  
  #creating the base linear regression (no removal of outliers, no transformations)
 # slr_dist_base  <- lm(size_metric ~ explanatory_var_metric, data = dataframe_metric) #generating the linear regression 
  slr_dist_base  <- lm(as.formula(paste0(size_variable_name, " ~ ", explanatory_var_name)), data = dataframe_metric) #generating the linear regression) #generating the linear regression 
  
  #linear regression with transformations
  
  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_base_lg  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric) #generating the linear regression
  
  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_base_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric) #generating the linear regression 
  
  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_base_inv  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric) #generating the linear regression 
  
  #linear regression with transformations and removal of outliers
  
  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_no_out_lg  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers) #generating the linear regression 
  
  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_no_out_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers) #generating the linear regression 
  
  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_no_out_inv  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers) #generating the linear regression 
  
  #finding and storing the results of the best performing model
  
  #making a list of all of the potential models
  linear_regressions_models <- list("Base Model" = slr_dist_base, # base model slope test p-value
                                    "Base Model with Log Transformation" = slr_dist_base_lg, 
                                    "Base Model with Square Root Transformation" = slr_dist_base_sqrt, 
                                    "Base Model with Inverse Transformation" = slr_dist_base_inv, #slope test p-values with transformations of the response variable
                                    "No Outliers Model with Log Transformation" = slr_dist_no_out_lg, 
                                    "No Outliers Model with Square Root Transformation" =slr_dist_no_out_sqrt, 
                                    "No Outliers Model with Inverse Transformation" =slr_dist_no_out_inv)  #slope test p-values with no outliers and transformations of the response variable
  
  #creating an empty list to store Shapiro-Wilks Test
  shapiro_wilk_list <- c()
  
  #creating an empty list to store Breusch-Pagan Test
  breusch_pagan_test_list <- c()
  
  #assigning a variable to "unmet" until one model meets the conditions
  conditions <- "unmet"
  
  #iterating over each model to test if it meets the normality and equal variance conditions and stopping when a model meets the conditions
  for (k in 1:length(linear_regressions_models)){
    #assigning the correct model based on the iteration
    test_model <- linear_regressions_models[[k]] 
    
    #checking the conditions, specifically the normality and equal variance conditions, assuming the linearity, independence, and simple random sample conditions are met
    
    #checking normality with Shapiro-Wilk Test
    shapiro_wilk <- shapiro.test(test_model$residuals)
    
    #appending the Shapiro-Wilks results in a list
    shapiro_wilk_list <- c(shapiro_wilk_list, shapiro_wilk$p.value)
    
    #Equal variance tests (Bartlett's Test)
    
    # Getting the data used in the model for the Bartlett's Test
    model_data <- model.frame(test_model)
    
    # Getting the response and explanatory variable names for the Bartlett's Test
    response_variable <- names(model_data)[1] #all.vars(formula(chosen_model))[2]
    explanatory_variable <- names(model_data)[2] #all.vars(formula(chosen_model))[3]
    
    #Breusch-Pagan Test for equal variances when data is normal, if significant, residuals do not have equal variance
    breusch_pagan_test <- bptest(as.formula(paste(response_variable, " ~ ", explanatory_variable)), # This will cause an error,
                                 data = model_data)
    
    #appending the Breusch-Pagan Test results in a list
    breusch_pagan_test_list <- c(breusch_pagan_test_list, breusch_pagan_test$p.value)
    
    #return test results of the model that meets the conditions
    
    #if the normal residuals and equal variance of residuals conditions ARE met
    if (shapiro_wilk$p.value > 0.05 & breusch_pagan_test$p.value > 0.05) { 
      
      #printing whether the normality and equal variance conditions were met
      print(paste(names(linear_regressions_models)[k], "Data met the Normality and Equal Variance Conditions."))
      
      #Storing the results of the slope test from the summary of the linear model
      slope_test_results <- test_model$coefficients #slope test results
      
      #assigning the chosen model index to a variable
      chose_model_index <- k
      
      #assigning a variable to "met" when one model met the conditions
      conditions <- "met"
      
      break #exiting the loop when both conditions have been met
      
    } else if (shapiro_wilk$p.value < 0.05 & breusch_pagan_test$p.value > 0.05) {
      print(paste(names(linear_regressions_models)[k], "Data did not meet the Normality Condition but met the Equal Variance Condition.")) #printing whether the normality and equal variance conditions were met
      next #move on to testing if the next model meets both conditions
    } else if (shapiro_wilk$p.value > 0.05 & breusch_pagan_test$p.value < 0.05) {
      print(paste(names(linear_regressions_models)[k], "Data met the Normality Condition but not the Equal Variance Condition.")) #printing whether the normality and equal variance conditions were met
      next #move on to testing if the next model meets both conditions
    }
  }
  
  #if the conditions are not met (such as the residuals were not normal and the variance not equal) we use the model with the closest to normal residuals
  if (conditions != "met") {
    for (k in 1:length(linear_regressions_models)){
      
      #none of the models met the conditions for a slope test, so we stored the results of the model that met the slope test the clostest and the results of the non-parametric Mann-Kendall's Tau Test
      
      #we assign the test model be the one with the highest Shapriro-Wilks Test p-value, so the model with the residuals closest to normal
      test_model <- linear_regressions_models[[which.max(shapiro_wilk_list)]]
      
      #Storing the results of the slope test from the summary of the linear model
      slope_test_results <- test_model$coefficients #slope test results
      
      #assigning the chosen model index to a variable
      chose_model_index <- which.max(shapiro_wilk_list)
      
      #extracting the best model based on which one has the smallest p-value
      chosen_model <- test_model 
      
      #printing the chosen model
      chosen_model_statement <- print(paste("The chosen model is ", names(linear_regressions_models)[chose_model_index]))
      chosen_model_statement
      
      #because the conditions for the slope test were not met, we can use the non-parametric Mann-Kendall's test to check for significant associations
      
      #non-parametric Mann-Kendall Test, non-parametric test
      tau_result <- cor.test(explanatory_var_metric, size_metric,  method = "kendall")
      
      #printing that none of the models met the conditions
      lack_conditions_statement <- print(paste("None of the models met the normality and equal variance conditions, so we ran the non-parametric Mann-Kendall Test and stored the slope test results of the model that met the conditions the best"))
      lack_conditions_statement
      
    }
  } else { #if the condition is met
    
    #extracting the best model based on which one has the smallest p-value
    chosen_model <- test_model 
    
    #printing the chosen model
    chosen_model_statement <- print(paste("The chosen model is ", names(linear_regressions_models)[chose_model_index]))
    chosen_model_statement
    
    #printing that one of the models met the conditions
    lack_conditions_statement <- print(paste("Atleast one of the models met the normality and equal variance conditions"))
    lack_conditions_statement
    
    #if conditions are met, we do not need to run the non-parametric Mann-Kendall Test, so we can set those results to NA
    tau_result <- NA
    
  }
  
  
  #storing the model summary to extract values
  slope_test_summary <- summary(chosen_model) 
  
  #Storing the slope test t-value from the summary of the linear model
  slope_test_t_value <- slope_test_summary$coefficients[explanatory_var_name, "t value"]
  
  #Storing the slope test p-value from the summary of the linear model
  slope_test_p_value <- slope_test_summary$coefficients[explanatory_var_name, "Pr(>|t|)"]
  
  print(paste("Slope Test Results:", 
              "t =", slope_test_t_value, 
              "p =", slope_test_p_value))
  
  print(slope_test_summary)
  
  #returns the name of the model that performed the best, normality results, equal variance results, and the slope test results
  return(list("chosen_model" = chosen_model, #returning the chosen model
              "shapiro_wilk" = shapiro_wilk, #returning the Shapiro-Wilks Test for normal residual results
              "breusch_pagan_test" = breusch_pagan_test,#returning the Breusch-Pagan Test for equal variance of residuals results
              "slope_test_summary" = slope_test_summary, #returning the linear regression slope test results
              "slope_test_t_value" = slope_test_t_value, #returning the linear regression slope test results t-value
              "slope_test_p_value" = slope_test_p_value, #returning the linear regression slope test results p-value
              "chosen_model_statement" = chosen_model_statement, #returning a statement of which model was chosen
              "lack_conditions_statement" = lack_conditions_statement, #returning a statement if none of the models met the conditions
              "tau_result" = tau_result)) #returning the non-parametric Mann-Kendall Test results
  
}

#### Sizes vs. Elevation ####

# removing NAs
LM_fixed_field_data_processed_terrain_dist <- LM_fixed_field_data_processed_terrain_dist %>%
  drop_na(Elevation..m.FIXED) #removing NAs in elevation 

# LM

#SCA

#running the simple linear regression function
simple_linear_regressions_LM_SCA_elevation <- simple_linear_regressions("LM", "SCA", "Elevation")
simple_linear_regressions_LM_SCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_SCA_elevation$chosen_model, aes(x= simple_linear_regressions_LM_SCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_SCA_elevation$chosen_model, aes(sample = simple_linear_regressions_LM_SCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_SCA_elevation$chosen_model, aes(x = simple_linear_regressions_LM_SCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LM_SCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#LCA

#running the simple linear regression function
simple_linear_regressions_LM_LCA_elevation <- simple_linear_regressions("LM", "LCA", "Elevation")
simple_linear_regressions_LM_LCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_LCA_elevation$chosen_model, aes(x= simple_linear_regressions_LM_LCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_LCA_elevation$chosen_model, aes(sample = simple_linear_regressions_LM_LCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_LCA_elevation$chosen_model, aes(x = simple_linear_regressions_LM_LCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LM_LCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#CA

#running the simple linear regression function
simple_linear_regressions_LM_CA_elevation <- simple_linear_regressions("LM", "CA", "Elevation")
simple_linear_regressions_LM_CA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CA_elevation$chosen_model, aes(x= simple_linear_regressions_LM_CA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CA_elevation$chosen_model, aes(sample = simple_linear_regressions_LM_CA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CA_elevation$chosen_model, aes(x = simple_linear_regressions_LM_CA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LM_CA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

# Calculating the trend line for plotting
LM_trend_line_CA_elevation <- predict(loess(LM_fixed_field_data_processed_terrain_dist$Canopy_area ~ LM_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = (LM_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = LM_trend_line_CA_elevation), color = "red") +
  labs(x = "Elevation", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_LM_CS_elevation <- simple_linear_regressions("LM", "CS", "Elevation")
simple_linear_regressions_LM_CS_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CS_elevation$chosen_model, aes(x= simple_linear_regressions_LM_CS_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CS_elevation$chosen_model, aes(x= simple_linear_regressions_LM_CS_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CS_elevation$chosen_model, aes(x = simple_linear_regressions_LM_CS_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LM_CS_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#DBH

#running the simple linear regression function
simple_linear_regressions_LM_DBH_elevation <- simple_linear_regressions("LM", "DBH", "Elevation")
simple_linear_regressions_LM_DBH_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=Elevation..m.FIXED, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_LM_DBH_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_LM_DBH_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_DBH_elevation$chosen_model, aes(x = simple_linear_regressions_LM_DBH_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LM_DBH_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")


# LC


# removing NAs
LC_fixed_field_data_processed_terrain_dist <- LC_fixed_field_data_processed_terrain_dist %>%
  drop_na(Elevation..m.FIXED) #removing NAs in elevation 

#SCA

#running the simple linear regression function
simple_linear_regressions_LC_SCA_elevation <- simple_linear_regressions("LC", "SCA", "Elevation")
simple_linear_regressions_LC_SCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_SCA_elevation$chosen_model, aes(x= simple_linear_regressions_LC_SCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_SCA_elevation$chosen_model, aes(sample = simple_linear_regressions_LC_SCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_SCA_elevation$chosen_model, aes(x = simple_linear_regressions_LC_SCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LC_SCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#LCA

#running the simple linear regression function
simple_linear_regressions_LC_LCA_elevation <- simple_linear_regressions("LC", "LCA", "Elevation")
simple_linear_regressions_LC_LCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_LCA_elevation$chosen_model, aes(x= simple_linear_regressions_LC_LCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_LCA_elevation$chosen_model, aes(sample = simple_linear_regressions_LC_LCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_LCA_elevation$chosen_model, aes(x = simple_linear_regressions_LC_LCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LC_LCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#CA

#running the simple linear regression function
simple_linear_regressions_LC_CA_elevation <- simple_linear_regressions("LC", "CA", "Elevation")
simple_linear_regressions_LC_CA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CA_elevation$chosen_model, aes(x= simple_linear_regressions_LC_CA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CA_elevation$chosen_model, aes(sample = simple_linear_regressions_LC_CA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CA_elevation$chosen_model, aes(x = simple_linear_regressions_LC_CA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LC_CA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

# Calculating the trend line for plotting
LC_trend_line_CA_elevation <- predict(loess(LC_fixed_field_data_processed_terrain_dist$Canopy_area ~ LC_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = (LC_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = LC_trend_line_CA_elevation), color = "red") +
  labs(x = "Elevation", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_LC_CS_elevation <- simple_linear_regressions("LC", "CS", "Elevation")
simple_linear_regressions_LC_CS_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CS_elevation$chosen_model, aes(x= simple_linear_regressions_LC_CS_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CS_elevation$chosen_model, aes(x= simple_linear_regressions_LC_CS_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CS_elevation$chosen_model, aes(x = simple_linear_regressions_LC_CS_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LC_CS_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#DBH

#running the simple linear regression function
simple_linear_regressions_LC_DBH_elevation <- simple_linear_regressions("LC", "DBH", "Elevation")
simple_linear_regressions_LC_DBH_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=Elevation..m.FIXED, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_LC_DBH_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_LC_DBH_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_DBH_elevation$chosen_model, aes(x = simple_linear_regressions_LC_DBH_elevation$chosen_model$fitted.values, y = simple_linear_regressions_LC_DBH_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")


# SD

# removing NAs
SD_fixed_field_data_processed_terrain_dist <- SD_fixed_field_data_processed_terrain_dist %>%
  drop_na(Elevation..m.FIXED) #removing NAs in elevation 

#SCA

#running the simple linear regression function
simple_linear_regressions_SD_SCA_elevation <- simple_linear_regressions("SD", "SCA", "Elevation")
simple_linear_regressions_SD_SCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_SCA_elevation$chosen_model, aes(x= simple_linear_regressions_SD_SCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_SCA_elevation$chosen_model, aes(sample = simple_linear_regressions_SD_SCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_SCA_elevation$chosen_model, aes(x = simple_linear_regressions_SD_SCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_SD_SCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#LCA

#running the simple linear regression function
simple_linear_regressions_SD_LCA_elevation <- simple_linear_regressions("SD", "LCA", "Elevation")
simple_linear_regressions_SD_LCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_LCA_elevation$chosen_model, aes(x= simple_linear_regressions_SD_LCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_LCA_elevation$chosen_model, aes(sample = simple_linear_regressions_SD_LCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_LCA_elevation$chosen_model, aes(x = simple_linear_regressions_SD_LCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_SD_LCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#CA

#running the simple linear regression function
simple_linear_regressions_SD_CA_elevation <- simple_linear_regressions("SD", "CA", "Elevation")
simple_linear_regressions_SD_CA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CA_elevation$chosen_model, aes(x= simple_linear_regressions_SD_CA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CA_elevation$chosen_model, aes(sample = simple_linear_regressions_SD_CA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CA_elevation$chosen_model, aes(x = simple_linear_regressions_SD_CA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_SD_CA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

# Calculating the trend line for plotting
SD_trend_line_CA_elevation <- predict(loess(SD_fixed_field_data_processed_terrain_dist$Canopy_area ~ SD_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = (SD_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, y = SD_trend_line_CA_elevation), color = "red") +
  labs(x = "Elevation", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_SD_CS_elevation <- simple_linear_regressions("SD", "CS", "Elevation")
simple_linear_regressions_SD_CS_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CS_elevation$chosen_model, aes(x= simple_linear_regressions_SD_CS_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CS_elevation$chosen_model, aes(x= simple_linear_regressions_SD_CS_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CS_elevation$chosen_model, aes(x = simple_linear_regressions_SD_CS_elevation$chosen_model$fitted.values, y = simple_linear_regressions_SD_CS_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#DBH

#running the simple linear regression function
simple_linear_regressions_SD_DBH_elevation <- simple_linear_regressions("SD", "DBH", "Elevation")
simple_linear_regressions_SD_DBH_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=Elevation..m.FIXED, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_SD_DBH_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_SD_DBH_elevation$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_DBH_elevation$chosen_model, aes(x = simple_linear_regressions_SD_DBH_elevation$chosen_model$fitted.values, y = simple_linear_regressions_SD_DBH_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")


# All Points 

# removing NAs
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Elevation..m.FIXED) #removing NAs in elevation 

#SCA

#running the simple linear regression function
simple_linear_regressions_All_SCA_elevation <- simple_linear_regressions("All Points", "SCA", "Elevation")
simple_linear_regressions_All_SCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=Elevation..m.FIXED, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_SCA_elevation$chosen_model, aes(x= simple_linear_regressions_All_SCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_SCA_elevation$chosen_model, aes(sample = simple_linear_regressions_All_SCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_SCA_elevation$chosen_model, aes(x = simple_linear_regressions_All_SCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_All_SCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Elevation")

#LCA

#running the simple linear regression function
simple_linear_regressions_All_LCA_elevation <- simple_linear_regressions("All Points", "LCA", "Elevation")
simple_linear_regressions_All_LCA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=Elevation..m.FIXED, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_LCA_elevation$chosen_model, aes(x= simple_linear_regressions_All_LCA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_LCA_elevation$chosen_model, aes(sample = simple_linear_regressions_All_LCA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_LCA_elevation$chosen_model, aes(x = simple_linear_regressions_All_LCA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_All_LCA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Elevation")

#CA

#running the simple linear regression function
simple_linear_regressions_All_CA_elevation <- simple_linear_regressions("All Points", "CA", "Elevation")
simple_linear_regressions_All_CA_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=Elevation..m.FIXED, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Canopy Area (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_CA_elevation$chosen_model, aes(x= simple_linear_regressions_All_CA_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_CA_elevation$chosen_model, aes(sample = simple_linear_regressions_All_CA_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_CA_elevation$chosen_model, aes(x = simple_linear_regressions_All_CA_elevation$chosen_model$fitted.values, y = simple_linear_regressions_All_CA_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Elevation")

#CS

#running the simple linear regression function
simple_linear_regressions_All_CS_elevation <- simple_linear_regressions("All Points", "CS", "Elevation")
simple_linear_regressions_All_CS_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=Elevation..m.FIXED, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_CS_elevation$chosen_model, aes(x= simple_linear_regressions_All_CS_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_CS_elevation$chosen_model, aes(sample = simple_linear_regressions_All_CS_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_CS_elevation$chosen_model, aes(x = simple_linear_regressions_All_CS_elevation$chosen_model$fitted.values, y = simple_linear_regressions_All_CS_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Elevation")

#DBH

#running the simple linear regression function
simple_linear_regressions_All_DBH_elevation <- simple_linear_regressions("All Points", "DBH", "Elevation")
simple_linear_regressions_All_DBH_elevation

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=Elevation..m.FIXED, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Elevation (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_DBH_elevation$chosen_model, aes(x= simple_linear_regressions_All_DBH_elevation$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Elevation")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_DBH_elevation$chosen_model, aes(sample = simple_linear_regressions_All_DBH_elevation$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_DBH_elevation$chosen_model, aes(x = simple_linear_regressions_All_DBH_elevation$chosen_model$fitted.values, y = simple_linear_regressions_All_DBH_elevation$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Elevation")

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

# removing NAs
LM_fixed_field_data_processed_terrain_dist <- LM_fixed_field_data_processed_terrain_dist %>%
  drop_na(LM_slope_raster_15_data_pts) #removing NAs in elevation 

# LM

#SCA

#running the simple linear regression function
simple_linear_regressions_LM_SCA_slope <- simple_linear_regressions("LM", "SCA", "Slope")
simple_linear_regressions_LM_SCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square Root of Short Canopy Axis (º)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_SCA_slope$chosen_model, aes(x= simple_linear_regressions_LM_SCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_SCA_slope$chosen_model, aes(sample = simple_linear_regressions_LM_SCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_SCA_slope$chosen_model, aes(x = simple_linear_regressions_LM_SCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LM_SCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope")

#LCA

#running the simple linear regression function
simple_linear_regressions_LM_LCA_slope <- simple_linear_regressions("LM", "LCA", "Slope")
simple_linear_regressions_LM_LCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Long Canopy Axis")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_LCA_slope$chosen_model, aes(x= simple_linear_regressions_LM_LCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_LCA_slope$chosen_model, aes(sample = simple_linear_regressions_LM_LCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_LCA_slope$chosen_model, aes(x = simple_linear_regressions_LM_LCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LM_LCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope")

#CA

#running the simple linear regression function
simple_linear_regressions_LM_CA_slope <- simple_linear_regressions("LM", "CA", "Slope")
simple_linear_regressions_LM_CA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CA_slope$chosen_model, aes(x= simple_linear_regressions_LM_CA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CA_slope$chosen_model, aes(sample = simple_linear_regressions_LM_CA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CA_slope$chosen_model, aes(x = simple_linear_regressions_LM_CA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LM_CA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope")

# Calculating the trend line for plotting
LM_trend_line_CA_slope <- predict(loess(LM_fixed_field_data_processed_terrain_dist$Canopy_area ~ LM_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LM_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = (LM_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = LM_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = LM_trend_line_CA_slope), color = "red") +
  labs(x = "Slope", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_LM_CS_slope <- simple_linear_regressions("LM", "CS", "Slope")
simple_linear_regressions_LM_CS_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (m)")+
  ylab("Square root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CS_slope$chosen_model, aes(x= simple_linear_regressions_LM_CS_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CS_slope$chosen_model, aes(x= simple_linear_regressions_LM_CS_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CS_slope$chosen_model, aes(x = simple_linear_regressions_LM_CS_slope$chosen_model$fitted.values, y = simple_linear_regressions_LM_CS_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope")

#DBH

#running the simple linear regression function
simple_linear_regressions_LM_DBH_slope <- simple_linear_regressions("LM", "DBH", "Slope")
simple_linear_regressions_LM_DBH_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=LM_slope_raster_15_data_pts, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_DBH_slope$chosen_model, aes(x= simple_linear_regressions_LM_DBH_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_DBH_slope$chosen_model, aes(x= simple_linear_regressions_LM_DBH_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_DBH_slope$chosen_model, aes(x = simple_linear_regressions_LM_DBH_slope$chosen_model$fitted.values, y = simple_linear_regressions_LM_DBH_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope")


# LC

# removing NAs
LC_fixed_field_data_processed_terrain_dist <- LC_fixed_field_data_processed_terrain_dist %>%
  drop_na(LM_slope_raster_15_data_pts) #removing NAs in Slope 

#SCA

#running the simple linear regression function
simple_linear_regressions_LC_SCA_slope <- simple_linear_regressions("LC", "SCA", "Slope")
simple_linear_regressions_LC_SCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=log(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Log of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_SCA_slope$chosen_model, aes(x= simple_linear_regressions_LC_SCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Log of Short Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_SCA_slope$chosen_model, aes(sample = simple_linear_regressions_LC_SCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_SCA_slope$chosen_model, aes(x = simple_linear_regressions_LC_SCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LC_SCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope")

#LCA

#running the simple linear regression function
simple_linear_regressions_LC_LCA_slope <- simple_linear_regressions("LC", "LCA", "Slope")
simple_linear_regressions_LC_LCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=log(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Log of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_LCA_slope$chosen_model, aes(x= simple_linear_regressions_LC_LCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_LCA_slope$chosen_model, aes(sample = simple_linear_regressions_LC_LCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_LCA_slope$chosen_model, aes(x = simple_linear_regressions_LC_LCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LC_LCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope")

#CA

#running the simple linear regression function
simple_linear_regressions_LC_CA_slope <- simple_linear_regressions("LC", "CA", "Slope")
simple_linear_regressions_LC_CA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CA_slope$chosen_model, aes(x= simple_linear_regressions_LC_CA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CA_slope$chosen_model, aes(sample = simple_linear_regressions_LC_CA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CA_slope$chosen_model, aes(x = simple_linear_regressions_LC_CA_slope$chosen_model$fitted.values, y = simple_linear_regressions_LC_CA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope")

# Calculating the trend line for plotting
LC_trend_line_CA_slope <- predict(loess(LC_fixed_field_data_processed_terrain_dist$Canopy_area ~ LC_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = LC_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = (LC_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = LC_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = LC_trend_line_CA_slope), color = "red") +
  labs(x = "Slope", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_LC_CS_slope <- simple_linear_regressions("LC", "CS", "Slope")
simple_linear_regressions_LC_CS_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Log of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CS_slope$chosen_model, aes(x= simple_linear_regressions_LC_CS_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CS_slope$chosen_model, aes(x= simple_linear_regressions_LC_CS_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CS_slope$chosen_model, aes(x = simple_linear_regressions_LC_CS_slope$chosen_model$fitted.values, y = simple_linear_regressions_LC_CS_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope")

#DBH

#running the simple linear regression function
simple_linear_regressions_LC_DBH_slope <- simple_linear_regressions("LC", "DBH", "Slope")
simple_linear_regressions_LC_DBH_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=LM_slope_raster_15_data_pts, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_DBH_slope$chosen_model, aes(x= simple_linear_regressions_LC_DBH_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_DBH_slope$chosen_model, aes(x= simple_linear_regressions_LC_DBH_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_DBH_slope$chosen_model, aes(x = simple_linear_regressions_LC_DBH_slope$chosen_model$fitted.values, y = simple_linear_regressions_LC_DBH_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope")


# SD

# removing NAs
SD_fixed_field_data_processed_terrain_dist <- SD_fixed_field_data_processed_terrain_dist %>%
  drop_na(LM_slope_raster_15_data_pts) #removing NAs in Slope 

#SCA

#running the simple linear regression function
simple_linear_regressions_SD_SCA_slope <- simple_linear_regressions("SD", "SCA", "Slope")
simple_linear_regressions_SD_SCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_SCA_slope$chosen_model, aes(x= simple_linear_regressions_SD_SCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_SCA_slope$chosen_model, aes(sample = simple_linear_regressions_SD_SCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_SCA_slope$chosen_model, aes(x = simple_linear_regressions_SD_SCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_SD_SCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope")

#LCA

#running the simple linear regression function
simple_linear_regressions_SD_LCA_slope <- simple_linear_regressions("SD", "LCA", "Slope")
simple_linear_regressions_SD_LCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Log of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_LCA_slope$chosen_model, aes(x= simple_linear_regressions_SD_LCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_LCA_slope$chosen_model, aes(sample = simple_linear_regressions_SD_LCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_LCA_slope$chosen_model, aes(x = simple_linear_regressions_SD_LCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_SD_LCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope")

#CA

#running the simple linear regression function
simple_linear_regressions_SD_CA_slope <- simple_linear_regressions("SD", "CA", "Slope")
simple_linear_regressions_SD_CA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CA_slope$chosen_model, aes(x= simple_linear_regressions_SD_CA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CA_slope$chosen_model, aes(sample = simple_linear_regressions_SD_CA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CA_slope$chosen_model, aes(x = simple_linear_regressions_SD_CA_slope$chosen_model$fitted.values, y = simple_linear_regressions_SD_CA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope")

# Calculating the trend line for plotting
SD_trend_line_CA_slope <- predict(loess(SD_fixed_field_data_processed_terrain_dist$Canopy_area ~ SD_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts))

# Creating a trend line plot
ggplot() +
  geom_point(aes(x = SD_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = (SD_fixed_field_data_processed_terrain_dist$Canopy_area), color = "blue")) +
  geom_line(aes(x = SD_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts, y = SD_trend_line_CA_slope), color = "red") +
  labs(x = "Slope", y = "Canopy Area", title = "Trend Line Plot") +
  theme_minimal()

#CS

#running the simple linear regression function
simple_linear_regressions_SD_CS_slope <- simple_linear_regressions("SD", "CS", "Slope")
simple_linear_regressions_SD_CS_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_terrain_dist, (aes(x=LM_slope_raster_15_data_pts, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("Square root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CS_slope$chosen_model, aes(x= simple_linear_regressions_SD_CS_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CS_slope$chosen_model, aes(x= simple_linear_regressions_SD_CS_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CS_slope$chosen_model, aes(x = simple_linear_regressions_SD_CS_slope$chosen_model$fitted.values, y = simple_linear_regressions_SD_CS_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope")

#DBH

#running the simple linear regression function
simple_linear_regressions_SD_DBH_slope <- simple_linear_regressions("SD", "DBH", "Slope")
simple_linear_regressions_SD_DBH_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=LM_slope_raster_15_data_pts, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_DBH_slope$chosen_model, aes(x= simple_linear_regressions_SD_DBH_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_DBH_slope$chosen_model, aes(x= simple_linear_regressions_SD_DBH_slope$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_DBH_slope$chosen_model, aes(x = simple_linear_regressions_SD_DBH_slope$chosen_model$fitted.values, y = simple_linear_regressions_SD_DBH_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope")


#all points 

#removing NAs in SCA and slope from all points to run tests
all_points_fixed_field_data_processed_terrain <- all_points_fixed_field_data_processed_terrain %>%
  drop_na(Canopy_short) %>%
  drop_na(all_points_slope_raster_15_data_pts)

#SCA

#running the simple linear regression function
simple_linear_regressions_All_SCA_slope <- simple_linear_regressions("All Points", "SCA", "Slope")
simple_linear_regressions_All_SCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("SCA")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_SCA_slope$chosen_model, aes(x= simple_linear_regressions_All_SCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for SCA vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_SCA_slope$chosen_model, aes(sample = simple_linear_regressions_All_SCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_SCA_slope$chosen_model, aes(x = simple_linear_regressions_All_SCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_All_SCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Slope")

#LCA

#running the simple linear regression function
simple_linear_regressions_All_LCA_slope <- simple_linear_regressions("All Points", "LCA", "Slope")
simple_linear_regressions_All_LCA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("LCA")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_LCA_slope$chosen_model, aes(x= simple_linear_regressions_All_LCA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for LCA vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_LCA_slope$chosen_model, aes(sample = simple_linear_regressions_All_LCA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_LCA_slope$chosen_model, aes(x = simple_linear_regressions_All_LCA_slope$chosen_model$fitted.values, y = simple_linear_regressions_All_LCA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Slope")

#CA

#running the simple linear regression function
simple_linear_regressions_All_CA_slope <- simple_linear_regressions("All Points", "CA", "Slope")
simple_linear_regressions_All_CA_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("CA")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_CA_slope$chosen_model, aes(x= simple_linear_regressions_All_CA_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for CA vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_CA_slope$chosen_model, aes(sample = simple_linear_regressions_All_CA_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_CA_slope$chosen_model, aes(x = simple_linear_regressions_All_CA_slope$chosen_model$fitted.values, y = simple_linear_regressions_All_CA_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Slope")

#CS

#running the simple linear regression function
simple_linear_regressions_All_CS_slope <- simple_linear_regressions("All Points", "CS", "Slope")
simple_linear_regressions_All_CS_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("CS")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_CS_slope$chosen_model, aes(x= simple_linear_regressions_All_CS_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for CS vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_CS_slope$chosen_model, aes(sample = simple_linear_regressions_All_CS_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_CS_slope$chosen_model, aes(x = simple_linear_regressions_All_CS_slope$chosen_model$fitted.values, y = simple_linear_regressions_All_CS_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Slope")

#DBH

#running the simple linear regression function
simple_linear_regressions_All_DBH_slope <- simple_linear_regressions("All Points", "DBH", "Slope")
simple_linear_regressions_All_DBH_slope

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = all_points_fixed_field_data_processed_terrain, (aes(x=all_points_slope_raster_15_data_pts, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Slope (º)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_All_DBH_slope$chosen_model, aes(x= simple_linear_regressions_All_DBH_slope$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Slope")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_All_DBH_slope$chosen_model, aes(sample = simple_linear_regressions_All_DBH_slope$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_All_DBH_slope$chosen_model, aes(x = simple_linear_regressions_All_DBH_slope$chosen_model$fitted.values, y = simple_linear_regressions_All_DBH_slope$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Slope")

#### Sizes vs. Distance to River ####

#LM

#SCA

#running the simple linear regression function
simple_linear_regressions_LM_SCA_distance <- simple_linear_regressions("LM", "SCA", "d")
simple_linear_regressions_LM_SCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_SCA_distance$chosen_model, aes(x= simple_linear_regressions_LM_SCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_SCA_distance$chosen_model, aes(sample = simple_linear_regressions_LM_SCA_distance$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_SCA_distance$chosen_model, aes(x = simple_linear_regressions_LM_SCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LM_SCA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LM_LCA_distance <- simple_linear_regressions("LM", "LCA", "d")
simple_linear_regressions_LM_LCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_LCA_distance$chosen_model, aes(x= simple_linear_regressions_LM_LCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_LCA_distance$chosen_model, aes(sample = simple_linear_regressions_LM_LCA_distance$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_LCA_distance$chosen_model, aes(x = simple_linear_regressions_LM_LCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LM_LCA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LM_CA_distance <- simple_linear_regressions("LM", "CA", "d")
simple_linear_regressions_LM_CA_distance 

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CA_distance$chosen_model, aes(x= simple_linear_regressions_LM_CA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CA_distance$chosen_model, aes(sample = simple_linear_regressions_LM_CA_distance$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CA_distance$chosen_model, aes(x = simple_linear_regressions_LM_CA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LM_CA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LM_CS_distance <- simple_linear_regressions("LM", "CS", "d")
simple_linear_regressions_LM_CS_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CS_distance$chosen_model, aes(x= simple_linear_regressions_LM_CS_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CS_distance$chosen_model, aes(x= simple_linear_regressions_LM_CS_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CS_distance$chosen_model, aes(x = simple_linear_regressions_LM_CS_distance$chosen_model$fitted.values, y = simple_linear_regressions_LM_CS_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LM_DBH_distance <- simple_linear_regressions("LM", "DBH", "d")
simple_linear_regressions_LM_DBH_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_DBH_distance$chosen_model, aes(x= simple_linear_regressions_LM_DBH_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_DBH_distance$chosen_model, aes(x= simple_linear_regressions_LM_DBH_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_DBH_distance$chosen_model, aes(x = simple_linear_regressions_LM_DBH_distance$chosen_model$fitted.values, y = simple_linear_regressions_LM_DBH_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance")

#LC

#SCA

#running the simple linear regression function
simple_linear_regressions_LC_SCA_distance <- simple_linear_regressions("LC", "SCA", "d")
simple_linear_regressions_LC_SCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=log(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_SCA_distance$chosen_model, aes(x= simple_linear_regressions_LC_SCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_SCA_distance$chosen_model, aes(x= simple_linear_regressions_LC_SCA_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_SCA_distance$chosen_model, aes(x = simple_linear_regressions_LC_SCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LC_SCA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LC_LCA_distance <- simple_linear_regressions("LC", "LCA", "d")
simple_linear_regressions_LC_LCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=log(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_LCA_distance$chosen_model, aes(x= simple_linear_regressions_LC_LCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_LCA_distance$chosen_model, aes(x= simple_linear_regressions_LC_LCA_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_LCA_distance$chosen_model, aes(x = simple_linear_regressions_LC_LCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LC_LCA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LC_CA_distance <- simple_linear_regressions("LC", "CA", "d")
simple_linear_regressions_LC_CA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=log(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CA_distance$chosen_model, aes(x= simple_linear_regressions_LC_CA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CA_distance$chosen_model, aes(x= simple_linear_regressions_LC_CA_distance$chosen_model$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CA_distance$chosen_model, aes(x = simple_linear_regressions_LC_CA_distance$chosen_model$fitted.values, y = simple_linear_regressions_LC_CA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LC_CS_distance <- simple_linear_regressions("LC", "CS", "d")
simple_linear_regressions_LC_CS_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=log(Crown_spread))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CS_distance$chosen_model, aes(x= simple_linear_regressions_LC_CS_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CS_distance$chosen_model, aes(x= simple_linear_regressions_LC_CS_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CS_distance$chosen_model, aes(x = simple_linear_regressions_LC_CS_distance$chosen_model$fitted.values, y = simple_linear_regressions_LC_CS_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LC_DBH_distance <- simple_linear_regressions("LC", "DBH", "d")
simple_linear_regressions_LC_DBH_distance 

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_DBH_distance$chosen_model, aes(x= simple_linear_regressions_LC_DBH_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_DBH_distance$chosen_model, aes(x= simple_linear_regressions_LC_DBH_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_DBH_distance$chosen_model, aes(x = simple_linear_regressions_LC_DBH_distance$chosen_model$fitted.values, y = simple_linear_regressions_LC_DBH_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance")

#SD

#SCA

#running the simple linear regression function
simple_linear_regressions_SD_SCA_distance <- simple_linear_regressions("SD", "SCA", "d")
simple_linear_regressions_SD_SCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_SCA_distance$chosen_model, aes(x= simple_linear_regressions_SD_SCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_SCA_distance$chosen_model, aes(x= simple_linear_regressions_SD_SCA_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_SCA_distance$chosen_model, aes(x = simple_linear_regressions_SD_SCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_SD_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_SD_LCA_distance <- simple_linear_regressions("SD", "LCA", "d")
simple_linear_regressions_SD_LCA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=log(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_LCA_distance$chosen_model, aes(x= simple_linear_regressions_SD_LCA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_LCA_distance$chosen_model, aes(x= simple_linear_regressions_SD_LCA_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_LCA_distance$chosen_model, aes(x = simple_linear_regressions_SD_LCA_distance$chosen_model$fitted.values, y = simple_linear_regressions_SD_LCA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_SD_CA_distance <- simple_linear_regressions("SD", "CA", "d")
simple_linear_regressions_SD_CA_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CA_distance$chosen_model, aes(x= simple_linear_regressions_SD_CA_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CA_distance$chosen_model, aes(x= simple_linear_regressions_SD_CA_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CA_distance$chosen_model, aes(x = simple_linear_regressions_SD_CA_distance$chosen_model$fitted.values, y = simple_linear_regressions_SD_CA_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_SD_CS_distance <- simple_linear_regressions("SD", "CS", "d")
simple_linear_regressions_SD_CS_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Crown_spread))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CS_distance$chosen_model, aes(x= simple_linear_regressions_SD_CS_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CS_distance$chosen_model, aes(x= simple_linear_regressions_SD_CS_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CS_distance$chosen_model, aes(x = simple_linear_regressions_SD_CS_distance$chosen_model$fitted.values, y = simple_linear_regressions_SD_CS_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_SD_DBH_distance <- simple_linear_regressions("SD", "DBH", "d")
simple_linear_regressions_SD_DBH_distance

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_DBH_distance$chosen_model, aes(x= simple_linear_regressions_SD_DBH_distance$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance (m)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_DBH_distance$chosen_model, aes(x= simple_linear_regressions_SD_DBH_distance$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_DBH_distance$chosen_model, aes(x = simple_linear_regressions_SD_DBH_distance$chosen_model$fitted.values, y = simple_linear_regressions_SD_DBH_distance$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance (m)")


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






