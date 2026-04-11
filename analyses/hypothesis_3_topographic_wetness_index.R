# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei's size/shape is influenced by either elevation, slope, aspect, distance to rivers, and/or topographic wetness index %%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

#the purpose of the script is to determine how the elevation, slope, and distance to river each 
#relate to size/shape of the trees using single linear regressions and how aspect relates to the size/shape
#of the trees using difference in means tests/ANOVAs/Kruskal-Wallis Tests

# A significant relationship could indicate how elevation, slope, aspect, distance to river, and topographic wetness index may influence the size and shape of the 
# trees and potentially help explain their distribution. 

#The script is broken up into these sections:

# 1) Loading and processing the packages and processed data for the trees, topography, soil metrics, and distance to river in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
# Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
# Extracting and processing slope, elevation, and aspect (4 and 8 cardinal directions) data using 15 m res rasters,
# Extracting the topographic wetness index
# Extracting the distance to the river of each tree for each population,
# Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
# 2) Comparing size of the trees to their elevation to look for any relationships using single linear regression
# 3) Comparing size of the trees to their slope to look for any relationships using single linear regression
# 4) Comparing size of the trees to their distance to a river to look for any relationships using single linear regression
# 5) Comparing size of the trees to their aspect to look for a relationship using ANOVA / Kruskal-Wallis Tests
# 6) Comparing size of the trees to their topographic wetness index to look for any relationships using single linear regression


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
library(lmtest) #to be able to run the Breusch-Pagan Test

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

population = "LC"
size_variable = "SCA"
explanatory_var = "d"



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
  
  #removing NAs 
  dataframe_metric <- dataframe_metric %>%
    filter(!is.na(Elevation..m.FIXED)) %>%
    filter(!is.na(DBH_ag)) %>%
    filter(!is.na(Canopy_short)) %>%
    filter(!is.na(Canopy_long)) %>%
    filter(!is.na(Crown_spread)) %>%
    filter(!is.na(Canopy_area)) 
  
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
  slr_dist_base_lg  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric, subset = is.finite(dataframe_metric[[transformed_variable]])) #generating the linear regression
  
  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_base_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric, subset = is.finite(dataframe_metric[[transformed_variable]])) #generating the linear regression 
  
  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_base_inv  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric, subset = is.finite(dataframe_metric[[transformed_variable]])) #generating the linear regression 
  
  #linear regression with transformations and removal of outliers
  
  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_no_out_lg  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers, subset = is.finite(dataframe_metric_no_outliers[[transformed_variable]])) #generating the linear regression 
  
  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_no_out_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers, subset = is.finite(dataframe_metric_no_outliers[[transformed_variable]])) #generating the linear regression 
  
  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_no_out_inv  <- lm(as.formula(paste0(transformed_variable, " ~ ", explanatory_var_name)), data = dataframe_metric_no_outliers, subset = is.finite(dataframe_metric_no_outliers[[transformed_variable]])) #generating the linear regression 
  
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

