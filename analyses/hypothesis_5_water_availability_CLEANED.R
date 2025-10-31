# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Investigating to see if Quercus brandegeei tree shapes are influenced by Water Availability Proxies%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to determine if the trees' distance to river (a proxy for water availability)
# has a significant relationship with the size/shape of trees (short canopy axis, long canopy axis, canopy area, crown spread, DBH)
# A significant relationship could indicate how distance to river may influence the size and shape of the 
# trees and potentially help explain their distribution. 

# To test this, we used Simple Linear Regressions for each population, using distance to river as the explanatory 
# variables and each size metric as the response variable. We used a slope test to check for a significant relationship
# between distance the size/shape of trees.

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
        # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
                #generating river and point buffers and bounding boxes,
        # Extracting the distance to the river of each tree for each population,
# 3) Making a function for creating the single linear regressions, finding the best models (with transformations and removal of outliers or not) 
# that meet the conditions for SLRs
# 4) Running the function and storing the outputs

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
library(lmtest) #to use the Breuch-Pagan Test

# loading in the processed tree data 
source("./analyses/Data_Processing_Script.R")

#### Creating the Simple Linear Regression function ####

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

simple_linear_regressions <- function(population, size_variable){ #input the population name and the size variable/response variable
  
  #assigning the population based on the inputted population
  if (population == "LM"){  #LM trees
    dataframe_distance = LM_fixed_field_data_processed_distance #assigning the LM tree/distance dataframe
  } else if (population == "LC"){ #LC trees
    dataframe_distance = LC_fixed_field_data_processed_distance #assigning the LC tree/distance dataframe
  } else if (population == "SD"){ #SD trees
    dataframe_distance = SD_fixed_field_data_processed_distance #assigning the SD tree/distance dataframe
  }
  
  #assigning the size/response variable based on the user input
  if (size_variable == "SCA"){  #Short Canopy Axis
    size_variable_name = "Canopy_short" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_short #assigning the short canopy axis variable to the size metric variable
  } else if (size_variable == "LCA"){ #Long Canopy Axis
    size_variable_name = "Canopy_long" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_long #assigning the long canopy axis variable to the size metric variable
  } else if (size_variable == "CA"){ #Canopy Area
    size_variable_name = "Canopy_area" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Canopy_area #assigning the canopy area variable to the size metric variable
  } else if (size_variable == "CS"){ #Crown Spread
    size_variable_name = "Crown_spread" #storing the name of the size variable we are using
    size_metric = dataframe_distance$Crown_spread #assigning the crown spread variable to the size metric variable
  } else if (size_variable == "DBH"){ #DBH
    size_variable_name = "DBH_ag" #storing the name of the size variable we are using
    size_metric = dataframe_distance$DBH_ag #assigning the DBH variable to the size metric variable
  }
  
  #creating a dataframe with influential/outlier points removed 
  
  #using Cook's D to check for highly influential points that may skew the linear model results
  slr <- lm(size_metric ~ dataframe_distance$d) #creating a linear regression to use to calculate the Cook's D
  slr_cooks <- cooks.distance(slr) #creating a linear regression to use to calculate the Cook's D
  plot(slr_cooks, type = 'h') #checking to see which Cook's D are unusually high
  influential <- slr_cooks[(slr_cooks > (3 * mean(slr_cooks, na.rm = TRUE)))] #remove points with Cook's D that are bigger than 3 times the mean Cook's D (the influential points)
  influential
  
  #removing points that were deemed too influential on the linear model fit
  dataframe_distance_no_outliers <- dataframe_distance[-c(as.numeric(names(influential))),]
  
  #creating the linear regressions
  
  #creating the base linear regression (no removal of outliers, no transformations)
  slr_dist_base  <- lm(size_metric ~ dataframe_distance$d) #generating the linear regression 

  #linear regression with transformations

  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_base_lg  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression

  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_base_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression 

  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_base_inv  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance) #generating the linear regression 

  #linear regression with transformations and removal of outliers
  
  #linear regression with a log transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_lg") #storing the name of the transformed variable
  slr_dist_no_out_lg  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 

  #linear regression with a square root transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_sqrt") #storing the name of the transformed variable
  slr_dist_no_out_sqrt  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 

  #linear regression with a inverse transformation of the response variable
  transformed_variable <- paste0(size_variable_name, "_inv") #storing the name of the transformed variable
  slr_dist_no_out_inv  <- lm(as.formula(paste0(transformed_variable, " ~ d")), data = dataframe_distance_no_outliers) #generating the linear regression 
 
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
    response_var <- names(model_data)[1] #all.vars(formula(chosen_model))[2]
    explanatory_var <- names(model_data)[2] #all.vars(formula(chosen_model))[3]
    
    #Breusch-Pagan Test for equal variances when data is normal, if significant, residuals do not have equal variance
    breusch_pagan_test <- bptest(as.formula(paste(response_var, " ~ ", explanatory_var)), # This will cause an error,
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
      
      #printing that none of the models met the conditions
      lack_conditions_statement <- print(paste("None of the models met the normality and equal variance conditions, so we used the one that met the conditions the best"))
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
    
  }
  
  #storing the model summary to extract values
  slope_test_summary <- summary(chosen_model) 
  
  #Storing the slope test t-value from the summary of the linear model
  slope_test_t_value <- slope_test_summary$coefficients["d", "t value"]
  
  #Storing the slope test p-value from the summary of the linear model
  slope_test_p_value <- slope_test_summary$coefficients["d", "Pr(>|t|)"]
  
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
              "lack_conditions_statement" = lack_conditions_statement)) #returning a statement if none of the models met the conditions
  
}


#### Running the Simple Linear Regressions ####

#LM

#SCA

#running the simple linear regression function
simple_linear_regressions_LM_SCA <- simple_linear_regressions("LM", "SCA")
simple_linear_regressions_LM_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_short))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_SCA$chosen_model, aes(x= simple_linear_regressions_LM_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Square Root of Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_SCA$chosen_model, aes(sample = simple_linear_regressions_LM_SCA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_SCA$chosen_model, aes(x = simple_linear_regressions_LM_SCA$chosen_model$fitted.values, y = simple_linear_regressions_LM_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LM_LCA <- simple_linear_regressions("LM", "LCA")
simple_linear_regressions_LM_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_long))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("Square Root of Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_LCA$chosen_model, aes(x= simple_linear_regressions_LM_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_LCA$chosen_model, aes(sample = simple_linear_regressions_LM_LCA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_LCA$chosen_model, aes(x = simple_linear_regressions_LM_LCA$chosen_model$fitted.values, y = simple_linear_regressions_LM_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LM_CA <- simple_linear_regressions("LM", "CA")
simple_linear_regressions_LM_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=sqrt(Canopy_area))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Square Root of Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CA$chosen_model, aes(x= simple_linear_regressions_LM_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CA$chosen_model, aes(sample = simple_linear_regressions_LM_CA$chosen_model$residual))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CA$chosen_model, aes(x = simple_linear_regressions_LM_CA$chosen_model$fitted.values, y = simple_linear_regressions_LM_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LM_CS <- simple_linear_regressions("LM", "CS")
simple_linear_regressions_LM_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_CS$chosen_model, aes(x= simple_linear_regressions_LM_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_CS$chosen_model, aes(x= simple_linear_regressions_LM_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_CS$chosen_model, aes(x = simple_linear_regressions_LM_CS$chosen_model$fitted.values, y = simple_linear_regressions_LM_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LM_DBH <- simple_linear_regressions("LM", "DBH")
simple_linear_regressions_LM_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LM_fixed_field_data_processed_distance, (aes(x=d, y=log(DBH_ag))))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LM_DBH$chosen_model, aes(x= simple_linear_regressions_LM_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LM_DBH$chosen_model, aes(x= simple_linear_regressions_LM_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LM_DBH$chosen_model, aes(x = simple_linear_regressions_LM_DBH$chosen_model$fitted.values, y = simple_linear_regressions_LM_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#LC

#SCA

#running the simple linear regression function
simple_linear_regressions_LC_SCA <- simple_linear_regressions("LC", "SCA")
simple_linear_regressions_LC_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_SCA$chosen_model, aes(x= simple_linear_regressions_LC_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_SCA$chosen_model, aes(x= simple_linear_regressions_LC_SCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_SCA$chosen_model, aes(x = simple_linear_regressions_LC_SCA$chosen_model$fitted.values, y = simple_linear_regressions_LC_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_LC_LCA <- simple_linear_regressions("LC", "LCA")
simple_linear_regressions_LC_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_LCA$chosen_model, aes(x= simple_linear_regressions_LC_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_LCA$chosen_model, aes(x= simple_linear_regressions_LC_LCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_LCA$chosen_model, aes(x = simple_linear_regressions_LC_LCA$chosen_model$fitted.values, y = simple_linear_regressions_LC_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Inverse Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_LC_CA <- simple_linear_regressions("LC", "CA")
simple_linear_regressions_LC_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CA$chosen_model, aes(x= simple_linear_regressions_LC_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CA$chosen_model, aes(x= simple_linear_regressions_LC_CA$chosen_model$residuals))+
  geom_qq()

#Shapiro-Wilk Test
shapiro.test(LC_slr_dist_ca$residuals) #significantly not normal, except when outliers are removed

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CA$chosen_model, aes(x = simple_linear_regressions_LC_CA$chosen_model$fitted.values, y = simple_linear_regressions_LC_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_LC_CS <- simple_linear_regressions("LC", "CS")
simple_linear_regressions_LC_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_CS$chosen_model, aes(x= simple_linear_regressions_LC_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_CS$chosen_model, aes(x= simple_linear_regressions_LC_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_CS$chosen_model, aes(x = simple_linear_regressions_LC_CS$chosen_model$fitted.values, y = simple_linear_regressions_LC_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_LC_DBH <- simple_linear_regressions("LC", "DBH")
simple_linear_regressions_LC_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = LC_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_LC_DBH$chosen_model, aes(x= simple_linear_regressions_LC_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_LC_DBH$chosen_model, aes(x= simple_linear_regressions_LC_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_LC_DBH$chosen_model, aes(x = simple_linear_regressions_LC_DBH$chosen_model$fitted.values, y = simple_linear_regressions_LC_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Inverse Distance")

#SD

#SCA

#running the simple linear regression function
simple_linear_regressions_SD_SCA <- simple_linear_regressions("SD", "SCA")
simple_linear_regressions_SD_SCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Short Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_SCA$chosen_model, aes(x= simple_linear_regressions_SD_SCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_SCA$chosen_model, aes(x= simple_linear_regressions_SD_SCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_SCA$chosen_model, aes(x = simple_linear_regressions_SD_SCA$chosen_model$fitted.values, y = simple_linear_regressions_SD_SCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for SCA and Inverse Distance")

#LCA

#running the simple linear regression function
simple_linear_regressions_SD_LCA <- simple_linear_regressions("SD", "LCA")
simple_linear_regressions_SD_LCA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_long)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Long Canopy Axis (m)")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_LCA$chosen_model, aes(x= simple_linear_regressions_SD_LCA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Long Canopy Axis vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_LCA$chosen_model, aes(x= simple_linear_regressions_SD_LCA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_LCA$chosen_model, aes(x = simple_linear_regressions_SD_LCA$chosen_model$fitted.values, y = simple_linear_regressions_SD_LCA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for LCA and Distance")

#CA

#running the simple linear regression function
simple_linear_regressions_SD_CA <- simple_linear_regressions("SD", "CA")
simple_linear_regressions_SD_CA

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Canopy_area)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Canopy Area")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CA$chosen_model, aes(x= simple_linear_regressions_SD_CA$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Canopy Area vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CA$chosen_model, aes(x= simple_linear_regressions_SD_CA$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CA$chosen_model, aes(x = simple_linear_regressions_SD_CA$chosen_model$fitted.values, y = simple_linear_regressions_SD_CA$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CA and Inverse Distance")

#CS

#running the simple linear regression function
simple_linear_regressions_SD_CS <- simple_linear_regressions("SD", "CS")
simple_linear_regressions_SD_CS

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=Crown_spread)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Inverse Distance (m)")+
  ylab("Crown Spread")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_CS$chosen_model, aes(x= simple_linear_regressions_SD_CS$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Crown Spread vs. Inverse Distance")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_CS$chosen_model, aes(x= simple_linear_regressions_SD_CS$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_CS$chosen_model, aes(x = simple_linear_regressions_SD_CS$chosen_model$fitted.values, y = simple_linear_regressions_SD_CS$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for CS and Inverse Distance")

#DBH

#running the simple linear regression function
simple_linear_regressions_SD_DBH <- simple_linear_regressions("SD", "DBH")
simple_linear_regressions_SD_DBH

#checking linearity 

#plotting the scatterplot and linear model in ggplot
ggplot(data = SD_fixed_field_data_processed_distance, (aes(x=d, y=DBH_ag)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Distance (m)")+
  ylab("DBH")

#looking at the normality of residuals with a histogram and qqnorm plot

#histogram
ggplot(simple_linear_regressions_SD_DBH$chosen_model, aes(x= simple_linear_regressions_SD_DBH$chosen_model$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for DBH vs. Distance (m)")+
  xlab("Residuals")+
  ylab("Frequency")

#qqnorm plot
ggplot(simple_linear_regressions_SD_DBH$chosen_model, aes(x= simple_linear_regressions_SD_DBH$chosen_model$residuals))+
  geom_qq()

#looking at equal variance of residuals with a residuals vs. fitted values plot with a residuals vs. fitted values plot
ggplot(data = simple_linear_regressions_SD_DBH$chosen_model, aes(x = simple_linear_regressions_SD_DBH$chosen_model$fitted.values, y = simple_linear_regressions_SD_DBH$chosen_model$residual))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for DBH and Distance (m)")


