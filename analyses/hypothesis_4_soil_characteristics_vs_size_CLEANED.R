# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Quercus brandegeei tree shape are Influenced by Their Soil Metrics%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to determine whether there is a significant relationship between the 
# size (short canopy axis, DBH, etc.) of Q. brandegeei trees at each population and their soil
# characteristics (clay content, pH, etc.).
# A significant relationship could indicate which soil characteristics may influence the size and shape of the 
# trees and potentially help explain their distribution. 

# To test this, we used Simple Linear Regressions for each population, using the soil metrics as the explanatory 
# variables and each size metric as the response variable. We used a permutation test, comparing the slope of the 
# real linear regressions to a distribution of slopes for the linear regression if the size and shape metrics 
# were randomly shuffled, to support if the relationship is significant without assumption a normal distribution.

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
          # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
          # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
#boxs around the  rivers, stacking the rasters for each population, and processing them into one dataframe for all and each population
          # Creating the four new soil metrics: Sandy available Water (0-5 and 100-200 cm) and Clay/Loam Available Water (0-5 and 100-200 cm)
# 2) Making a function for calculating the real and randomized single linear regressions/slopes comparing tree 
# size to soil metrics for each population and calculating the p-values. 
# 3) Running the function and storing the output, graphing the histograms/slopes, turning the results into a dataframe
# 4) Using tables and heat maps to summarize the results (the slopes and their significances)

#### Loading libraries and Relevant Data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) #ggplot extention
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

# loading in the processed tree data 
source("./analyses/Data_Processing_Script.R")

#### Creating the Function Comparing the soil vs. size values ####

# This function that generates real and randomly generated slopes from single linear regressions comparing tree 
# sizes to soil metrics for each population to see if there is a significant relationship between tree size and soil characteristics. 

# The function components are 
    # 1) generating an empty list for the slopes and the p-values from comparing the real and randomized slopes. 
    # 2) a loop to iterate over each population, over each soil metric, and over each size characteristics 
          # in which we generate the list of randomized slopes (1,000 permutations, shuffling the size metrics 
              # across the trees) and real slope, 
          # plot the histogram of randomized slopes and real slope,
          # generate the p-values (using significantly different),
          # and store the real slopes and p-values in their respective arrays

# The function returns the slope array, p-value array, and stored histogram plot list.

slopes_simulations <- function () {
  
  # Creating the slope and p-value arrays 
  
  Populations <- c("LM", "LC", "SD") #creating a list of population names
  Size.metrics <- c("SCA", "LCA", "CA", "CS", "DBH") #creating a list of size metric names
  Soil.metrics <- c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                     "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                     "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                    "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                    "Volume of water content -1500 kpa 100-200", 
                    "Nitrogen 0-5", "Nitrogen 100-200", 
                    "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                    "Sand Available Water 0-5", "Sand Available Water 100-200",
                    "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200") #creating a list of soil metric names
  #creating the slopes and p-value arrays [soil metrics x size metrics x populations]
  slopes_array <- array(NA, dim = c(length(Soil.metrics), length(Size.metrics), length(Populations)), 
                        list(Soil.metrics, Size.metrics, Populations)) #real slopes array
  p_values_array <- array(NA, dim = c(length(Soil.metrics), length(Size.metrics), length(Populations)), 
                        list(Soil.metrics, Size.metrics, Populations)) # p-value array

  #to store generated histogram plots
  plot_list <- list()   
  
  #starting a loop iterating over populations
  for (i in 1:length(Populations)) { 
   
     #assigning the population based on the current tree
     if (Populations[i] == "LM"){  #LM trees
       Population = "LM"
       dataframe_soils = LM_fixed_field_data_processed_soils #assigning the LM tree/soils dataframe
     } else if (Populations[i] == "LC"){ #LC trees
       Population = "LC"
       dataframe_soils = LC_fixed_field_data_processed_soils #assigning the LC tree/soils dataframe
     } else if (Populations[i] == "SD"){ #SD trees
       Population = "SD"
       dataframe_soils = SD_fixed_field_data_processed_soils #assigning the SD tree/soils dataframe
     }

    #creating a dataframe of just the soil characteristics
    fixed_field_data_processed_soils.condensed <-  st_drop_geometry(dataframe_soils) #dropping the spatial geometry to easily remove columns
    fixed_field_data_processed_soils.condensed <- fixed_field_data_processed_soils.condensed[, c(c("clay.content.0.5", "clay.content.100.200", "silt.0.5", "silt.100.200", "sand.0.5", "sand.100.200",
                                                                                                   "ph_0.5", "ph_100.200",  "vol_water_.10_0.5",
                                                                                                   "vol_water_.10_100.200", "vol_water_0.5",
                                                                                                   "vol_water_100.200", "vol_water_.1500kPa_0.5", 
                                                                                                   "vol_water_.1500_100.200", 
                                                                                                   "nitrogen.0.5", "nitrogen.100.200", 
                                                                                                   "SOC.0.5", "SOC.100.200",
                                                                                                   "sandy_avail_water_0.5", "sandy_avail_water_100.200",
                                                                                                   "clay_loam_avail_water_0.5", "clay_loam_avail_water_100.200"))] #dropping all columns except the relevant soil metrics
   
    
    #creating a dataframe of just the size variables
    fixed_field_data_processed_size_variables <- st_drop_geometry(dataframe_soils) #dropping the spatial geometry to easily remove columns
    fixed_field_data_processed_size_variables <- fixed_field_data_processed_size_variables[, c("Canopy_short", "Canopy_long",
                                                                                               "Canopy_area","Crown_spread", 
                                                                                               "DBH_ag")] #dropping all columns except the relevant size metrics
    
    #iterating over the number of soil variables
    for (j in 1:(ncol(fixed_field_data_processed_soils.condensed))){ 
      
      #iterating over the number of soil variables
      for (k in 1:(ncol(fixed_field_data_processed_size_variables))) {
        
        #extracting slopes from comparing soil values with randomized shape/size values with linear regressions
        
        #creating empty list to collect randomized slope values
        slopes <- c() 
        
        #creating the list of randomized slopes
        set.seed(21) #setting the seed
        for (y in 1:1000){ #for 1000 permutations
          fixed_field_data_processed_soils_shuffled <- transform(dataframe_soils, variable.shuffled = sample(fixed_field_data_processed_size_variables[,k])) #create a data frame with shuffled soil metrics across trees for the respective population
          lm <- lm(fixed_field_data_processed_soils_shuffled$variable.shuffled ~ fixed_field_data_processed_soils.condensed[,j]) #creating the single linear regression of the shuffle size metric vs. soil metrics for the population's trees
          lm_sum <- summary(lm) #extracting the linear regression information
          slopes <- c(slopes, lm_sum$coefficients[2]) #add the current p-value from the randomized single linear regressions to the list of stored slopes
        }
        
        #extracting the real slope from our points
        lm_real <- lm(fixed_field_data_processed_size_variables[,k] ~ fixed_field_data_processed_soils.condensed[,j]) #creating the linear regression of the size metric vs. soil metrics for the respective population
        lm_real_sum <- summary(lm_real) #extract the summary 
        lm_real_slope <- lm_real_sum$coefficients[2] #storing the slope from the lienar regression summary
        
        #creating the histogram plot of the randomly distributed p-values and our real slope
        Size_Variable = colnames(fixed_field_data_processed_soils.condensed)[j] #extracting the size variable name
        plot_out <- ggplot()+
          geom_histogram(aes(x=slopes),  fill = "dodgerblue1", color = "black", bins = 50 )+
          geom_vline(xintercept=lm_real_slope, col = "red") + #line of our real slope
          xlab(paste(Population, "Slopes of Shuffled", Size_Variable, "vs. our", Size_Variable)) +
          ylab("Frequency") +
          theme_classic()
        
        # store the histogram plot list with a descriptive name 
        plot_name <- paste(Populations[i], #population name
                           colnames(fixed_field_data_processed_soils.condensed)[j], #soil metric name
                           colnames(fixed_field_data_processed_size_variables)[k], #size metric name
                           sep = "_")
        plot_list[[plot_name]] <- plot_out
        
        #calculating pseudo p-value based on where real slope overlaps with histogram of randomized slopes
        
        # if using greater than hypothesis (which we did not use in this case)
        
          # total = 0  #set empty value
          # for (p in 1:length(slopes)){ #loop that adds 1 to the value total if the simulated slope value is less than our average value for our trees
          #   if (slopes[p] > lm_real_slope){
          #     total = total + 1
          #   }
          # } #add number of values of in the random set of slope values that are less than our mean slope
          # p_value <- (total / length(slopes)) #the proportion of random ANNs that are less than our slope
        
        # using the significantly different alternative hypothesis 
        
        p_value_greater_than <- sum(slopes >= lm_real_slope)/length(slopes)   # proportion of simulated slopes higher than our real slope
        p_value_less_than <- sum(slopes <= lm_real_slope)/length(slopes)   # proportion of simulated slopes lower than our real slope
        p_value <- min(1, 2 * min(p_value_greater_than, p_value_less_than)) # find twice the smaller tail (the "more extreme" one) and if the value is over 1, use 1
        
        #storing the results
        
        slopes_array[j, k, i] = lm_real_slope #assigning the each index in the array with the new slope value
        p_values_array[j, k, i] = p_value #assigning the each index in the array with the new p value

       # print(paste("Updating:", j, k, i)) #soil metric, size metric, population
       # print(paste("Real slope:", lm_real_slope))
       # print(paste("P-value:", p_value))
        
      }
    }
    
  }
  
  #returning the slope array, p-value array, and histogram plot list
  return(list(slopes_array = slopes_array, p_values_array = p_values_array, plot_list = plot_list))
  
}

#### Running the Function and Storing the Results ####

#running the simulation function
slope_simultations <- slopes_simulations()

#Example of extracting one of the histograms comparing the slopes for our original soil vs. size metrics to the shuffled ones
slope_simultations$plot_list$LM_clay.content.0.5_Canopy_short
plot <- slope_simultations$plot_list$LM_clay.content.0.5_Canopy_short #storing the example plot
plot +  #example of changing the x-axis label
  xlab("LM Slopes of Linear Regressions with Shuffled Short Canopy Axis vs Known Clay Content (0-5 cm)")

#if you want to see all of the plots at once run: 
    #slope_simultations$plot_list

#turning the simulation results into a data frame

slopes_df <- as.data.frame.table(slope_simultations$slopes_array, 
                                 responseName = "slope") #creating a dataframe from the slope array
pvals_df <- as.data.frame.table(slope_simultations$p_values_array, 
                                responseName = "p_value") #creating a dataframe from the p value array

# merging the slopes and p-values dataframes 
size.pop.slopes.df <- merge(slopes_df, pvals_df, by = c("Var3", "Var2", "Var1")) #merging the two dataframes into one
names(size.pop.slopes.df) <- c("Population", "Size.Variable", "Soil.Metric", "Slope", "P.value") #re-naming the columns to be more appropriate

# Bonferroni correction of the p-values because of multiple testing
size.pop.slopes.df <- size.pop.slopes.df %>%
  mutate(p_bonf_corrected = as.numeric(p.adjust(size.pop.slopes.df$P.value, method = "bonferroni")))

#creating a column for whether the p values are significant or not
size.pop.slopes.df <- size.pop.slopes.df %>%
  mutate(Significance.reg = case_when(P.value < 0.05 ~ "Y",
                                  P.value >= 0.05 ~ "N")) %>% #using uncorrected p-values
  mutate(Significance.bonf = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                       p_bonf_corrected >= 0.05 ~ "N")) #using Bonferonni corrected p-values

#ensuring the categorical variables have the desired order of levels
size.pop.slopes.df$Population <- factor(size.pop.slopes.df$Population, levels = c("LM", "LC", "SD"))
size.pop.slopes.df$Size.Variable <- factor(size.pop.slopes.df$Size.Variable, levels = c("SCA", "LCA", "CA", "CS", "DBH"))
size.pop.slopes.df$Soil.Metric <- factor(size.pop.slopes.df$Soil.Metric, levels = c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                                                                                    "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                                                                                    "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                                                                                    "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                                                                                    "Volume of water content -1500 kpa 100-200", 
                                                                                    "Nitrogen 0-5", "Nitrogen 100-200", 
                                                                                    "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                                                                                    "Sand Available Water 0-5", "Sand Available Water 100-200",
                                                                                    "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200"))
#rearranging the rows based on the specified levels for the categorical variables
size.pop.slopes.df <- size.pop.slopes.df %>%
  arrange(Population, Size.Variable, Soil.Metric)
  
#### Summarizing the results ####

#summarizing the results data frame 
summary(size.pop.slopes.df)

## Tables

#table of number of yes or no by population 
table(size.pop.slopes.df$Population, size.pop.slopes.df$Significance.bonf)

# table of number of yes or nos for each size variables 

#LM
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"])

# table of number of yes or nos for each soil variables 

#LM
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"])

#LC
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"])

#SD
table(size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"])


# table of number of yes or nos for each size variables and each soil variable 

#LM 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LM"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LM"])

#LC 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="LC"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="LC"])

#SD 
table(size.pop.slopes.df$Size.Variable[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Significance.bonf[size.pop.slopes.df$Population=="SD"], size.pop.slopes.df$Soil.Metric[size.pop.slopes.df$Population=="SD"])

## plots to visualize the results

#plot of population, yes/no, and the interaction between size/shape and soil values
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Size.Variable,Soil.Metric)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, andsize/shape 
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Size.Variable)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()

#plot of population, yes/no, and soil values
ggplot(size.pop.slopes.df, aes(x = Significance.bonf, fill = interaction(Soil.Metric)))+
  geom_bar(position = 'stack')+
  facet_wrap(~Population)+
  theme_minimal()


## Heat Maps

#setting population to a factor and ensuring the appropriate levels
size.pop.slopes.df$Population <- as.factor(size.pop.slopes.df$Population)
levels(size.pop.slopes.df$Population)

#changing the number of digits included to be 3
options(digits=3) 
size.pop.slopes.df$P.value
#across all populations heat map of Bonferonni corrected P-values 
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Values",
       title = "All Sites Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) 
  #theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=0.5)) +
  #coord_flip()

#across LM, heat map of Bonferonni corrected P-values 
size.pop.slopes.df.LM <- size.pop.slopes.df[size.pop.slopes.df$Population == "LM", ] #isolating the LM p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(P.value < 0.05, P.value, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P.value < 0.05, P.value, NA))) 

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LM) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value", 
       title = "Las Matancitas Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, round(Slope, 3), NA)), color = "black")



#across LC, heat map of Bonferonni corrected P-values 
size.pop.slopes.df.LC <- size.pop.slopes.df[size.pop.slopes.df$Population == "LC", ] #isolating the LC p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P-values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)))

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.LC) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", fill = "P Value",  
       title = "La Cobriza Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1)  + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, round(Slope, 3), NA)))


#across SD
size.pop.slopes.df.SD <- size.pop.slopes.df[size.pop.slopes.df$Population == "SD", ] #isolating the SD p-values/slopes

#with labeled p-values
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant P Values Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, p_bonf_corrected, NA)))

#with labeled slopes
ggplot(aes(x = Size.Variable, y = Soil.Metric, fill = ifelse(P.value < 0.05, P.value, NA)), data = size.pop.slopes.df.SD) +
  geom_tile() + 
  labs(x = "Size Characteristic", y = "Soil Characteristic", 
       fill = "P Value",  
       title = "San Dionisio Significant Associations between Soil and Size/Shape",
       subtitle = "Significant Slopes Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P.value < 0.05, round(Slope, 3), NA))) 

