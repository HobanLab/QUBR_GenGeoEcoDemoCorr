# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Quercus brandegeei Population Locations Are Influenced by Their Soil Metrics%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to determine if the soil of the Quercus brandegeei locations are 
# significantly different than would be expected if the populations were randomly distributed around Baja California Sur.
# A significant difference could indicate which soil metrics may explain why the trees are found only in these
# distinct riverbed pockets. 

# To test this, we used soil metric spatial null model analysis whereby we compared the mean soil metric of the 
# 20 known populations to a distribution of the mean soil metric values for randomly permuted populations across BCS,
# to support the populations have significantly different soil metrics compared to other parts of BCS without assuming
# a normal distribution. 

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
        # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
        # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
#boxs around the  rivers, stacking the rasters for each population, and processing them into one dataframe for all and each population
        # Creating the four new soil metrics: Sandy available Water (0-5 and 100-200 cm) and Clay/Loam Available Water (0-5 and 100-200 cm)
        # Loading and processing a dataframe of the 20 known QUBR population locations and creating cropped soil metric rasters for the 7km radius around the populations in BCS,
#and extracting the known population soil metric values
# 2) Making a function for calculating the real and randomized mean soil metrics, comparing the real
# and randomized mean metrics for each soil metric, and calculating the p-values. 
# 3) Running the function and storing the output, graphing the histograms/slopes, turning the results into a dataframe
# 4) Using a heat map to summarize the results (the mean soil metrics and their significances)

# NOTE: Uncomment and run line 48, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) #ggplot extension
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

# loading in the processed tree data 
# NOTE: Uncomment and run line 48, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
#source("./analyses/Data_Processing_Script.R")

## Finalizing the tree soil metric dataframe

#combining the LM, LC, and SD tree dataframes with the soil metrics and randomly chosen points within each grid cell

fixed_field_data_processed_soils <- rbind(LM_fixed_field_data_processed_soils, LC_fixed_field_data_processed_soils) #combining the LM and LC soil and randomly chosen tree data
fixed_field_data_processed_soils <- rbind(fixed_field_data_processed_soils, SD_fixed_field_data_processed_soils) #combining the SD tree point data to the LM and LC soil and randomly chosen tree point data

#### Creating the Function Comparing Average Soil Values from Inside Populations to Outside Populations ####

# This function iterates over each soil metric and compares the average soil metric across the known populations
# to a distribution of randomly generated distribution populations means to see if the known populations have significantly
# different soil metrics compared at random. 

# The components of the function are
  # 1) generating an empty list for the soil means and the p-values from comparing the real and randomized population soil means.
  # 2) a loop to iterate over each soil metric in which 
        # generate a list of random soil means from randomly distributed populations (n=20, 1000 permutations),
        # plots of the cropped Baja California Sur polygon and the randomly selected points,
        # store the soil mean for the known population,
        # plot and stored the figure showing the real soil mean and a histogram of the randomly generated soil means,
        # generate the p-values (using significantly different),
        # and store the real population soil eans and p-values in their respective lists

#creating the list of soil metrics to iterate over
Soil.metrics <- c("Clay 0-5", "Clay 100-200", "Silt 0-5", "Silt 100-200", "Sand 0-5", "Sand 100-200",
                  "Ph 0-5", "Ph 100-200",  "Volume of water content -10 kpa 0-5",
                  "Volume of water content -10 kpa 100-200", "Volume of water content -33 kpa 0-5",
                  "Volume of water content -33 kpa 100-200", "Volume of water content -1500 kpa 0-5", 
                  "Volume of water content -1500 kpa 100-200", 
                  "Nitrogen 0-5", "Nitrogen 100-200", 
                  "Soil Organic Carbon 0-5", "Soil Organic Carbon 100-200",
                  "Sand Available Water 0-5", "Sand Available Water 100-200",
                  "Clay/Loam Available Water 0-5", "Clay/Loam Available Water 100-200")

# The function the soil mean list and p-value list, the randomly generated point plots, and stored histogram plot list.

random_pop_soils <- function(){
  

  #creating empty list to collect know population soil means and p-values
  known_soil_means <- c()  #for the p values of the known population means
  random_soil_p_values <- c()  #for the p values of the randomly generated population means compared to our known soil mean
  
  #to store generated histogram plots
  plot_list <- list()   
  
  #loop iterating over each soil metric
  for (i in 1:length(Soil.metrics)){
    
    #assigning the population based on the current soil metric in the list
    if (Soil.metrics[i] == "Clay 0-5"){ 
      soil_stack = soil_stack_clay
      soil_metric = all_known_pop_soils$clay.content.0.5
    } else if (Soil.metrics[i] == "Clay 100-200"){
      soil_stack = soil_stack_clay
      soil_metric = all_known_pop_soils$clay.content.100.200
    } else if (Soil.metrics[i] == "Silt 0-5"){
      soil_stack = soil_stack_silt
      soil_metric = all_known_pop_soils$silt.0.5
    } else if (Soil.metrics[i] == "Silt 100-200"){
      soil_stack = soil_stack_silt
      soil_metric = all_known_pop_soils$silt.100.200
    } else if (Soil.metrics[i] == "Sand 0-5"){
      soil_stack = soil_stack_sand
      soil_metric = all_known_pop_soils$sand.0.5
    } else if (Soil.metrics[i] == "Sand 100-200"){
      soil_stack = soil_stack_sand
      soil_metric = all_known_pop_soils$sand.100.200
    } else if (Soil.metrics[i] == "Ph 0-5"){
      soil_stack = soil_stack_ph
      soil_metric = all_known_pop_soils$ph_0.5
    } else if (Soil.metrics[i] == "Ph 100-200"){
      soil_stack = soil_stack_ph
      soil_metric = all_known_pop_soils$ph_100.200
    } else if (Soil.metrics[i] == "Volume of water content -10 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_10kpa
      soil_metric = all_known_pop_soils$vol_water_.10_0.5
    } else if (Soil.metrics[i] == "Volume of water content -10 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_10kpa
      soil_metric = all_known_pop_soils$vol_water_.10_100.200
    } else if (Soil.metrics[i] == "Volume of water content -33 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_33kpa
      soil_metric = all_known_pop_soils$vol_water_0.5
    } else if (Soil.metrics[i] == "Volume of water content -33 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_33kpa
      soil_metric = all_known_pop_soils$vol_water_100.200
    } else if (Soil.metrics[i] == "Volume of water content -1500 kpa 0-5"){
      soil_stack = soil_stack_vol_wat_1500kpa
      soil_metric = all_known_pop_soils$vol_water_.1500kPa_0.5
    } else if (Soil.metrics[i] == "Volume of water content -1500 kpa 100-200"){
      soil_stack = soil_stack_vol_wat_1500kpa
      soil_metric = all_known_pop_soils$vol_water_.1500_100.200
    } else if (Soil.metrics[i] == "Nitrogen 0-5"){
      soil_stack = soil_stack_nitrogen
      soil_metric = all_known_pop_soils$nitrogen.0.5
    } else if (Soil.metrics[i] == "Nitrogen 100-200"){
      soil_stack = soil_stack_nitrogen
      soil_metric = all_known_pop_soils$nitrogen.100.200
    } else if (Soil.metrics[i] == "Soil Organic Carbon 0-5"){
      soil_stack = soil_stack_soc
      soil_metric = all_known_pop_soils$SOC.0.5
    } else if (Soil.metrics[i] == "Soil Organic Carbon 100-200"){
      soil_stack = soil_stack_soc
      soil_metric = all_known_pop_soils$SOC.100.200
    } else if (Soil.metrics[i] == "Sand Available Water 0-5"){
      soil_stack = soil_stack_sandy_water
      soil_metric = all_known_pop_soils$sandy_avail_water_0.5
    } else if (Soil.metrics[i] == "Sand Available Water 100-200"){
      soil_stack = soil_stack_sandy_water
      soil_metric = all_known_pop_soils$sandy_avail_water_100.200
    } else if (Soil.metrics[i] == "Clay/Loam Available Water 0-5"){
      soil_stack = soil_stack_clay_loam_water
      soil_metric = all_known_pop_soils$clay_loam_avail_water_0.5
    } else if (Soil.metrics[i] == "Clay/Loam Available Water 100-200"){
      soil_stack = soil_stack_clay_loam_water
      soil_metric = all_known_pop_soils$clay_loam_avail_water_100.200
    } 
    
    #creating empty list to collect means
    random_soil_means <- c()  #for the means of the randomly generated population means
     
    
    set.seed(20) #setting a seed
    
    #looping for 1000 permutations
    for (y in 1:1000){ 
      
      #generating and extracting the randomly distributed population soil metric mean
      random_20 <- st_sample(BCS_polygon_box_sf_cropped, 20) #select random 20 points within the cropped BCS polygon
      random_20 <- random_20 %>%
        st_as_sf() #making sure the random points are stored as simple features
      random_20_pop_soil <- raster::extract(soil_stack, random_20) #extracting the soil metrics for the random points
      
      #storing the mean of the soil metric of the randomly generated populations depending on if it is the 0-5 or 100-200 cm version of the rasters
      if (i %% 2 == 1){ #if the iteration we are on is odd, then we use the 0-5 cm variable
        random_mean <- mean(random_20_pop_soil[,1]) #storing the mean soil metric, using the 0-5 cm raster
      } else {  #if the iteration we are on is odd, then we use the 100-200 cm variable
        random_mean <- mean(random_20_pop_soil[,2]) #storing the mean soil metric, using the 100-200 cm raster
      }
      
      random_soil_means <- c(random_soil_means, random_mean) #adding the 0-5 or 100-200 cm means to the list of means
      
    }
    
    
    #plotting the randomly selected points on the cropped, buffered, and full Baja California Sur polygons
    random_points_BCS <- ggplot()+
      geom_sf(data=BCS_polygon_UTM)+
      geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
      geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
      geom_sf(data=random_20, color ="blue")
    
    #plotting the randomly generated points just on the cropped polygon
    random_points_BCS_crop <- ggplot()+
      geom_sf(data=BCS_polygon_box_sf_cropped, color = "red")+
      geom_sf(data=all_pop_locations.df_sf_trans_coordinates)+
      geom_sf(data=random_20, color ="blue")
  
    #storing the real soil metric mean for our known populations
    all_known_mean <- mean(soil_metric)
    
    #adding the current known population mean soil metric to the list
    known_soil_means <- c(known_soil_means, all_known_mean)
    
    #plotting the histogram of the randomly distributed means and our real mean
    plot_out_histogram <- ggplot(data.frame(random_soil_means = random_soil_means))+ #data.frame(random_soil_means = random_soil_means)
      geom_histogram(aes(x=random_soil_means),  fill = "dodgerblue1", color = "black", bins = 50 )+
      geom_vline(xintercept=all_known_mean, col = "red")+ #line of our real slope
      xlab(paste0("Mean ", Soil.metrics[i], " of Random Populations vs. Known Populations (n=20)"))+
      theme_classic()
    
    # store the histogram in list with a descriptive name
    plot_name_histogram <- paste(Soil.metrics[i], "Histogram",
                       sep = "_")
    plot_list[[plot_name_histogram]] <- plot_out_histogram
    
    random_soil_means <- na.omit(random_soil_means) #removing NAs
    
    # if using greater than hypothesis
    
    #calculating pseudo p-value for 
    # total = 0  #set empty value
    # for (k in 1:length(random_soil_means)){ #loop that adds 1 to the value total if the simulated ANN value is less than our average value for our trees
    #   if (random_soil_means[k] < all_known_mean){
    #     total = total + 1
    #   }
    # } #add number of values of in the random set of means values that are less than our mean ANN
    # random_p.value <- 1 - (total / length(random_soil_means)) #the proportion of random ANNs that are greater than our ANN
    
    # using the significantly different alternative hypothesis 
    
    p_value_greater_than <- sum(random_soil_means >= all_known_mean)/length(random_soil_means)   # proportion of simulated slopes higher than our real slope
    p_value_less_than <- sum(random_soil_means <= all_known_mean)/length(random_soil_means)   # proportion of simulated slopes lower than our real slope
    random_p.value <- min(1, 2 * min(p_value_greater_than, p_value_less_than)) # take the smaller tail (the "more extreme" one), then double it
    
    #adding the p value to total list of p-values for all soil metrics
    random_soil_p_values <- c(random_soil_p_values, random_p.value)
    
    #print(paste("Updating:", i)) # printing the iteration number we are currently on
    
  }
  
  return(list(known_soil_means = known_soil_means, #soil mean list
         random_soil_p_values = random_soil_p_values, #p-value list
         random_points_BCS = random_points_BCS, #randomly generated points plot
         random_points_BCS_crop = random_points_BCS_crop, #randomly generated points plot with the crop
         plot_list = plot_list)) #stored histogram plot list

}

#### Running and Storing the Function and its Results ####
random_pop_soils_function <- random_pop_soils()

#Example of extracting one of the histograms comparing the slopes for our original soil vs. size metrics to the shuffled ones
random_pop_soils_function$plot_list$`Clay 0-5_Histogram`
plot <- random_pop_soils_function$plot_list$`Clay 0-5_Histogram`

#if you want to see all of the plots at once run: 
#random_pop_soils_function$plot_list

# Bonferroni correction for the p-values to protect against issues due to multiple testing
p_bonf_corrected <- p.adjust(random_pop_soils_function$random_soil_p_values, method = "bonferroni")
  
#making a dataframe from the function output
random_pop.df <- data.frame("Soil.metrics" = Soil.metrics, 
                            "P_values" = random_pop_soils_function$random_soil_p_values, 
                            "P_values_bonf_corrected" = p_bonf_corrected,
                            "Significance" = c(rep(NA, 22)))

#creating the significance column for the p-values (p<0.05 is significant)
random_pop.df <- random_pop.df %>%
  mutate(Significance = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                  p_bonf_corrected >= 0.05 ~ "N"))


#### Generating a Heat Map ####


#heat map of Bonferonni corrected P-values with labeled p-values
ggplot(aes(x = fct_reorder(Soil.metrics, p_bonf_corrected), y = Significance, fill = p_bonf_corrected), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Association Between Soil Metrics and Population Locations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(p_bonf_corrected < 0.05, round(p_bonf_corrected, 4), NA)), col = "white") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))

#heat map of P-values without correction with labeled p-values
ggplot(aes(x = fct_reorder(Soil.metrics, P_values), y = Significance, fill = P_values), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Association Between Soil Metrics and Population Locations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + 
  geom_text(aes(label = ifelse(P_values < 0.05, round(P_values, 4), NA)), col = "white") +
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))

#Histograms of the soil metric comparisons that were significiant
random_pop_soils_function$plot_list$`Sand 0-5_Histogram`
random_pop_soils_function$plot_list$`Volume of water content -33 kpa 0-5_Histogram`
random_pop_soils_function$plot_list$`Volume of water content -33 kpa 100-200_Histogram`

