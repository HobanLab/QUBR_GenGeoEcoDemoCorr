# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if the Quercus brandegeei populations have different mean soil metrics%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is determine whether there are significantly different mean
# soil metrics (clay content, pH, etc.) between the three Quercus brandegeei populations.
# A lack of differences could indicate which characteristics may explain why the trees are found at
# these distinct sites versus other locations. 
# Observed differences could may explain other differences in tree growth/survival between the sites, 
# as well as relate to environmental/topographic characteristics. 

# To test this, we used difference in means tests based on which conditions were met for each soil metric.
    # If the residuals were not normal, we used a Kruskal-Wallis test and Post-Hoc Wilcoxon Rank Sum Test
    # If the residuals were normal but the variance was NOT equal, we used a Welch's ANOVA test and Post-Hoc Tamhane's T2 Test
    # If the residuals were normal and the variance was equal, we used a Traditional ANOVA test and Post-Hoc Pairwise T-test

# The script is broken into sections of 
# 1) Loading and processing the packages and processed data for the trees, topography, and soil metrics in the Las Matancitas,
#San Dionisio, and La Cobriza populations. The processed data used in this script includes:
        # Processing the tree spatial/size data and river outline shapefiles to be in UTM 12 N Equal Area Projection, fixing errors in elevation,
#generating river and point buffers and bounding boxes,
        # Processing the soil raster data: loading the data in projecting the data, cropping them to the bounding 
#boxs around the  rivers, stacking the rasters for each population, and processing them into one dataframe for all and each population
        # Creating the four new soil metrics: Sandy available Water (0-5 and 100-200 cm) and Clay/Loam Available Water (0-5 and 100-200 cm)
# 2) Choosing random trees per grid cell for each population to avoid independence issues during tests. 
# 3) Making the function for checking conditions and running the appropriate difference in means test
# 4) Running the function and storing the outputs for each soil metric
# 5) Presenting the results as a heat map

# NOTE: Uncomment and run line 48, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

#### Loading libraries and Relevant Data ####

library(tidyverse) # for graphing and data organization
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc) # ggplot extension
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(raster) #to plot rasters
library(rstatix) #to run the Games-Howell Test
library(ggnewscale) #to be able to assign different colors to different layered rasters

# loading in the processed tree data 
# NOTE: Uncomment and run the line below, line 48, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
# source("./analyses/Data_Processing_Script.R")

#### Choosing a Random Tree per Grid Cell ####

# For each population, randomly selecting one tree from each grid cell to avoid issues because
# the rasters are only 250 m resolution and certain populations might have more points in specific
# cells, skewing the mean results away from the potentially real mean

#LM

#creating a grid over the soil cells
LM_tree_grid_cropped <- st_make_grid(soil_stack_LM_soil_text, cellsize = c(230, 265)) #230, 265 were used for cellsize because that is the dimensions of soil stack

#plotting the grid over an example soil raster to ensure they were made properly
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_LM_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = LM_tree_grid_cropped, fill = NA)

#selecting a point from each grid cell with trees within them
LM_list_grids_and_points <- st_contains(LM_tree_grid_cropped, LM_fixed_field_data_processed_sf, sparse =T) #making sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
LM_list_grids_and_trees <- lapply(LM_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})


#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
LM_list_grids_and_point_trees_df <- as.data.frame(unlist(LM_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(LM_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name

#filters out grid cells that do not have trees within them
LM_list_grids_and_trees_fixed <- LM_list_grids_and_point_trees_df %>%
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  mutate(data_row = LM_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
LM_fixed_field_data_processed_trees_soils <- LM_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% LM_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with the row numbers that match between the overall tree points dataframe and the focal tree points dataframe

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LM_tree_grid_cropped)+
  geom_sf(data= LM_fixed_field_data_processed_sf)+
  geom_sf(data = LM_fixed_field_data_processed_trees_soils, color = "red")


#LC

#creating a grid over the soil cells
LC_tree_grid_cropped <- st_make_grid(soil_stack_LC_soil_text, cellsize = c(230, 265)) #230, 265 were used for cellsize because that is the dimensions of soil stack

#plotting the grid over an example soil raster to ensure they were made properly
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_LC_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = LC_tree_grid_cropped, fill = NA)

#selecting a point from each grid cell with trees within them
LC_list_grids_and_points <- st_contains(LC_tree_grid_cropped, LC_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
LC_list_grids_and_trees <- lapply(LC_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})

#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
LC_list_grids_and_point_trees_df <- as.data.frame(unlist(LC_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(LC_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name

#filters out grid cells that do not have trees within them
LC_list_grids_and_trees_fixed <- LC_list_grids_and_point_trees_df %>%
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  mutate(data_row = LC_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
LC_fixed_field_data_processed_trees_soils <- LC_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% LC_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with the row numbers that match between the overall tree points dataframe and the focal tree points dataframe

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = LC_tree_grid_cropped)+
  geom_sf(data= LC_fixed_field_data_processed_sf)+
  geom_sf(data = LC_fixed_field_data_processed_trees_soils, color = "red")

#SD

#creating a grid over the soil cells
SD_tree_grid_cropped <- st_make_grid(soil_stack_SD_soil_text, cellsize = c(230, 265)) #230, 265 were used for cellsize because that is the dimensions of soil stack

#plotting the grid over an example soil raster to ensure they were made properly
ggplot()+
  geom_raster(data= as.data.frame(soil_stack_SD_soil_text, xy = T), aes(x=x, y=y, fill = clay.content.0.5))+
  geom_sf(data = SD_tree_grid_cropped, fill = NA)

#selecting a point from each grid cell with trees within them
SD_list_grids_and_points <- st_contains(SD_tree_grid_cropped, SD_fixed_field_data_processed_sf, sparse =T) #make sure row number in the data frame of grid cells corresponds to the order of what is in the points dataframe within st_contains
set.seed(24) #setting the seed
SD_list_grids_and_trees <- lapply(SD_list_grids_and_points, function(cell){ #iterates over the list of each grid cell with what row of points is within that grid cell made by st_contains
  if(length(cell) > 1){ #for each grid cell, if there is more than one tree in each cell
    tree_pt <- sample(cell, size = 1, replace = F) #randomly select a row from the row of trees within that polygon
  }
  else if(length(cell) == 1) { #for each grid cell, if there is exactly one tree in each cell
    tree_pt <- cell #set the selected tree point to be the tree that is within the cell
  } else { # if there are no trees
    tree_pt <- NA # set the focal tree point to be NA
  }
  return(tree_pt)
})


#creating a dataframe of all of the trees with their row number in the overall tree point dataframe and in which grid cell they are in
SD_list_grids_and_point_trees_df <- as.data.frame(unlist(SD_list_grids_and_trees)) #turns the list of grid cells and what focal trees were within them into a dataframe
colnames(SD_list_grids_and_point_trees_df) <- c("tree_row_num") #changes the column name
#filters out grid cells that do not have trees within them
SD_list_grids_and_trees_fixed <- SD_list_grids_and_point_trees_df %>%
  mutate(cell_num = row_number()) %>% #assigns the cell number to each row/tree
  mutate(data_row = SD_fixed_field_data_processed$X[tree_row_num]) %>% #adding a column that writes the real row number the focal tree is in the overall data
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the focal points
SD_fixed_field_data_processed_trees_soils <- SD_fixed_field_data_processed_soils %>%
  filter(X_sequential %in% SD_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with row numbers that match between the overall tree points dataframe and the focal tree points dataframe

#plotting the points, grid, and randomly selected points from each grid
ggplot()+
  geom_sf(data = SD_tree_grid_cropped)+
  geom_sf(data= SD_fixed_field_data_processed_sf)+
  geom_sf(data = SD_fixed_field_data_processed_trees_soils, color = "red")

## Finalizing the tree soil metric dataframe

#combining the LM, LC, and SD tree dataframes with the soil metrics and randomly chosen points within each grid cell

fixed_field_data_processed_trees_soils <- rbind(LM_fixed_field_data_processed_trees_soils, LC_fixed_field_data_processed_trees_soils) #combining the LM and LC soil and randomly chosen tree data
fixed_field_data_processed_trees_soils <- rbind(fixed_field_data_processed_trees_soils, SD_fixed_field_data_processed_trees_soils) #combining the SD tree point data to the LM and LC soil and randomly chosen tree point data

#creating a column/variable with locality as a factor to be able to use it in the Tamhane's T2 Test later

fixed_field_data_processed_trees_soils$Locality_Factor <- as.factor(fixed_field_data_processed_trees_soils$Locality)

#### Making Function for Differences in Means ####

# Function that checks the conditions and then runs the appropriate difference in means test
    # If the residuals were not normal, we used a Kruskal-Wallis test and Post-Hoc Wilcoxon Rank Sum Test.
        # We also performed this test for every soil metric because it is a non-parametric test allowing for comparisons across the soil metrics. 
    # If the residuals were normal but the variance was NOT equal, we used a Welch's ANOVA test and Post-Hoc Tamhane's T2 Test
    # If the residuals were normal and the variance was equal, we used a Traditional ANOVA test and Post-Hoc Pairwise T-test

# For checking the residuals were normal, we used a Shapiro-Wilks Test
# For checking equal variance, we check that test outcomes were both true, 1) rule of thumb (< 2 means equal variance) test and 
    # if the residuals were normal, 
          #we used the 2) Bartlett's test (p>0.05 means equal variance) to test if there was equal variance because it works better with normal residuals
    # if the residuals were NOT normal,
          #we used the 2) Fligner-Killeen test to test (p>0.05 means equal variance) if there was equal variance because it works better with non-normal residuals
# We also computed the Levene's Test, as an additional, less robust test for additional evidence, but we did not factor it into our difference in means test decision.

# The function returns:
  # a) an ANOVA model, 
  # b) the Shapiro-Wilks Test results, 
  # c) Fligner-Killeen Test Results,
  # d) Bartlett's Test results,
  # e) Levene's Test results,
  # f) Rule of Thumb Test results,
  # g) the chosen Difference in Means Test results,
  # h) the chosen Difference in Means Test Post-Hoc Test results,
  # i) a print of the difference in means tests chosen,
  # j) Kruskal-Wallis Test results,
  # k) and Post-Hoc Wilcoxon Rank Sum Test results.

mean_soil_function <- function(soil_group, data = fixed_field_data_processed_trees_soils, Populations = "Locality") {
  
  # Building the formula for the difference in means tests
  formula <- as.formula(paste(soil_group, "~", Populations)) #creating the formula of the soil group as the y variable and the locality as the x variable for expediting the following tests
  
  # Initial ANOVA test 
  anova_model <- aov(formula, data = data) 
  
  # checking to see if the residuals are normal
  shapiro_test <- shapiro.test(anova_model$residuals) #Shapiro-Wilks Test
  
  #Equal variance tests
  
  #Fligner-Killeen, more useful when data is not normal or there are outliers 
  fligner_test <- fligner.test(formula, data = fixed_field_data_processed_trees_soils)
  #bartlett's test for equal variances when data is normal
  bartlett_test <- bartlett.test(formula, data = fixed_field_data_processed_trees_soils)
  #levene's test, not super robust to strong differences to normality
  levenes_test <- car::leveneTest(formula, data = fixed_field_data_processed_trees_soils)
  #rule of thumb test
  thumb_test <- tapply(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality, sd)
  thumb_test_results <- max(thumb_test, na.rm = T) / min(thumb_test, na.rm = T) # if the max sd divided by the min sd is greater than two,the test did not pass
  
  # checking conditions to choose with difference in means test to use
  if (shapiro_test$p.value < 0.05) { #if the residuals are NOT normally distributed
    #kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
    test <- kruskal.test(formula, data = fixed_field_data_processed_trees_soils)
    #post-hoc Wilcoxon rank sum tests
    post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality,
                         p.adjust.method = "fdr") #p value adjusted using false discovery rate method
    #storing the tests used
    test_type = "kruskal-Wallis + Wilcoxon Rank Sum Test"
    #printing out which test was used
    print(paste("kruskal-Wallis Test with a Post-Hoc Wilcoxon Rank Sum Test"))
    
  } else if (shapiro_test$p.value > 0.05) { #if the residuals are normally distributed
    if (bartlett_test$p.value < 0.05 & thumb_test_results > 2) { #if the equal variance of residuals condition is NOT met
      
      #Welch's ANOVA, does not assume equal variances, but does meet normality condition
      test <- oneway.test(formula, data = fixed_field_data_processed_trees_soils, var.equal = F)
      #post hoc Welch's ANOVA test: Tamhane's T2 Test
      post_hoc <- tamhaneT2Test(fixed_field_data_processed_trees_soils[[soil_group]]~fixed_field_data_processed_trees_soils$Locality_Factor, data = fixed_field_data_processed_trees_soils)
      #storing the tests used
      test_type = "Welch's ANOVA + Tamhane's Test"
      #printing out which test was used
      print(paste("Welch's ANOVA with a Post-Hoc Tamhane's Test"))
      
    } else if (bartlett_test$p.value > 0.05 & thumb_test_results < 2) { #if the equal variance of residuals condition IS met
      #traditional ANOVA because equal variance and normality were met
      test <- anova(anova_model)
      #post-hoc Tukey's HSD
      post_hoc <- TukeyHSD(anova_model)
      #storing the tests used
      test_type = "ANOVA + Tukey's HSD"
      #printing out which test was used
      print(paste("ANOVA Test with a Post-Hoc Tukey's HSD"))
    }
  }
  
  #kruskal-Wallis test because it is non-parametric
  kruskal_test <- kruskal.test(formula, data = fixed_field_data_processed_trees_soils)
  #post-hoc Wilcoxon rank sum tests
  kruskal_post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils[[soil_group]], fixed_field_data_processed_trees_soils$Locality,
                                   p.adjust.method = "fdr") #p value adjusted using false discovery rate method
  
  
  return(list(
    anova_model = summary(anova_model),
    shapiro_test = shapiro_test,
    fligner_test = fligner_test,
    bartlett_test = bartlett_test,
    levenes_test = levenes_test,
    thumb_test_results = thumb_test_results,
    final_test = test,
    posthoc = post_hoc,
    test_type = test_type,
    kruskal_test = kruskal_test, 
    kruskal_post_hoc = kruskal_post_hoc
  ))
}

#### Running the Function Results ####
  
##clay 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_0.5 <- mean_soil_function("clay.content.0.5")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_clay_0.5$final_test

#post hoc test of ANOVA test: Tamhane's Test
mean_soil_function_clay_0.5$posthoc

#Kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_0.5_mean_p.value <- mean_soil_function_clay_0.5$kruskal_test$p.value

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_0_5 <- aov(clay.content.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content vs. Population")


##clay 100-200 

#running the Difference in Means Analysis function
mean_soil_function_clay_100.200 <- mean_soil_function("clay.content.100.200")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_clay_100.200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_100.200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_100.200_mean_p.value <- mean_soil_function_clay_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay.content.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_100_200 <- aov(clay.content.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

#silt 0-5

#running the Difference in Means Analysis function
mean_soil_function_silt_0.5 <- mean_soil_function("silt.0.5")

#CHOSEN TEST: ANOVA test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_silt_0.5$final_test

#kruskal-Wallis test post-hoc Tukey's HSD sum tests
mean_soil_function_silt_0.5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_silt_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
silt_0.5_mean_p.value <- mean_soil_function_silt_0.5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_silt_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_silt_0_5 <- aov(silt.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_silt_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

##silt 100-200

#running the Difference in Means Analysis function
mean_soil_function_silt_100.200 <- mean_soil_function("silt.100.200")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_silt_100.200$final_test

#post-hoc Pairwise Tukey's HSD
mean_soil_function_silt_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_silt_100.200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
silt_100.200_mean_p.value <- mean_soil_function_silt_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_silt_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, silt.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_silt_100_200 <- aov(silt.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_silt_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

##sand  0-5 

#running the Difference in Means Analysis function
mean_soil_function_sand_0.5 <- mean_soil_function("sand.0.5")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_sand_0.5$final_test

#post-hoc Pairwise Tukey's HSD Tests
mean_soil_function_sand_0.5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sand_0.5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sand_0.5_mean_p.value <- mean_soil_function_sand_0.5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sand_0.5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_sand_0_5 <- aov(sand.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sand_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

## sand 100-200

#running the Difference in Means Analysis function
mean_soil_function_sand_100.200 <- mean_soil_function("sand.100.200")

#CHOSEN TEST: Kruskal-Wallis
mean_soil_function_sand_100.200$final_test

#post hoc Kruskal-Wallis test: pairwise wilcoxon test, fdr adjustment
mean_soil_function_sand_100.200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sand_100.200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sand_100.200_mean_p.value <- mean_soil_function_sand_100.200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sand_100.200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sand.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_sand_100_200 <- aov(sand.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sand_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay Content at 100-200 cm vs. Population")

## ph 0-5

#running the Difference in Means Analysis function
mean_soil_function_ph_0_5 <- mean_soil_function("ph_0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_ph_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_ph_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
ph_0.5_mean_p.value <- mean_soil_function_ph_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_ph_0_5 <- aov(ph_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_ph_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for pH at 0-5 cm vs. Population")

##ph 100-200

#running the Difference in Means Analysis function
mean_soil_function_ph_100_200 <- mean_soil_function("ph_100.200")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_ph_100_200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_ph_100_200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
ph_100.200_mean_p.value <- mean_soil_function_ph_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_ph_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, ph_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_ph_100_200 <- aov(ph_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_ph_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for pH at 100-200 cm vs. Population")

##soil organic carbon 0-5

#running the Difference in Means Analysis function
mean_soil_function_SOC_0_5 <- mean_soil_function("SOC.0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_SOC_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_SOC_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
soc_0.5_mean_p.value <- mean_soil_function_SOC_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_soc_0_5 <- aov(SOC.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_soc_0_5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Organic Carbon at 0-5 cm vs. Population")

#soil organic carbon 100-200

#running the Difference in Means Analysis function
mean_soil_function_SOC_100_200 <- mean_soil_function("SOC.100.200")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_SOC_100_200$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_SOC_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_SOC_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
soc_100.200_mean_p.value <- mean_soil_function_SOC_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_SOC_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, SOC.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_soc_100_200 <- aov(SOC.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_soc_100_200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Soil Organic Carbon at 100-200 cm vs. Population")

#volume of water content at -10 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_10_0_5 <- mean_soil_function("vol_water_.10_0.5")

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_10_0_5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_10_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_10kpa_0.5_mean_p.value <- mean_soil_function_vol_water_10_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_10_0.5 <- aov(vol_water_.10_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_10_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -10 kpa at 0-5 cm vs. Population")

#volume of water content at -10 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_10_100_200 <- mean_soil_function("vol_water_.10_100.200")

#CHOSEN TEST: ANOVA test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_10_100_200$final_test

#ANOVA test post-hoc Tukey's HSD tests
mean_soil_function_vol_water_10_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_10_100_200$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_10kpa_100.200_mean_p.value <- mean_soil_function_vol_water_10_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_10_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.10_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_10_100.200 <- aov(vol_water_.10_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_10_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -10 kpa at 100-200 cm vs. Population")

#volume of water content at -33 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_33_0_5 <- mean_soil_function("vol_water_0.5")

#CHOSEN TEST: kruskal-Wallis, does not assume equal variances 
mean_soil_function_vol_water_33_0_5$final_test

#post hoc kruskal-Wallis  test: Wilcoxon Rank Sum Test
mean_soil_function_vol_water_33_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_33_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_33kpa_0.5_mean_p.value <- mean_soil_function_vol_water_33_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_33_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_33_0.5 <- aov(vol_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_33_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -33 kpa at 0-5 cm vs. Population")

#volume of water content at -33 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_33_100_200 <- mean_soil_function("vol_water_100.200")

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_33_100_200$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
mean_soil_function_vol_water_33_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_33_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_33kpa_100.200_mean_p.value <- mean_soil_function_vol_water_33_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_33_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_33_100.200 <- aov(vol_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_33_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -33 kpa at 100-200 cm vs. Population")

#volume of water content at -1500 kpa 0-5

#running the Difference in Means Analysis function
mean_soil_function_vol_water_1500_0_5 <- mean_soil_function("vol_water_.1500kPa_0.5")

#CHOSEN TEST: kruskal-Wallis Test, does not assume equal variances 
mean_soil_function_vol_water_1500_0_5$final_test

#post hoc kruskal-Wallis test: Wilcoxon Rank Sum Test
mean_soil_function_vol_water_1500_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_1500_0_5$kruskal_test

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_1500kpa_0.5_mean_p.value <- mean_soil_function_vol_water_1500_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500kPa_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_1500_0.5 <- aov(vol_water_.1500kPa_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_1500_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -1500 kpa at 0-5 cm vs. Population")

#volume of water content at -1500 kpa 100-200

#running the Difference in Means Analysis function
mean_soil_function_vol_water_1500_100_200 <- mean_soil_function("vol_water_.1500_100.200")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_1500_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_vol_water_1500_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
vol_wat_1500kpa_100.200_mean_p.value <- mean_soil_function_vol_water_1500_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_vol_water_1500_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, vol_water_.1500_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_vol_water_1500_100.200 <- aov(vol_water_.1500_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_vol_water_1500_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Volume of Water Content at -1500 kpa at 100-200 cm vs. Population")

#nitrogen 0-5

#running the Difference in Means Analysis function
mean_soil_function_nitrogen_0_5 <- mean_soil_function("nitrogen.0.5")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_nitrogen_0_5$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_nitrogen_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
nitrogen_0.5_mean_p.value <- mean_soil_function_nitrogen_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_nitrogen_0.5 <- aov(nitrogen.0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_nitrogen_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Nitrogen Content at 0-5 cm vs. Population")

# nitrogen 100-200

#running the Difference in Means Analysis function
mean_soil_function_nitrogen_100_200 <- mean_soil_function("nitrogen.100.200")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_nitrogen_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_nitrogen_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
nitrogen_100.200_mean_p.value <- mean_soil_function_nitrogen_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_nitrogen_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, nitrogen.100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_nitrogen_100.200 <- aov(nitrogen.100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_nitrogen_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Nitrogen Content at 100-200 cm vs. Population")

# sandy available water 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_sandy_avail_water_0_5 <- mean_soil_function("sandy_avail_water_0.5")

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_sandy_avail_water_0_5$final_test

#post-hoc Pairwise T-Tests
mean_soil_function_sandy_avail_water_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sandy_avail_water_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sandy_avail_water_0.5_mean_p.value <- mean_soil_function_sandy_avail_water_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sandy_avail_water_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sandy_avail_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_sandy_avail_water_0.5 <- aov(sandy_avail_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sandy_avail_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Sand Available Water at 0-5 cm vs. Population")

# sandy available water 100-200 cm

#running the Difference in Means Analysis function
mean_soil_function_sandy_avail_water_100_200 <- mean_soil_function("sandy_avail_water_100.200")

#CHOSEN TEST: Welch's ANOVA test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_sandy_avail_water_100_200$final_test

#Welch's ANOVA test post-hoc Tamhane's tests
mean_soil_function_sandy_avail_water_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_sandy_avail_water_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
sandy_avail_water_100.200_mean_p.value <- mean_soil_function_sandy_avail_water_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_sandy_avail_water_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, sandy_avail_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_sandy_avail_water_100.200 <- aov(sandy_avail_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_sandy_avail_water_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Sand Available Water at 100-200 cm vs. Population")

# clay loam available water 0-5 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_loam_avail_water_0_5 <- mean_soil_function("clay_loam_avail_water_0.5")

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_clay_loam_avail_water_0_5$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_0_5$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_loam_avail_water_0_5$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_loam_avail_water_0.5_mean_p.value <- mean_soil_function_clay_loam_avail_water_0_5$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_0_5$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay_loam_avail_water_0.5))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_loam_avail_water_0.5 <- aov(clay_loam_avail_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_loam_avail_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

# clay loam available water 100-200 cm

#running the Difference in Means Analysis function
mean_soil_function_clay_loam_avail_water_100_200 <- mean_soil_function("clay_loam_avail_water_100.200")

#CHOSEN TEST: Kruskal Wallis test, assumes equal variance and normal residuals
mean_soil_function_clay_loam_avail_water_100_200$final_test

#post-hoc Wilcoxon Rank Sum Tests
mean_soil_function_clay_loam_avail_water_100_200$posthoc

#kruskal-Wallis test because it is non-parametric and comparable across soil metrics
mean_soil_function_clay_loam_avail_water_100_200$kruskal_test 

#storing the Kruskal-Wallis Test result p-values for a heat map
clay_loam_avail_water_100.200_mean_p.value <- mean_soil_function_clay_loam_avail_water_100_200$kruskal_test$p.value

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
mean_soil_function_clay_loam_avail_water_100_200$kruskal_post_hoc

#boxplots to show the spread of data
ggplot()+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, clay_loam_avail_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_clay_loam_avail_water_100.200 <- aov(clay_loam_avail_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_clay_loam_avail_water_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 100-200 cm vs. Population")


#### Presenting Results ####

#Heat Map 

#Storing the p-values from the chosen Difference in Means test in a vector
p_value_mean <- c(clay_0.5_mean_p.value, clay_100.200_mean_p.value, silt_0.5_mean_p.value,
                    silt_100.200_mean_p.value,
                    sand_0.5_mean_p.value,
                    sand_100.200_mean_p.value,
                    ph_0.5_mean_p.value,
                    ph_100.200_mean_p.value,
                    soc_0.5_mean_p.value,
                    soc_100.200_mean_p.value,
                    vol_wat_10kpa_0.5_mean_p.value,
                    vol_wat_10kpa_100.200_mean_p.value,
                    vol_wat_33kpa_0.5_mean_p.value,
                    vol_wat_33kpa_100.200_mean_p.value,
                    vol_wat_1500kpa_0.5_mean_p.value,
                    vol_wat_1500kpa_100.200_mean_p.value,
                    nitrogen_0.5_mean_p.value,
                    nitrogen_100.200_mean_p.value,
                    sandy_avail_water_0.5_mean_p.value,
                    sandy_avail_water_100.200_mean_p.value,
                    clay_loam_avail_water_0.5_mean_p.value,
                    clay_loam_avail_water_100.200_mean_p.value
)



# Bonferroni correction of the p-values because of multiple testing
p_bonf_corrected <- p.adjust(p_value_mean, method = "bonferroni")

#creating empty dataframe for inputting the function into
random_pop.df <- data.frame("Shape.Size" = rep(c("Clay 0-5 cm", "Clay 100-200 cm", "Silt 0-5 cm", "Silt 100-200 cm", "Sand 0-5 cm", "Sand 100-200 cm", #column of the soil metric names
                                                 "pH 0-5 cm", "pH 100-200 cm", "Soil Organic Carbon 0-5 cm", "Soil Organic Carbon 100-200 cm", 
                                                 "Volume of water content -10 kPa 0-5 cm", "Volume of water content -10 kPa 100-200 cm",
                                                 "Volume of water content -33 kPa 0-5 cm", "Volume of water content -33 kPa 100-200 cm",
                                                 "Volume of water content -1500 kPa 0-5 cm", "Volume of water content -1500 kPa 100-200 cm", 
                                                 "Nitrogen 0-5 cm", "Nitrogen 100-200 cm", "Sand Available Water 0-5 cm", "Sand Available Water 100-200 cm",
                                                 "Clay/Loam Available Water 0-5 cm", "Clay/Loam Available Water 100-200 cm")),
                            "P_Value" = p_bonf_corrected, #Bonferonni-corrected p-values
                            "Significance" = c(rep(NA, 22))) #whether the p-values are significant (p<0.05) or not

#creating the significance column based on the significance of the p-values
random_pop.df <- random_pop.df %>%
  mutate(Significance = case_when(p_bonf_corrected < 0.05 ~ "Y",
                                  p_bonf_corrected >= 0.05 ~ "N"))

#Turning off scientific notation format
options(scipen=999)

#Creating a heat map in ggplot of the p-values for each soil metric (y-axis) and whether they are 
#significant or not (x-axis), with significant p-values labeled
ggplot(aes(x = fct_reorder(Shape.Size, P_Value), y = Significance, fill = P_Value), data = random_pop.df) +
  geom_tile() + 
  labs(y = "Significant P-Value", x  = "Soil Characteristic", 
       fill = "P-Value",  
       title = "Difference Between Mean Soil Metrics Between Populations",
       subtitle = "P-Values Below 0.5 Labeled") + 
  scale_fill_distiller(palette = "RdPu", direction = -1) + #setting color pallete
  geom_text(aes(label = ifelse(P_Value < 0.001, "< 0.001", NA)), col = "white") + #labeling significant p-values if less than 0.001 as "<0.001"
  geom_text(aes(label = ifelse(P_Value < 0.5 & P_Value > 0.001, round(P_Value, 8), NA)), col = "white") + #labeling cells with significant p-values less then 0.05 and greater than 0.001, rounding to 8 decimal places
  coord_flip() +
  theme_classic() +
  theme(axis.text = element_text(size = 13),
        axis.title = element_text(size=13),
        title = element_text(size = 13),
        legend.title = element_text(size = 13),
        plot.subtitle = element_text(size = 12))


#### Comparing if there is a significant difference between sand available water for and LM and clay/loam available water for SD ####

#creating available water column that allows for a mix of sand and clay/loam available water
fixed_field_data_processed_trees_soils <- fixed_field_data_processed_trees_soils %>%
  mutate(available_water_0.5 = case_when(Locality == "LM" ~ sandy_avail_water_0.5,
                                         Locality == "LC" ~ sandy_avail_water_0.5,
                                         Locality == "SD" ~ clay_loam_avail_water_0.5),
         available_water_100.200 = case_when(Locality == "LM" ~ sandy_avail_water_100.200,
                                             Locality == "LC" ~ sandy_avail_water_100.200,
                                             Locality == "SD" ~ clay_loam_avail_water_100.200))


#0.5 cm available water
ggplot()+
  labs(y = "Available Water", title = "Blue is 0.5 cm, SD is clay/Loam and LM/LC is sand")+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, available_water_0.5), fill = "light blue")+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "LM"), aes(Locality, sandy_avail_water_0.5), fill = "light blue")+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "LC"), aes(Locality, sandy_avail_water_0.5), fill = "light blue")+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "SD"), aes(Locality, clay_loam_avail_water_0.5), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(available_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(available_water_0.5 ~ Locality, data = fixed_field_data_processed_trees_soils)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils$available_water_0.5, fixed_field_data_processed_trees_soils$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc

#100-200 cm boxplots to show the spread of data
ggplot()+
  labs(y = "Available Water", title = "100.200 Available Water, SD is clay/Loam and LM/LC is sand")+
  geom_boxplot(data = fixed_field_data_processed_trees_soils, aes(Locality, available_water_100.200), fill = "light blue")+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "LM"), aes(Locality, sandy_avail_water_100.200))+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "LC"), aes(Locality, sandy_avail_water_100.200))+
  # geom_boxplot(data = subset(fixed_field_data_processed_trees_soils, Locality == "SD"), aes(Locality, clay_loam_avail_water_100.200))+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_100.200 <- aov(available_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
hist(anova_available_water_100.200$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(available_water_100.200 ~ Locality, data = fixed_field_data_processed_trees_soils)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(fixed_field_data_processed_trees_soils$available_water_100.200, fixed_field_data_processed_trees_soils$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc


#### Plotting the Soil Classifications of the Different Populations ####

#installing/loading in the package 
devtools::install_github("Saryace/ggsoiltexture")
library(ggsoiltexture)

#creating a dataframe with a copy of the fixed_field_data_processed_trees_soils dataframe
fixed_field_data_processed_trees_soils.texture.triangle <- fixed_field_data_processed_trees_soils

#converting them from g/kg to percentage (%)
fixed_field_data_processed_trees_soils.texture.triangle <- fixed_field_data_processed_trees_soils.texture.triangle %>%
  mutate(clay.content.0.5.perc = round(clay.content.0.5/10, digits = 0)) %>%
  mutate(silt.0.5.perc = round(silt.0.5/10, digits = 0)) %>%
  mutate(sand.0.5.perc = round(sand.0.5/10, digits = 0)) %>%
  mutate(clay.content.100.200.perc = round(clay.content.100.200/10, digits = 0)) %>%
  mutate(silt.100.200.perc = round(silt.100.200/10, digits = 0)) %>%
  mutate(sand.100.200.perc = round(sand.100.200/10, digits = 0))

## 0-5 cm

#fixing column names so they will be usable for the soil texture triangle 
fixed_field_data_processed_trees_soils.texture.triangle.0.5 <- fixed_field_data_processed_trees_soils.texture.triangle %>%
  rename(clay = clay.content.0.5.perc) %>%
  rename(silt = silt.0.5.perc) %>%
  rename(sand = sand.0.5.perc)

#make a list of the rows that have silt, clay, and sand percentages that do not add to 100%
remove <- which(fixed_field_data_processed_trees_soils.texture.triangle.0.5$clay+fixed_field_data_processed_trees_soils.texture.triangle.0.5$silt+fixed_field_data_processed_trees_soils.texture.triangle.0.5$sand != 100)


#replacing sand, silt, clay values to get a total of 100

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[1],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[1],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = 70.4)

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[2],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[2],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = 17.7) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[3],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[3],] %>%
  mutate(clay = 11) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[4],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[4],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = 70.7)

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[5],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[5],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = 18.2) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[6],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[6],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[7],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[7],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[8],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[8],] %>%
  mutate(clay = 16.2) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[9],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[9],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = 66.6)

fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[10],] <- fixed_field_data_processed_trees_soils.texture.triangle.0.5[remove[10],] %>%
  mutate(clay = round(clay.content.0.5/10, digits = 1)) %>%
  mutate(silt = round(silt.0.5/10, digits = 1)) %>%
  mutate(sand = round(sand.0.5/10, digits = 1))

#plotting the soil texture triangle
plot.0.5 <- 
  ggsoiltexture(fixed_field_data_processed_trees_soils.texture.triangle.0.5, class = "USDA") +
  geom_point(aes(color = Locality), size = 1) +
  scale_color_discrete(type = "viridis")+
  labs(title = "Soi Texture 0-5 cm")

plot.0.5

## 100-200 cm

fixed_field_data_processed_trees_soils.texture.triangle.100.200 <- st_drop_geometry(fixed_field_data_processed_trees_soils.texture.triangle.100.200)

#fixing column names so they will be usable for the soil texture triangle 
fixed_field_data_processed_trees_soils.texture.triangle.100.200 <- fixed_field_data_processed_trees_soils.texture.triangle %>%
  rename(clay = clay.content.100.200.perc) %>%
  rename(silt = silt.100.200.perc) %>%
  rename(sand = sand.100.200.perc)

remove <- which(fixed_field_data_processed_trees_soils.texture.triangle.100.200$clay+fixed_field_data_processed_trees_soils.texture.triangle.100.200$silt+fixed_field_data_processed_trees_soils.texture.triangle.100.200$sand != 100)

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[1],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[1],] %>%
  mutate(clay = round(clay.content.100.200/10, digits = 1)) %>%
  mutate(silt = round(silt.100.200/10, digits = 1)) %>%
  mutate(sand = round(sand.100.200/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[2],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[2],] %>%
  mutate(clay = 17) %>%
  mutate(silt = 15) %>%
  mutate(sand = 68)

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[3],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[3],] %>%
  mutate(clay = 16.3) %>%
  mutate(silt = 15.5) %>%
  mutate(sand = 68.2)

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[4],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[4],] %>%
  mutate(clay = 15.5) %>%
  mutate(silt = 15.7) %>%
  mutate(sand = 68.8)

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[5],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[5],] %>%
  mutate(clay = round(clay.content.100.200/10, digits = 1)) %>%
  mutate(silt = round(silt.100.200/10, digits = 1)) %>%
  mutate(sand = round(sand.100.200/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[6],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[6],] %>%
  mutate(clay = round(clay.content.100.200/10, digits = 1)) %>%
  mutate(silt = round(silt.100.200/10, digits = 1)) %>%
  mutate(sand = round(sand.100.200/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[7],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[7],] %>%
  mutate(clay = round(clay.content.100.200/10, digits = 1)) %>%
  mutate(silt = 17.5) %>%
  mutate(sand = round(sand.100.200/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[8],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[8],] %>%
  mutate(clay = round(clay.content.100.200/10, digits = 1)) %>%
  mutate(silt = round(silt.100.200/10, digits = 1)) %>%
  mutate(sand = round(sand.100.200/10, digits = 1))

fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[9],] <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[remove[9],] %>%
  mutate(clay = 22.3) %>%
  mutate(silt = 17.3) %>%
  mutate(sand = 60.4)



#fixed_field_data_processed_trees_soils.texture.triangle.100.200 <- fixed_field_data_processed_trees_soils.texture.triangle.100.200[-c(remove),]


plot.100.200 <- 
  ggsoiltexture(fixed_field_data_processed_trees_soils.texture.triangle.100.200, class = "USDA") +
  geom_point(aes(color = Locality), size = 1) +
  scale_color_discrete(type = "viridis")+
  labs(title = "Soil Texture 100-200 cm")

plot.100.200

#### Comparing the means of the explanatory variables ####

#creating a dataframe with the explanatory variables across all the populations for the comparison of means
means_fixed_field_data_processed_terrain_dist <- tibble(Elevation..m.FIXED = c(LM_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                               LC_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                               SD_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        slope_raster_15_data_pts = c(LM_fixed_field_data_processed_terrain_dist$LM_slope_raster_15_data_pts[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                                     LC_fixed_field_data_processed_terrain_dist$LC_slope_raster_15_data_pts[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                                     SD_fixed_field_data_processed_terrain_dist$SD_slope_raster_15_data_pts[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        Eastness = c(LM_fixed_field_data_processed_terrain_dist$LM_Eastness[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                     LC_fixed_field_data_processed_terrain_dist$LC_Eastness[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                     SD_fixed_field_data_processed_terrain_dist$SD_Eastness[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        Northness = c(LM_fixed_field_data_processed_terrain_dist$LM_Northness[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                      LC_fixed_field_data_processed_terrain_dist$LC_Northness[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                      SD_fixed_field_data_processed_terrain_dist$SD_Northness[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        TWI_values = c(LM_fixed_field_data_processed_terrain_dist$LM_TWI_values[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                       LC_fixed_field_data_processed_terrain_dist$LC_TWI_values[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                       SD_fixed_field_data_processed_terrain_dist$SD_TWI_values[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        heat.load = c(LM_fixed_field_data_processed_terrain_dist$heat.load[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                      LC_fixed_field_data_processed_terrain_dist$heat.load[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                      SD_fixed_field_data_processed_terrain_dist$heat.load[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        d = c(LM_fixed_field_data_processed_terrain_dist$d[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                              LC_fixed_field_data_processed_terrain_dist$d[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                              SD_fixed_field_data_processed_terrain_dist$d[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]),
                                                        Locality = c(LM_fixed_field_data_processed_terrain_dist$Locality[LM_fixed_field_data_processed_terrain_dist$Locality == "LM"],
                                                                     LC_fixed_field_data_processed_terrain_dist$Locality[LC_fixed_field_data_processed_terrain_dist$Locality == "LC"],
                                                                     SD_fixed_field_data_processed_terrain_dist$Locality[SD_fixed_field_data_processed_terrain_dist$Locality == "SD"]))


#elevation
ggplot()+
  labs(title = "Elevation Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  # geom_boxplot(data = LM_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  # geom_boxplot(data = LC_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  # geom_boxplot(data = SD_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(Elevation..m.FIXED ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(Elevation..m.FIXED ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc



#slope
ggplot()+
  labs(title = "Slope Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, slope_raster_15_data_pts), fill = "light blue")+
  # geom_boxplot(data = LM_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  # geom_boxplot(data = LC_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  # geom_boxplot(data = SD_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(slope_raster_15_data_pts ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(slope_raster_15_data_pts ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$slope_raster_15_data_pts, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc


#Eastness

ggplot()+
  labs(title = "Eastness Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, Eastness), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(Eastness ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(Eastness ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$Eastness, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc

#northness
ggplot()+
  labs(title = "Northness Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, Northness), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(Northness ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(Northness ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$Northness, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc

#TWI
ggplot()+
  labs(title = "TWI Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, TWI_values), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(TWI_values ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(TWI_values ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$TWI_values, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc


#distance

ggplot()+
  labs(title = "Distance Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, d), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(d ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(d ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$d, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc


#Heat load index

ggplot()+
  labs(title = "HLI Across Population")+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, heat.load), fill = "light blue")+
  theme_minimal()

# checking to see if residuals are normal
anova_available_water_0.5 <- aov(heat.load ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
anova_available_water_0.5
hist(anova_available_water_0.5$residuals, xlab = "Residuals", main = "Distribution of Residuals for Clay/Loam Available Water at 0-5 cm vs. Population")

#kruskal-Wallis test because the data does not have normally-distributed residuals and equal variance of residuals
test <- kruskal.test(heat.load ~ Locality, data = means_fixed_field_data_processed_terrain_dist)
test

#post-hoc Wilcoxon rank sum tests
post_hoc <- pairwise.wilcox.test(means_fixed_field_data_processed_terrain_dist$heat.load, means_fixed_field_data_processed_terrain_dist$Locality,
                                 p.adjust.method = "fdr") #p value adjusted using false discovery rate method
post_hoc



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

