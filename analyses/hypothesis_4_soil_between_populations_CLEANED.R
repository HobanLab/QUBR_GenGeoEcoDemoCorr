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
source("./analyses/Data_Processing_Script.R")

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
  filter(!is.na(tree_row_num)) #filters out the grids without trees inside of them

#filtering out point data to be just the trees within the grids
LM_fixed_field_data_processed_trees_soils <- LM_fixed_field_data_processed_soils %>%
  filter(X %in% LM_list_grids_and_trees_fixed$tree_row_num)  #creating a dataframe with the row numbers that match between the overall tree points dataframe and the focal tree points dataframe

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
        # We also performed this test for ever soil metric because it is a non-parametric test allowing for comparisons across the soil metrics. 
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

#post hoc test of Welch's ANOVA test: Tamhane's T2 Test
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

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_silt_0.5$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
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

#post-hoc Pairwise T-Tests
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

#post-hoc Pairwise T-Tests
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

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_sand_100.200$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
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
mean_soil_function_ph_100_200 <- mean_soil_function("ph_0.5")

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

#CHOSEN TEST: kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_vol_water_10_100_200$final_test

#kruskal-Wallis test post-hoc Wilcoxon rank sum tests
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

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_33_0_5$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
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

#CHOSEN TEST: Welch's ANOVA, does not assume equal variances 
mean_soil_function_vol_water_1500_0_5$final_test

#post hoc Welch's ANOVA test: Tamhane's T2 Test
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

#CHOSEN TEST: Kruskal-Wallis test, non-parametric for non-normal residuals and non-equal variance
mean_soil_function_sandy_avail_water_100_200$final_test

#Kruskal-Wallis test post-hoc Wilcoxon rank sum tests
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

#CHOSEN TEST: ANOVA test, assumes equal variance and normal residuals
mean_soil_function_clay_loam_avail_water_100_200$final_test

#post-hoc Pairwise T-Tests
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






