#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(spdep) # to use lag.listw
library(ape) # for computing the Moran's I stat

fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

#transforming the shapefile of trees from WGS84 into equal area projection UTM 12N
fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

#### Computing Average Nearest Neighbors for each tree ####

#create dataframe with X and Y UTM coordinates
fixed_field_data_processed_sf_trans_coords <- st_coordinates(fixed_field_data_processed_sf_transformed) #creates a dataframe with seperate x and y columns from the UTM 12N transformation
fixed_field_data_processed_sf_trans_coordinates <- fixed_field_data_processed_sf_transformed %>%
  cbind(fixed_field_data_processed_sf_trans_coords) #combines the x and y coordinate data frame with the transformed sf dataframe

#add average nearest neighbor for each individual column
fixed_field_data_processed_NN_UTM <- fixed_field_data_processed_sf_trans_coordinates %>%  #creates a dataframe with the ANN of the closest 5 individual trees for each individual
  mutate(dist1 = nndist(X = X.1, Y= Y, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
  mutate(dist2 = nndist(X = X.1, Y= Y, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
  mutate(dist3 = nndist(X = X.1, Y= Y, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
  mutate(dist4 = nndist(X = X.1, Y= Y, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
  mutate(dist5 = nndist(X = X.1, Y= Y, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
  rowwise()%>% #so that in the next part we take the averages across rows
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5))) %>% #creates a column of the average distances (1-5) of each individual
  dplyr::select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances

mean(c(1.405577,3.354128,8.840866,25.245919,25.470333))
View(fixed_field_data_processed_NN_UTM)


#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN_UTM %>%
  filter(Locality == "SD")

#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN_UTM) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)

#checking for duplicates
duplicates <- fixed_field_data_processed_NN_UTM %>% #creates a dataframe called duplicates that filters out X.1 and Y if they have duplicates
  filter(duplicated(X.1) == TRUE) %>%
  filter(duplicated(Y) == TRUE)
View(duplicates)
which(duplicated(fixed_field_data_processed_NN_UTM$X.1)) #finds which rows have duplicates and returns the row number
which(duplicated(fixed_field_data_processed_NN_UTM$Y))  #finds which rows have duplicates and returns the row number

#### Linear Model ####

#creating columns with transformations: logged all of the variables
fixed_field_data_processed_NN_UTM_log <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = log(Canopy_short))%>%
  mutate(Canopy_long_lg = log(Canopy_long))%>%
  mutate(Canopy_area_lg = log(Canopy_area))%>%
  mutate(Crown_spread_lg = log(Crown_spread))%>%
  mutate(DBH_ag_lg = log(DBH_ag))%>%
  mutate(ANN = log(ANN))
View(fixed_field_data_processed_NN_UTM_log)

#creating columns with transformations: square root all of the variables
fixed_field_data_processed_NN_UTM_sqrt <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = sqrt(Canopy_short))%>%
  mutate(Canopy_long_lg = sqrt(Canopy_long))%>%
  mutate(Canopy_area_lg = sqrt(Canopy_area))%>%
  mutate(Crown_spread_lg = sqrt(Crown_spread))%>%
  mutate(DBH_ag_lg = sqrt(DBH_ag))%>%
  mutate(ANN = sqrt(ANN))
View(fixed_field_data_processed_NN_UTM_sqrt)

#creating columns with transformations: square root all of the variables
fixed_field_data_processed_NN_UTM_inverse <- fixed_field_data_processed_NN_UTM %>%
  mutate(Canopy_short_lg = 1/(Canopy_short))%>%
  mutate(Canopy_long_lg = 1/(Canopy_long))%>%
  mutate(Canopy_area_lg = 1/(Canopy_area))%>%
  mutate(Crown_spread_lg = 1/(Crown_spread))%>%
  mutate(DBH_ag_lg = 1/(DBH_ag))%>%
  mutate(ANN = 1/(ANN))
View(fixed_field_data_processed_NN_UTM_inverse)


fixed_field_data_processed_NN_UTM$


#Linear Model for all points

#conditions are lINES: linearity, independence, normal distribution of residuals, equal variance, simple random sample

#checking linearity 

#plotting the linear model in ggplot for SCA, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN_UTM, (aes(x=ANN, y=Canopy_short)))+ 
  geom_smooth(method='lm')+
  geom_point()+
  xlab("ANN")+
  ylab("Short Canopy Axis")

#creating the linear regression
lm_ANN_Canopy_short <- lm(fixed_field_data_processed_NN_UTM_log$Canopy_short ~ fixed_field_data_processed_NN_UTM_log$ANN)


#checking normality of residuals with a histogram and qqnorm plot
ggplot(lm_ANN_Canopy_short, aes(x= lm_ANN_Canopy_short$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Short Canopy Axis vs. ANN")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(lm_ANN_Canopy_short, aes(sample = lm_ANN_Canopy_short$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = lm_ANN_Canopy_short, aes(x = lm_ANN_Canopy_short$fitted.values, y = lm_ANN_Canopy_short$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for ANN and SCA")
  
#Slope Test visible in summary of the lm
summary(lm_ANN_Canopy_short)

#plotting the linear model in ggplot for LCA

#checking linearity 

#plotting the linear model in ggplot for LCA, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN_UTM, (aes(x=ANN, y=Canopy_long)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("ANN")+
  ylab("Long Canopy Axis")

#creating the linear regression
lm_ANN_Canopy_long <- lm(fixed_field_data_processed_NN_UTM$Canopy_long ~ fixed_field_data_processed_NN_UTM$ANN)

#checking normality of residuals with the histogram and a qqnorm plot
ggplot(lm_ANN_Canopy_long, aes(x= lm_ANN_Canopy_long$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for ANN vs. Long Canopy Axis")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(lm_ANN_Canopy_long, aes(sample = lm_ANN_Canopy_long$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = lm_ANN_Canopy_long, aes(x = lm_ANN_Canopy_long$fitted.values, y = lm_ANN_Canopy_long$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for ANN and LCA")

#Slope Test visible in summary of the lm
summary(lm_ANN_Canopy_long)

#plotting the linear model in ggplot for CA

#checking linearity 

#plotting the linear model in ggplot for LCA, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN_UTM, (aes(x=ANN, y=Canopy_area)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("ANN")+
  ylab("Canopy Area")

#creating the linear regression
lm_ANN_Canopy_Area <- lm(fixed_field_data_processed_NN_UTM$Canopy_area ~ fixed_field_data_processed_NN_UTM$ANN)

#checking normality of residuals with a histogram and qqnorm plot
ggplot(lm_ANN_Canopy_Area, aes(x= lm_ANN_Canopy_Area$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for ANN vs. Canopy Area")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(lm_ANN_Canopy_Area, aes(sample = lm_ANN_Canopy_Area$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = lm_ANN_Canopy_Area, aes(x = lm_ANN_Canopy_Area$fitted.values, y = lm_ANN_Canopy_Area$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for ANN and CA")

#Slope Test visible in summary of the LM
summary(lm_ANN_Canopy_Area)

#plotting the linear model in ggplot for CS

#checking linearity 

#plotting the linear model in ggplot for CS, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN_UTM, (aes(x=ANN, y=Crown_spread)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("ANN")+
  ylab("Crown Spread")

#creating the linear regression
lm_ANN_Crown_Spread <- lm(fixed_field_data_processed_NN_UTM$ANN ~ fixed_field_data_processed_NN_UTM$Crown_spread)

#checking normality of residuals with histograms and qq norm plots
ggplot(lm_ANN_Crown_Spread, aes(x= lm_ANN_Crown_Spread$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for ANN vs. Crown Spread")+
  xlab("Residuals")+
  ylab("Frequency")


ggplot(lm_ANN_Crown_Spread, aes(sample = lm_ANN_Crown_Spread$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = lm_ANN_Crown_Spread, aes(x = lm_ANN_Crown_Spread$fitted.values, y = lm_ANN_Crown_Spread$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for ANN and CS")

#Slope Test visible in summary of the lm
summary(lm_ANN_Crown_Spread)

#plotting the linear model in ggplot for DBH_ag

#checking linearity 

#plotting the linear model in ggplot for DBH_ag, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN_UTM, (aes(x=ANN, y=DBH_ag)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("ANN")+
  ylab("DBH_ag")

#creating the linear regression
lm_ANN_DBH_ag <- lm(fixed_field_data_processed_NN_UTM$ANN ~ fixed_field_data_processed_NN_UTM$DBH_ag)

#checking normality of residuals with histogram and qq norm plot
ggplot(lm_ANN_DBH_ag, aes(x= lm_ANN_DBH_ag$residuals))+
  geom_histogram()+
  labs(title = "Distribution of Residuals for Aggregated ANN vs. DBH")+
  xlab("Residuals")+
  ylab("Frequency")

ggplot(lm_ANN_DBH_ag, aes(sample = lm_ANN_DBH_ag$residuals))+
  geom_qq()

#checking equal variance
ggplot(data = lm_ANN_DBH_ag, aes(x = lm_ANN_DBH_ag$fitted.values, y = lm_ANN_DBH_ag$residuals))+
  geom_point()+
  geom_abline(intercept = 0, slope = 0)+
  xlab("Fitted Values")+
  ylab("Residuals")+
  labs(title = "Residuals vs. Fitted Values for ANN and Aggregated DBH")

#Slope Test visible in summary of the lm
summary(lm_ANN_DBH_ag)


#### Moran's I ####

tree.dist <- as.matrix(dist(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))) #making a matrix of the distances between trees
View(tree.dist)
tree.dist.inv <- 1/tree.dist #makes it so closer trees are higher in the matrix
diag(tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
View(tree.dist.inv)

tree.dist.inv[is.infinite(tree.dist.inv)] <- 0 
Moran.I(fixed_field_data_processed_NN_UTM$Canopy_short, tree.dist.inv)

#Monte Carlo

MC<- moran.mc(fixed_field_data_processed_NN_UTM$Canopy_short, tree.dist.inv, nsim = 999)
