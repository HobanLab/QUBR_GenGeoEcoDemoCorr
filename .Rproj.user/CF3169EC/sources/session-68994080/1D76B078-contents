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


fixed_field_data_processed <- read.csv("./analyses/fixed_field_data_processed.csv") #imports the csv created from analyzing_morpho_data_cleaned.R

#transforming the data into shapefiles with either WGS84 or equal area projection UTM 12N
fixed_field_data_processed_sf <- st_as_sf(fixed_field_data_processed, 
                                          coords = c("long", "lat"), crs = 4326)

fixed_field_data_processed_sf_transformed <- st_transform(fixed_field_data_processed_sf, crs = 26912) # this in UTM 12 N an equal area projection

#### Computing Average Nearest Neighbors for each tree ####

#add average nearest neighbor for each individual column
fixed_field_data_processed_NN <- fixed_field_data_processed %>%
  mutate(dist1 = nndist(X = long, Y= lat, k = 1))%>% #creates column for the distances of each tree to their 1st nearest neighbor
  mutate(dist2 = nndist(X = long, Y= lat, k = 2)) %>% #creates column for the distances of each tree to their 2nd nearest neighbor
  mutate(dist3 = nndist(X = long, Y= lat, k = 3)) %>% #creates column for the distances of each tree to their 3rd nearest neighbor
  mutate(dist4 = nndist(X = long, Y= lat, k = 4)) %>% #creates column for the distances of each tree to their 4th nearest neighbor
  mutate(dist5 = nndist(X = long, Y= lat, k = 5)) %>% #creates column for the distances of each tree to their 5th nearest neighbor
  rowwise()%>% #so that in the next part we take the averages across rows
  mutate(ANN = mean(c(dist1, dist2, dist3, dist4, dist5)))%>% #creates a column of the average distances (1-5) of each individual
  select(!c(dist1, dist2, dist3, dist4, dist5)) #removes the excess columns with the 5 nearest neighbor distances

mean(c(1.374369e-05,3.064129e-05,8.661858e-05,2.288619e-04, 2.312526e-04)) #Average nearest neighbors for first row
View(fixed_field_data_processed_NN)


#### Creating fixed_field_data_processed dataframes for each population with the nearest neighbor columns ####

LM_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed_NN %>%
  filter(Locality == "SD")

#### Descriptive Summary ####

#histograms
ggplot(fixed_field_data_processed_NN) + # Generate the base plot
  geom_histogram(aes(x = Canopy_short))+
  xlab("Short Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN) + # Generate the base plot
  geom_histogram(aes(x = Canopy_long))+
  xlab("Long Canopy Axis")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN) + # Generate the base plot
  geom_histogram(aes(x = Crown_spread))+
  xlab("Canopy Spread")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN) + # Generate the base plot
  geom_histogram(aes(x = Canopy_area))+
  xlab("Canopy Area")+
  ylab("Frequency")

ggplot(fixed_field_data_processed_NN) + # Generate the base plot
  geom_histogram(aes(x = DBH_ag))+
  xlab("Aggregated DBH")+
  ylab("Frequency")

hist(fixed_field_data_processed_NN$Canopy_short, main = "Distribution of Short Canopy Axis")
hist(fixed_field_data_processed_NN$Canopy_long, main = "Distribution of Long Canopy Axis")
hist(fixed_field_data_processed_NN$Crown_spread, main = "Distribution of Canopy Spread")
hist(fixed_field_data_processed_NN$Canopy_area, main = "Distribution of Canopy Area")
hist(fixed_field_data_processed_NN$DBH_ag, main = "Distribution of DBH") # slight tail

#Summaries
# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- fixed_field_data_processed_NN %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, DBH_ag) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots
View(field_data_summarized)

#### Linear Model ####

#Linear Model for all points

#conditions are lINES: linearity, independence, normal distribution of residuals, equal variance, simple random sample

#checking linearity 

#plotting the linear model in ggplot for SCA, lineaerity condition is not well met
ggplot(data = fixed_field_data_processed_NN, (aes(x=Canopy_short, y=ANN)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Short Canopy Axis")+
  ylab("ANN")

plot(fixed_field_data_processed_NN$Canopy_short, fixed_field_data_processed_NN$ANN, xlab = "Short Canopy Axis", ylab = "ANN")
lm_LM_ANN <- lm(fixed_field_data_processed_NN$ANN ~ fixed_field_data_processed_NN$Canopy_short)
abline(lm_LM_ANN)


#checking residuals#checking normality
ggplot(lm_LM_ANN, aes(x= lm_LM_ANN$residuals))+
  geom_histogram()+
  labs("Distribution of Residuals for Short Canopy Axis vs. ANN")+
  xlab("Residuals")+
  ylab("Frequency")

qqnorm(lm_LM_ANN$residuals)

#checking equal variance
ggplot(data = lm_LM_ANN, aes(x = fitted.values, y = residuals))

plot(lm_LM_ANN$residuals ~ lm_LM_ANN$fitted.values, 
     xlab = "Fitted Values",
     ylab = "Residuals",
     main = "Residuals vs. Fitted Values")
abline(0,0)


#plotting the linear model in ggplot for LCA
ggplot(data = fixed_field_data_processed_NN, (aes(x=Canopy_short, y=ANN)))+
  geom_smooth(method='lm')+
  geom_point()+
  xlab("Short Canopy Axis")+
  ylab("ANN")

plot(fixed_field_data_processed_NN$Canopy_short, fixed_field_data_processed_NN$ANN, xlab = "Short Canopy Axis", ylab = "ANN")
lm_LM_ANN <- lm(fixed_field_data_processed_NN$ANN ~ fixed_field_data_processed_NN$Canopy_short)
abline(lm_LM_ANN)