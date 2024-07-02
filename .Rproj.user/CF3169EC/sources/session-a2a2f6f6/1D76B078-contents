#### Loading libraries and relevant data ####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the nndist function
library(spdep) # to use morna's I functions like lag.listw
library(ape) # for computing the Moran's I stat
library(raster) #to use point distance

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

#Test for all points

#regular Moran's I

#creating a matrix of distances between trees where the higher values are at the top
tree.dist <- as.matrix(dist(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))) #making a matrix of the distances between trees
tree.dist.inv <- 1/tree.dist #makes it so closer trees are higher in the matrix
diag(tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
tree.dist.inv[is.infinite(tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(fixed_field_data_processed_NN_UTM$Canopy_short, tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

tree.coord.matrix <- as.matrix(cbind(fixed_field_data_processed_NN_UTM$X.1, 
                                  fixed_field_data_processed_NN_UTM$Y))

max(fixed_field_data_processed_NN_UTM$X.1)

#creates nearest neighbor knn using a matrix of the tree coordinates and k = 1, means the distance to the nearbor is conputed only for the nearest one
knn <- knearneigh(tree.coord.matrix, k = 5) #I have playing around with the k, trying to include all or half of the trees for example
#turns knn into neighbors list
nb <- knn2nb(knn, row.names = NULL, sym = FALSE)


#assigning weights to each neighbor, W style assigns weight to be 1/# of neighbors
      #lw <- nb2listw(nb,zero.policy=TRUE, style="W")
#inverse distance weighting with raw distance-based weights without applying any normalisation
lw <- nb2listwdist(nb, fixed_field_data_processed_NN_UTM, type="idw", style="raw", 
                   alpha = 1, dmax = NULL, longlat = NULL, zero.policy=NULL)

#checks the neighbor weights for the first tree
lw$weights[1]

#creating lags, which computes the average neighboring short canopy axis for each tree
fixed_field_data_processed_NN_UTM$lag.canopy.short <- lag.listw(lw, fixed_field_data_processed_NN_UTM$Canopy_short)
# Create a regression model
M <- lm(lag.canopy.short ~ Canopy_short, fixed_field_data_processed_NN_UTM)

# Plot the lagged variable vs. the variable 
ggplot(data=fixed_field_data_processed_NN_UTM, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(fixed_field_data_processed_NN_UTM$Canopy_short, listw = lw, n = length(nb), S0 = Szero(lw))

#assessing statistical significance with a Monte-Carlo simulation
MC<- moran.mc(fixed_field_data_processed_NN_UTM$Canopy_short, lw, nsim = 999)
MC

#plot of simulated Moran's I values against our value
plot(MC, main="", las=1)
MC$p.value

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local <- localmoran_perm(fixed_field_data_processed_NN_UTM$Canopy_short, lw, nsim = 9999, alternative = "greater")
MC_local.df <- as.data.frame(MC_local)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")

#calculating the p-values for each individual tree Moran's I, observed vs. expected
fixed_field_data_processed_NN_UTM$p  <- MC_local.df$`Pr(folded) Sim`

###Test for LM###

#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
LM.tree.dist <- as.matrix(dist(cbind(LM_fixed_field_data_processed$X.1, 
                                  LM_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
LM.tree.dist.inv <- 1/LM.tree.dist #makes it so closer trees are higher in the matrix
diag(LM.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
LM.tree.dist.inv[is.infinite(LM.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_short, LM.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

LM.tree.coord.matrix <- as.matrix(cbind(LM_fixed_field_data_processed$X.1, 
                                     LM_fixed_field_data_processed$Y))

#creates nearest neighbor knn using a matrix of the tree coordinates and k = 1, means the distance to the nearbor is conputed only for the nearest one
knn.LM <- knearneigh(LM.tree.coord.matrix, k = 5, longlat = F) #I have playing around with the k, trying to include all or half of the trees for example
#turns knn into neighbors list
nb.LM <- knn2nb(knn.LM, row.names = NULL, sym = FALSE)

#assigning weights to each neighbor, W style assigns weight to be 1/# of neighbors
#lw <- nb2listw(nb,zero.policy=TRUE, style="W")
#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.LM <- nb2listwdist(nb.LM, LM_fixed_field_data_processed, type="idw", style="W", 
                   alpha = 1, dmax = NULL, longlat = F, zero.policy=NULL) #it converts it to latlon and then takes the distance because that is also more accurate

#checks the neighbor weights for the first tree
sum(lw.LM$weights[[1]])

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.LM, LM_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.LM <- lm(lag.canopy.short ~ Canopy_short, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_short, listw = lw.LM, n = length(nb.LM), S0 = Szero(lw.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.short <- moran.mc(LM_fixed_field_data_processed$Canopy_short, lw.LM, nsim = 999)
MC.LM.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.LM.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.short <- localmoran_perm(LM_fixed_field_data_processed$Canopy_short, lw.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.short.df <- as.data.frame(MC_local.LM.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.short  <- MC_local.LM.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_long, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.LM, LM_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.LM.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_long, listw = lw.LM, n = length(nb.LM), S0 = Szero(lw.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.long <- moran.mc(LM_fixed_field_data_processed$Canopy_long, lw.LM, nsim = 999)
MC.LM.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.LM.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.long <- localmoran_perm(LM_fixed_field_data_processed$Canopy_long, lw.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.long.df <- as.data.frame(MC_local.LM.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.long  <- MC_local.LM.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.long, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Crown_spread, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.LM, LM_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.LM.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Crown_spread, listw = lw.LM, n = length(nb.LM), S0 = Szero(lw.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.crown.spread <- moran.mc(LM_fixed_field_data_processed$Crown_spread, lw.LM, nsim = 999)
MC.LM.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.LM.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.LM.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.crown.spread <- localmoran_perm(LM_fixed_field_data_processed$Crown_spread, lw.LM, nsim = 9999, alternative = "greater")
MC_local.LM.crown.spread.df <- as.data.frame(MC_local.LM.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.crown.spread  <- MC_local.LM.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(LM_fixed_field_data_processed$p.crown.spread, 
                                                                 method = "fdr", n=length(LM_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$Canopy_area, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.LM, LM_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.LM.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$Canopy_area, listw = lw.LM, n = length(nb.LM), S0 = Szero(lw.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.canopy.area <- moran.mc(LM_fixed_field_data_processed$Canopy_area, lw.LM, nsim = 999)
MC.LM.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.LM.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.LM.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.canopy.area <- localmoran_perm(LM_fixed_field_data_processed$Canopy_area, lw.LM, nsim = 9999, alternative = "greater")
MC_local.LM.canopy.area.df <- as.data.frame(MC_local.LM.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.canopy.area <- MC_local.LM.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(LM_fixed_field_data_processed$p.canopy.area, 
                                                                  method = "fdr", n=length(LM_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LM_fixed_field_data_processed$DBH_ag, LM.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LM_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.LM, LM_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.LM.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LM_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LM_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(LM_fixed_field_data_processed$DBH_ag, listw = lw.LM, n = length(nb.LM), S0 = Szero(lw.LM))

#assessing statistical significance with a Monte-Carlo simulation
MC.LM.dbh.ag <- moran.mc(LM_fixed_field_data_processed$DBH_ag, lw.LM, nsim = 999)
MC.LM.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.LM.dbh.ag, main="", las=1, xlab = "DBH")
MC.LM.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LM.dbh.ag <- localmoran_perm(LM_fixed_field_data_processed$DBH_ag, lw.LM, nsim = 9999, alternative = "greater")
MC_local.LM.dbh.ag.df <- as.data.frame(MC_local.LM.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LM.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LM_fixed_field_data_processed$p.dbh.ag <- MC_local.LM.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LM_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(LM_fixed_field_data_processed$p.dbh.ag, 
                                                                 method = "fdr", n=length(LM_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
LM_box <- st_bbox(river_LM_trans)
LM_fixed_field_data_processed <- LM_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LM_trans) +
  geom_sf(data =LM_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LM_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LM_box[1], LM_box[3]), ylim = c(LM_box[2], LM_box[4]))+
  labs(color = "Adjusted P Value for DBH")

###Test for LC###

#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
LC.tree.dist <- as.matrix(dist(cbind(LC_fixed_field_data_processed$X.1, 
                                     LC_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
LC.tree.dist.inv <- 1/LC.tree.dist #makes it so closer trees are higher in the matrix
diag(LC.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
LC.tree.dist.inv[is.infinite(LC.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_short, LC.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

LC.tree.coord.matrix <- as.matrix(cbind(LC_fixed_field_data_processed$X.1, 
                                        LC_fixed_field_data_processed$Y))

#creates nearest neighbor knn using a matrix of the tree coordinates and k = 1, means the distance to the nearbor is conputed only for the nearest one
knn.LC <- knearneigh(LC.tree.coord.matrix, k = 5, longlat = F) #I have playing around with the k, trying to include all or half of the trees for example
#turns knn into neighbors list
nb.LC <- knn2nb(knn.LC, row.names = NULL, sym = FALSE)

#assigning weights to each neighbor, W style assigns weight to be 1/# of neighbors
#lw <- nb2listw(nb,zero.policy=TRUE, style="W")
#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.LC <- nb2listwdist(nb.LC, LC_fixed_field_data_processed, type="idw", style="W", 
                      alpha = 1, dmax = NULL, longlat = F, zero.policy=NULL) #it converts it to latlon and then takes the distance because that is also more accurate

#checks the neighbor weights for the first tree
sum(lw.LC$weights[[1]])

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.LC, LC_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.LC <- lm(lag.canopy.short ~ Canopy_short, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_short, listw = lw.LC, n = length(nb.LC), S0 = Szero(lw.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.short <- moran.mc(LC_fixed_field_data_processed$Canopy_short, lw.LC, nsim = 999)
MC.LC.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.LC.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.short <- localmoran_perm(LC_fixed_field_data_processed$Canopy_short, lw.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.short.df <- as.data.frame(MC_local.LC.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.short  <- MC_local.LC.canopy.short.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_long, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.LC, LC_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.LC.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_long, listw = lw.LC, n = length(nb.LC), S0 = Szero(lw.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.long <- moran.mc(LC_fixed_field_data_processed$Canopy_long, lw.LC, nsim = 999)
MC.LC.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.LC.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.long <- localmoran_perm(LC_fixed_field_data_processed$Canopy_long, lw.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.long.df <- as.data.frame(MC_local.LC.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.long  <- MC_local.LC.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.long, 
                                                                 method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Crown_spread, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.LC, LC_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.LC.crown.spread <- lm(lag.crown.spread ~ Crown_spread, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Crown_spread, listw = lw.LC, n = length(nb.LC), S0 = Szero(lw.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.crown.spread <- moran.mc(LC_fixed_field_data_processed$Crown_spread, lw.LC, nsim = 999)
MC.LC.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.LC.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.LC.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.crown.spread <- localmoran_perm(LC_fixed_field_data_processed$Crown_spread, lw.LC, nsim = 9999, alternative = "greater")
MC_local.LC.crown.spread.df <- as.data.frame(MC_local.LC.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.crown.spread  <- MC_local.LC.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(LC_fixed_field_data_processed$p.crown.spread, 
                                                                  method = "fdr", n=length(LC_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$Canopy_area, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.LC, LC_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.LC.canopy.area <- lm(lag.canopy.area ~ Canopy_area, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$Canopy_area, listw = lw.LC, n = length(nb.LC), S0 = Szero(lw.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.canopy.area <- moran.mc(LC_fixed_field_data_processed$Canopy_area, lw.LC, nsim = 999)
MC.LC.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.LC.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.LC.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.canopy.area <- localmoran_perm(LC_fixed_field_data_processed$Canopy_area, lw.LC, nsim = 9999, alternative = "greater")
MC_local.LC.canopy.area.df <- as.data.frame(MC_local.LC.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.canopy.area <- MC_local.LC.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(LC_fixed_field_data_processed$p.canopy.area, 
                                                                 method = "fdr", n=length(LC_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(LC_fixed_field_data_processed$DBH_ag, LC.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
LC_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.LC, LC_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.LC.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, LC_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=LC_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(LC_fixed_field_data_processed$DBH_ag, listw = lw.LC, n = length(nb.LC), S0 = Szero(lw.LC))

#assessing statistical significance with a Monte-Carlo simulation
MC.LC.dbh.ag <- moran.mc(LC_fixed_field_data_processed$DBH_ag, lw.LC, nsim = 999)
MC.LC.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.LC.dbh.ag, main="", las=1, xlab = "DBH")
MC.LC.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.LC.dbh.ag <- localmoran_perm(LC_fixed_field_data_processed$DBH_ag, lw.LC, nsim = 9999, alternative = "greater")
MC_local.LC.dbh.ag.df <- as.data.frame(MC_local.LC.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.LC.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
LC_fixed_field_data_processed$p.dbh.ag <- MC_local.LC.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
LC_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(LC_fixed_field_data_processed$p.dbh.ag, 
                                                            method = "fdr", n=length(LC_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
LC_box <- st_bbox(river_LC_trans)
LC_fixed_field_data_processed <- LC_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_LC_trans) +
  geom_sf(data =LC_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = LC_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(LC_box[1], LC_box[3]), ylim = c(LC_box[2], LC_box[4]))+
  labs(color = "Adjusted P Value for CA")

#Test for SD


#Short Canopy Axis

#global Moran's I

#creating a matrix of distances between trees where the higher values are at the top
SD.tree.dist <- as.matrix(dist(cbind(SD_fixed_field_data_processed$X.1, 
                                     SD_fixed_field_data_processed$Y))) #making a matrix of the distances between trees
SD.tree.dist.inv <- 1/SD.tree.dist #makes it so closer trees are higher in the matrix
diag(SD.tree.dist.inv) <- 0 #makes so trees have a 0 distance with themselves
SD.tree.dist.inv[is.infinite(SD.tree.dist.inv)] <- 0 # solves problem presented by duplicated GPS points for trees that were very close to one another

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_short, SD.tree.dist.inv)

#Moran's I and Monte Carlo, using Lags, requires package: spdep

SD.tree.coord.matrix <- as.matrix(cbind(SD_fixed_field_data_processed$X.1, 
                                        SD_fixed_field_data_processed$Y))

#creates nearest neighbor knn using a matrix of the tree coordinates and k = 1, means the distance to the nearbor is conputed only for the nearest one
knn.SD <- knearneigh(SD.tree.coord.matrix, k = 5, longlat = F) #I have playing around with the k, trying to include all or half of the trees for example
#turns knn into neighbors list
nb.SD <- knn2nb(knn.SD, row.names = NULL, sym = FALSE)

#assigning weights to each neighbor, W style assigns weight to be 1/# of neighbors
#lw <- nb2listw(nb,zero.policy=TRUE, style="W")
#inverse distance weighting with raw distance-based weights without applying any normalisation
lw.SD <- nb2listwdist(nb.SD, SD_fixed_field_data_processed, type="idw", style="W", 
                      alpha = 1, dmax = NULL, longlat = F, zero.policy=NULL) #it converts it to latlon and then takes the distance because that is also more accurate

#checks the neighbor weights for the first tree
sum(lw.SD$weights[[1]])

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.short <- lag.listw(lw.SD, SD_fixed_field_data_processed$Canopy_short)
# Create a regression model
M.SD <- lm(lag.canopy.short ~ Canopy_short, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_short, y=lag.canopy.short))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Short Canopy Axis")+
  ylab("Lagged Short Canopy Axis")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_short, listw = lw.SD, n = length(nb.SD), S0 = Szero(lw.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.short <- moran.mc(SD_fixed_field_data_processed$Canopy_short, lw.SD, nsim = 999)
MC.SD.canopy.short

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.short, main="", las=1, xlab = "Short Canopy Axis")
MC.SD.canopy.short$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.short <- localmoran_perm(SD_fixed_field_data_processed$Canopy_short, lw.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.short.df <- as.data.frame(MC_local.SD.canopy.short)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.short.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic")+
  ylab("Expected Moran's I Statistic")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.short  <- MC_local.SD.canopy.short.df$`Pr(folded) Sim`
#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.short.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.short, 
                                                                  method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.short))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.short.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.short.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for SCA")

###Long Canopy Axis

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_long, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.long <- lag.listw(lw.SD, SD_fixed_field_data_processed$Canopy_long)
# Create a regression model
M.SD.Canopy.Long <- lm(lag.canopy.long ~ Canopy_long, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_long, y=lag.canopy.long))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Long Canopy Axis")+
  ylab("Lagged Long Canopy Axis")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_long, listw = lw.SD, n = length(nb.SD), S0 = Szero(lw.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.long <- moran.mc(SD_fixed_field_data_processed$Canopy_long, lw.SD, nsim = 999)
MC.SD.canopy.long

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.long, main="", las=1, xlab = "Long Canopy Axis")
MC.SD.canopy.long$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.long <- localmoran_perm(SD_fixed_field_data_processed$Canopy_long, lw.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.long.df <- as.data.frame(MC_local.SD.canopy.long)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.long.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Long Canopy Axis")+
  ylab("Expected Moran's I Statistic for Long Canopy Axis")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.long  <- MC_local.SD.canopy.long.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.long.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.long, 
                                                                 method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.long))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.long.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.long.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for LCA")

###Crown Spread

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Crown_spread, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.crown.spread <- lag.listw(lw.SD, SD_fixed_field_data_processed$Crown_spread)
# Create a regression model
M.SD.crown.spread <- lm(lag.crown.spread ~ Crown_spread, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Crown_spread, y=lag.crown.spread))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Crown Spread")+
  ylab("Lagged Crown Spread")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Crown_spread, listw = lw.SD, n = length(nb.SD), S0 = Szero(lw.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.crown.spread <- moran.mc(SD_fixed_field_data_processed$Crown_spread, lw.SD, nsim = 999)
MC.SD.crown.spread

#plot of simulated Moran's I values against our value
plot(MC.SD.crown.spread, main="", las=1, xlab = "Crown Spread")
MC.SD.crown.spread$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.crown.spread <- localmoran_perm(SD_fixed_field_data_processed$Crown_spread, lw.SD, nsim = 9999, alternative = "greater")
MC_local.SD.crown.spread.df <- as.data.frame(MC_local.SD.crown.spread)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.crown.spread.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Crown Spread")+
  ylab("Expected Moran's I Statistic for Crown Spread")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.crown.spread  <- MC_local.SD.crown.spread.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.crown.spread.adjusted <- p.adjust(SD_fixed_field_data_processed$p.crown.spread, 
                                                                  method = "fdr", n=length(SD_fixed_field_data_processed$p.crown.spread))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.crown.spread.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.crown.spread.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CS")

###Canopy Area

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$Canopy_area, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.canopy.area <- lag.listw(lw.SD, SD_fixed_field_data_processed$Canopy_area)
# Create a regression model
M.SD.canopy.area <- lm(lag.canopy.area ~ Canopy_area, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=Canopy_area, y=lag.canopy.area))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("Canopy Area")+
  ylab("Lagged Canopy Area")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$Canopy_area, listw = lw.SD, n = length(nb.SD), S0 = Szero(lw.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.canopy.area <- moran.mc(SD_fixed_field_data_processed$Canopy_area, lw.SD, nsim = 999)
MC.SD.canopy.area

#plot of simulated Moran's I values against our value
plot(MC.SD.canopy.area, main="", las=1, xlab = "Canopy Area")
MC.SD.canopy.area$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.canopy.area <- localmoran_perm(SD_fixed_field_data_processed$Canopy_area, lw.SD, nsim = 9999, alternative = "greater")
MC_local.SD.canopy.area.df <- as.data.frame(MC_local.SD.canopy.area)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.canopy.area.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for Canopy Area")+
  ylab("Expected Moran's I Statistic for Canopy Area")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.canopy.area <- MC_local.SD.canopy.area.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.canopy.area.adjusted <- p.adjust(SD_fixed_field_data_processed$p.canopy.area, 
                                                                 method = "fdr", n=length(SD_fixed_field_data_processed$p.canopy.area))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.canopy.area.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.canopy.area.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CA")

###Aggregated dbh

#global Moran's I

#computing the Global Moran's I statistic
Moran.I(SD_fixed_field_data_processed$DBH_ag, SD.tree.dist.inv)

#creating lags, which computes the average neighboring short canopy axis for each tree
SD_fixed_field_data_processed$lag.dbh.ag <- lag.listw(lw.SD, SD_fixed_field_data_processed$DBH_ag)
# Create a regression model
M.SD.canopy.area <- lm(lag.dbh.ag ~ DBH_ag, SD_fixed_field_data_processed)

# Plot the lagged variable vs. the variable 
ggplot(data=SD_fixed_field_data_processed, aes(x=DBH_ag, y=lag.dbh.ag))+
  geom_point()+
  geom_smooth(method = lm, col="blue")+
  xlab("DBH")+
  ylab("Lagged DBH")

#computing the Moran's I statistic
moran(SD_fixed_field_data_processed$DBH_ag, listw = lw.SD, n = length(nb.SD), S0 = Szero(lw.SD))

#assessing statistical significance with a Monte-Carlo simulation
MC.SD.dbh.ag <- moran.mc(SD_fixed_field_data_processed$DBH_ag, lw.SD, nsim = 999)
MC.SD.dbh.ag

#plot of simulated Moran's I values against our value
plot(MC.SD.dbh.ag, main="", las=1, xlab = "DBH")
MC.SD.dbh.ag$p.value #extracting the pvalue

#Local Moran's I 

#using the weighted neighbors to simulate size values
MC_local.SD.dbh.ag <- localmoran_perm(SD_fixed_field_data_processed$DBH_ag, lw.SD, nsim = 9999, alternative = "greater")
MC_local.SD.dbh.ag.df <- as.data.frame(MC_local.SD.dbh.ag)

##Ii is local moran statistic, E.Ii is expected local moran statistic, Vari.Ii is variance of local moran statistic, Z. Ii standard deviation of local moran statistic  
#plotting the local moran's I values vs. the expected
ggplot(data=MC_local.SD.dbh.ag.df)+
  geom_point(aes(x=Ii, y=E.Ii), size = 0.01)+
  xlab("Local Moran's I Statistic for DBH")+
  ylab("Expected Moran's I Statistic for DBH")+
  theme_gray()

#calculating the p-values for each individual tree Moran's I, observed vs. expected
SD_fixed_field_data_processed$p.dbh.ag <- MC_local.SD.dbh.ag.df$`Pr(folded) Sim`

#adjusting the p-vlaues to take into account multiple tests
SD_fixed_field_data_processed$p.dbh.ag.adjusted <- p.adjust(SD_fixed_field_data_processed$p.dbh.ag, 
                                                            method = "fdr", n=length(SD_fixed_field_data_processed$p.dbh.ag))

#representing the p-values of the points on a map
SD_box <- st_bbox(river_SD_trans)
SD_fixed_field_data_processed <- SD_fixed_field_data_processed %>%
  mutate(pval_sig = p.dbh.ag.adjusted <= .05)

ggplot() +
  geom_sf(data =river_SD_trans) +
  geom_sf(data =SD_fixed_field_data_processed, aes(color = p.dbh.ag.adjusted)) +
  geom_sf(data = SD_fixed_field_data_processed %>% filter(pval_sig == T), color = "red") +
  coord_sf(xlim = c(SD_box[1], SD_box[3]), ylim = c(SD_box[2], SD_box[4]))+
  labs(color = "Adjusted P Value for CA")

