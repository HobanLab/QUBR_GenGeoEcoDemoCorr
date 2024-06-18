#for this analysis we will be using the dataframe 
#created in analyzing_morpho_data_cleaned: fixed_field_data_processed
#so to run this file you should first run "analyzing_morpho_data_cleaned"

#### Loading libraries and relevant data####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing
library(spatstat) # to run the Ripley's K function: Kest

#### Creating fixed_field_data_processed dataframes for each population ####

LM_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LC")

SD_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "SD")

#### Importing Shapefile and creating PPP file ####

# read in baja california sur polygon
BCS_polygon <- st_read("bcs_entidad.shp")

#turning the shapefile into a ppp file
BCS_polygon <- as_Spatial(BCS_polygon) #converts the shapefile into spatial object 
W <- owin(BCS_polygon, xrange = c(-115.22376, -109.4132), yrange = c(22.87195, 28.0000)) # convert the BCS sp into the owin class and setting the x and y ranges, the "bbox"
BCS_ppp <- ppp(x = fixed_field_data_processed$long, y = fixed_field_data_processed$lat, window = W) #creating the poisson point pattern
plot(BCS_ppp)

#### Ripley's K ####

#running the Ripley's K analysis
K <- Kest(BCS_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K <- Kest(BCS_ppp) #includes the edge corrections
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot


#### Ripley's L ####

#Ripley's L
L <- Lest(BCS_ppp, main=NULL)
L <- Lest(BCS_ppp, main=NULL, correction = "Ripley")
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))


#### ANN Analysis ####


