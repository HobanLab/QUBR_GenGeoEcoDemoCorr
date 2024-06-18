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

#downloading fixed_field_data_processed as csv to help import into google maps and google earth
write.csv(fixed_field_data_processed, "~/Morton Arboretum REU 2024/Untitled/QUBR_GenGeoEcoDemoCorr/analyses/fixed_field_data_processed.csv")

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

#turning the shapefile into a spatial object and then an owin (a window)
BCS_polygon <- as_Spatial(BCS_polygon) #converts the shapefile into spatial object 
W <- owin(BCS_polygon) # convert the BCS sp into the owin class and setting the x and y ranges, the "bbox"
plot(W)

BCS_polygon$co
View(BCS_polygon)


#creating the ppp for the entire extent of points
BCS_ppp <- ppp(x = fixed_field_data_processed$long, y = fixed_field_data_processed$lat, window = W) #creating the poisson point pattern
plot(BCS_ppp)

#creating the ppp for LM
LM_ppp <- ppp(x = LM_fixed_field_data_processed$long, y = LM_fixed_field_data_processed$lat, window = W) #creating the poisson point pattern for lm
plot(LM_ppp)

#creating the ppp for LC
LC_ppp <- ppp(x = LC_fixed_field_data_processed$long, y = LC_fixed_field_data_processed$lat, window = W) #creating the poisson point pattern for lm
plot(LC_ppp)

#creating the ppp for SD
SD_ppp <- ppp(x = SD_fixed_field_data_processed$long, y = SD_fixed_field_data_processed$lat, window = W) #creating the poisson point pattern for lm
plot(SD_ppp)

#### Ripley's K ####

#running the Ripley's K analysis for all points
K_BCS <- Kest(BCS_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_BCS <- Kest(BCS_ppp) #includes the edge corrections
plot(K_BCS, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_BCS, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot

#running the Ripley's K analysis for LM
K_LM <- Kest(LM_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_LM <- Kest(LM_ppp) #includes the edge corrections
plot(K_LM, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_LM, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot

#running the Ripley's K analysis for LC
K_LC <- Kest(LC_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_LC <- Kest(LC_ppp) #includes the edge corrections
plot(K_LC, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_LC, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot

#running the Ripley's K analysis for SD
K_SD <- Kest(SD_ppp, correction = "Ripley") #focuses on the K poisson value, the Ripley's K
K_SD <- Kest(SD_ppp) #includes the edge corrections
plot(K_SD, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) #legend inside of the plot
plot(K_SD, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) )) #legend outside of the plot


#### Ripley's L ####

#Ripley's L
L <- Lest(BCS_ppp, main=NULL)
L <- Lest(BCS_ppp, main=NULL, correction = "Ripley")
plot(L, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE, inset=c(1.01, 0) ))


#### ANN Analysis (test for clustering/dispersion) ####

ann.p


