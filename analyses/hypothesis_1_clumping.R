#### Loading libraries and relevant data####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing

# Read in the raw data from the Field_datasheets_filled_before_KA_check_copy.csv
field_data_raw <- read.csv("./data/Field_datasheets_filled_before_KA_check_copy.csv", na.strings = c("NA", "")) 

#for this analysis we will be using the dataframe 
#created in analyzing_morpho_data_cleaned: fixed_field_data_processed

View(fixed_field_data_processed)

#### Ripley's K ####

#ripley(fixed_field_data_processed[fixed_field_data_processed$Locality == "LM"], "lat", "long")

K <- Kest(fixed_field_data_processed)
plot(K, main=NULL, las=1, legendargs=list(cex=0.8, xpd=TRUE)) 


