# QUBR_GenGeoEcoDemoCorr

#### Authors: Rebecca Wanger and Ash Hamilton

Looking for correlations in the genetics, geography, ecology, and demographics of *Quercus brandegeei*

*** 

 The [data](./data) folder contains all data files used in all scripts in this project that are cleaned and well commented. 

 The [analyses](./analyses) folder contains all the R scripts split by unique analysis type used in this project (all done in R version 4.2.1)

* [analyzing_moprho_data_cleaned.R](./data/analyzing_moprho_data_cleaned.R) is a well commented R script which has been modified only slightly from the code Ash wrote for their morphometrics project in the Winter of 2024

* [hypothesis_1_clumping_CLEANED.R](./analyses/hypothesis_1_clumping_CLEANED.R) is a well commented R script to 
test this hypothesis: whether Quercus brandegeei seeds are predominantly distributed either by heavy rainfall or gravity, 
which may impact both the spatial distribution and the genetic structure of a given population. If they were distributed by 
heavy rainfall, we would expect the trees to be more dispersed than at random. If they were distributed by gravity, we would expect more clumping 
of trees than at random. To test this, we used **Ripley's Ks** analyses across all trees and for each population 
to determine whether the tree points were more clustered or dispersed than we would expect at 
random. Across all and for each population, we used three different types of windows (three Ripley's K test for each) to generate random point distributions inside, 
using either a bounding box around the known tree points, a convex hull around the known tree points, or a buffer around the river shapefiles (100 m for LM and LC and 70 for SD).
We then used an **Average Nearest Neighbor Analysis (ANN)** to support with p-values whether the points seem more clustered or 
dispersed than at random. We used three different types of windows for random point generation (three ANN tests for each),
using just a river raster outline, using an inverse distance to river outline, and using an on and inside the river raster to control
for the potential varying influences of the river on tree clustering/dispersal. Finally, we used **Poisson Point Model Analyses**, to see if models that take into account the effect that the trees' inverse 
distance to the river has on the placement of the trees better explains the distribution of the points than if they were 
distributed at random. We observed significant clustering of trees and significant influence of the rivers on the distribution of trees at each population.

* [hypothesis_2_size_and_comp.R](./analyses/hypothesis_2_size_and_comp.R) is well commented R script to test this hypothesis: The size and shape 
of Quercus brandegeei individuals across all sites is impacted by distance to other individuals of the same species either due to 
competition or facilitation. If they are impacted by facilitation, we would expect closer trees would be bigger. If they are impacted 
by competition, we would expect closer trees to be smaller. To test this, we used Global and Local Moran's I to determine whether values 
of SCA, LCA, CS, CA, and DBH that were closer together were more similar in value or not. The global Moran's I looked for general spatial
autocorrelation and the local Moran's I looked for areas were values were more similar than other areas. We used also performed a linear 
regression to see if for focal trees, there was a relationship between how much competition they face (based on the size of the neighbors 
over their distance to the focal trees) and the size of the focal trees. 

* [hypothesis_3_size_and_exposure.R](./analyses/hypothesis_3_size_and_exposure.R) is a well commented R script to test this hypothesis: is the size and shape of Quercus brandegeei individuals affected by the frequency of  high wind and hurricane events. For trees facing more exposure, higher elevations, steeper slopes, and south and east facing slopes, we expected them to be smaller. On the other hand, we expected trees facing less exposure, lower elevation, flatter slopes, and north and west facing aspects to be larger. To test this, we used linear models to see if there was a relationship between size characteristics (SCA, LCA, CS, CA, DBH) and elevation. We performed a similar linear model see if the size of trees were affected by their slope. Finally, we compared the average size values to their aspect (for both N,E,S,W and for N, NW, W, SW, S, SE, E, NE) with ANOVAs/Kruskal-Wallis Models. We used 15 m elevation/slope/aspect rasters. 

* [hypothesis_3_multiple_lm_size_and_exposure.R](./analyses/hypothesis_3_multiple_lm_size_and_exposure.R) is a well commented R script to test this hypothesis the same hypothesis as above: is the size and shape of Quercus brandegeei individuals affected by the frequency of  high wind and hurricane events. For trees facing more exposure, higher elevations, steeper slopes, and south and east facing slopes, we expected them to be smaller. On the other hand, we expected trees facing less exposure, lower elevation, flatter slopes, and north and west facing aspects to be larger. We wanted to see if we could find best fit models for predicting size characteristics (Short canopy axis, long canopy axis, crown spread, canopy area, and DBH) from elevation, aspect, and/or slope through finding the best multiple linear regression model. We created a base model without interactions and one with (which we found through recursive binary partitioning), and then we used dredging, AIC values, and backwards stepwise regression to narrow down each model. We used nested F tets to compare final models to see which ones would be best. We used 15 m elevation/slope/aspect rasters.

* [hypothesis_4_soil_characteristics.R](./analyses/hypothesis_4_soil_characteristics.R) is a well commented R script to test this hypothesis: The 16 known populations have more similar soil characteristics (texture, moisture, Ph, etc.) than areas where we know there are no Q. brandegeei populations. We used 250 m resolution soil rasters from Soil Grids. We predicted that 16 known sites do not have significantly different soil characteristics (texture, moisture, Ph). For this we performed three analyses. The first was a comparison of the mean soil values between our three populations: La Cobriza, San Dionisio, and Las Matancitas using ANOVA and non-parametric ANOVA subsitute tests. The second was to compare the slopes of size values vs. soil characteristics for each population to the slopes of the same comparison but size values were shuffled and randomized (reflected in histogram of slopes). Finally, we compared the mean soil values for our known populations to randomly selected populations.

* [hypothesis_5_water_availability.R](./analyses/hypothesis_5_water_availability.R) is a well commented R script to test the hypothesis is Q. brandegeei is predominantly restricted by water availability? We predicted the closer the individuals are to the river, the larger they would be and the further they were, the smaller they would be. For this, we did linear regressions for each population to see if the trees distance to the river had a relationship with their size.  

* The results are summarized in this google sheet, entiled [General_Test_Results_Wanger_REEF_2025](https://docs.google.com/spreadsheets/d/1BemVj7ev1UcTnCs2zXe8JUpJoGbSO78u0YJfAZuehLs/edit?gid=0#gid=0) 
