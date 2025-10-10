# QUBR_GenGeoEcoDemoCorr

#### Authors: Rebecca Wanger and Ash Hamilton

Looking for correlations in the genetics, geography, ecology, and demographics of *Quercus brandegeei*

*** 

 The [data](./data) folder contains all data files used in all scripts in this project that are cleaned and well commented. 

 The [analyses](./analyses) folder contains all the R scripts split by unique analysis type used in this project (all done in R version 2025.05.1+513)

* [analyzing_moprho_data_cleaned.R](./data/analyzing_moprho_data_cleaned.R) is a well-commented R script which has been modified only slightly from the code Ash wrote for their morphometrics project in the Winter of 2024

**Hypothesis 1:** Quercus brandegeei seeds are predominantly distributed either by heavy rainfall or gravity, 
which may impact both the spatial distribution and the genetic structure of a given population. If they were distributed by 
heavy rainfall, we would expect the trees to be more dispersed than at random. If they were distributed by gravity, we would expect more clumping 
of trees than at random. 

* [hypothesis_1_clumping_CLEANED.R](./analyses/hypothesis_1_clumping_CLEANED.R) is a well-commented R script to 
test hypothesis 1. To test this, we used **Ripley's Ks** analyses across all trees and for each population 
to determine whether the tree points were more clustered or dispersed than we would expect at 
random. Across all and for each population, we used three different types of windows (three Ripley's K test for each) to generate random point distributions inside, 
using either a bounding box around the known tree points, a convex hull around the known tree points, or a buffer around the river shapefiles (100 m for LM and LC and 70 for SD).
We then used an **Average Nearest Neighbor Analysis (ANN)** to support with p-values whether the points seem more clustered or 
dispersed than at random. We used three different types of windows for random point generation (three ANN tests for each),
using just a river raster outline, using an inverse distance to river outline, and using an on and inside the river raster to control
for the potential varying influences of the river on tree clustering/dispersal. Finally, we used **Poisson Point Model Analyses**, to see if models that take into account the effect that the trees' inverse 
distance to the river has on the placement of the trees better explains the distribution of the points than if they were 
distributed at random. We observed significant clustering of trees and significant influence of the rivers on the distribution of trees at and across each population.

  **Hypothesis 2:** The size and shape of *Quercus brandegeei* individuals across all sites is impacted by distance to other individuals of the 
  same species either due to competition or facilitation.  The size and shape 
of Quercus brandegeei individuals across all sites is impacted by distance to other individuals of the same species either due to 
competition or facilitation. If they are impacted by facilitation, we would expect closer trees would be bigger. If they are impacted 
by competition, we would expect closer trees to be smaller. This hypothesis is explored by three scripts:
    - hypothesis_2_moran_CLEANED.R
    - hypothesis_2_linear_model_CLEANED.R
    - hypothesis_2_linear_model_no_outliers_CLEANED.R 


* [hypothesis_2_moran_CLEANED.R](./analyses/hypothesis_2_moran_CLEANED.R) is a well-commented R script to test hypothesis 2 with Global and Local Moran's I to determine whether values 
of Short Canopy Axis, Long Canopy Axis, Crown Spread, Canopy Area, and DBH that were closer together were more similar in value or not. The global Moran's I looked for general spatial
autocorrelation, and the local Moran's I looked for areas were values were more similar than other areas. We observed that there is significant global spatial autocorrelation
at all combinations of populations and size and shape metrics, except for La Cobriza and San Dionisio with DBH. We also found varying levels of local spatial autocorrelation with 
every combination of populations and size and shape metrics, except for Las Matancitas with Crown Spread and DBH and San Dionisio with Short Canopy Axis which showed none.

* [hypothesis_2_linear_model_CLEANED.R](./analyses/hypothesis_2_linear_model_CLEANED.R) is well-commented R script to test this hypothesis 2 with generalized linear 
regressions to see if for focal trees (a subset of trees that are independent and randomly selected), there is a relationship between the size of the trees and how much competition they face 
(based on the size and distance of the neighbors from the focal trees). A negative relationship would suggest that there is competition between the trees. A positive relationship would suggest there 
is facilitation between the trees. We tested the usefulness of models with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratics) and checked they met the best 
model (based on AIC) met the conditions. If they did, we used a slope test, and if they did not, we used a non-parametric Kendall's Tau Test to check for any relationships. We observed that 
only San Dionisio with Canopy Area and competition (with no control for spatial autocorrelation or outliers) showed a significant negative relationship, demonstrating potential competition. 
We observed that only in San Dionisio's canopy area showed significant competition between trees.

* [hypothesis_2_linear_model_no_outliers_CLEANED.R](./analyses/hypothesis_2_linear_model_no_outliers_CLEANED.R) is well-commented R script to test hypothesis 2 with generalized linear 
regressions that uses the same methodologies as for [hypothesis_2_linear_model_CLEANED.R](./analyses/hypothesis_2_linear_model_CLEANED.R), except the outliers are removed from the dataframe when the Cook's D
values are greater than three times the mean Cook's D value. Removing the outliers allows the models to better meet the conditions but sacrifices potentially informative trends in the data.
We observed that only San Dionisio's short canopy axis, long canopy axis, canopy area, and crown spread showed significant competition between trees. 

* [hypothesis_2_linear_models_permutations_CLEANED.R](./analyses/hypothesis_2_linear_models_permutations_CLEANED) is well-commented R script to test hypothesis 2 with generalized linear 
regressions that uses the same methodologies as for [hypothesis_2_linear_model_no_outliers_CLEANED.R](./analyses/hypothesis_2_linear_model_no_outliers_CLEANED.R), except we permutate the randomly selected focal trees (500 times)
to test how robust are the slope test/Kruskal-Wallis test results. We observed that across all of the combinations of size and shape metrics had mean slopes demonstrating competition,
but only San Dionisio's short canopy axis and long canopy axis median p-values showed significant competition.   

**Hypothesis 3:** The size and shape of Quercus brandegeei individuals is affected by the frequency of high wind and hurricane events which can be proxied by elevation, aspect, and slope. For trees facing more exposure, higher elevations, 
steeper slopes, and south and east facing slopes, we expected them to be smaller. On the other hand, we expected trees facing less exposure, lower 
elevation, flatter slopes, and north and west facing aspects to be larger.
This hypothesis is explored by three scripts:
    - hypothesis_3_size_and_exposure_CLEANED.R
    - hypothesis_3_GAM_size_and_exposure_CLEANED.R 

* [hypothesis_3_size_and_exposure_CLEANED.R](./analyses/hypothesis_3_size_and_exposure_CLEANED.R) is a well commented R script to test hypothesis 3. To test it, we explored using single linear models
to see if there was a relationship between size characteristics (SCA, LCA, CS, CA, DBH) and elevation. We performed a similar linear model see if the size of trees were affected by their 
slope. Finally, we compared the average size values to their aspect (for both N,E,S,W and for N, NW, W, SW, S, SE, E, NE) with ANOVAs/Kruskal-Wallis 
Models. We used 15 m elevation/slope/aspect rasters. We observed that for elevation vs. size/shape metrics, Las Matancitas showed significant positive relationships between short canopy axis,
long canopy axis, canopy area, and crown spread and elevation, La Cobriza showed a significant negative relationship between DBH and elevation, San Dionisio showed significant negative relationships between SCA, CA, CS,
and DBH and elevation, and across all populations, there was a significant positive relationship between DBH and elevation. As for slope vs. size/shape metrics, Las Matancitas showed a significant negative relationship between LCA and slope,
San Dionisio showed significant negative relationships between SCA, LCA, CA, and CS and slope, and across all populations, LCA, CA, and CS had significant negative relationships with slope. As for comparing the 8 and 4 categories of 
aspect with the size/shape of trees, in general, we observed the north and northeast facing aspects tending to have larger mean sizes compared to more southern and western facing aspects. 

* [hypothesis_3_GAM_size_and_exposure_CLEANED.R](./analyses/hypothesis_3_GAM_size_and_exposure_CLEANED.R) is a well-commented R script to test hypothesis 3 using generalized additive models
as a non-parametric and more flexible method to explore relationships between elevation, slope, and aspect and size characteristics (SCA, LCA, CS, CA, DBH), allowing for non-linear relationships. 

**Hypothesis 4:** Quercus brandegeei persist in a similar, narrow niche breadth of soil characteristics.
This hypothesis is explored by three scripts:
    - hypothesis_3_size_and_exposure_CLEANED.R
    - hypothesis_3_GAM_size_and_exposure_CLEANED.R 

* [hypothesis_4_soil_characteristics.R](./analyses/hypothesis_4_soil_characteristics.R) is a well commented R script to test this hypothesis: The 16 known populations have more similar soil characteristics (texture, moisture, Ph, etc.) than areas where we know there are no Q. brandegeei populations. We used 250 m resolution soil rasters from Soil Grids. We predicted that 16 known sites do not have significantly different soil characteristics (texture, moisture, Ph). For this we performed three analyses. The first was a comparison of the mean soil values between our three populations: La Cobriza, San Dionisio, and Las Matancitas using ANOVA and non-parametric ANOVA subsitute tests. The second was to compare the slopes of size values vs. soil characteristics for each population to the slopes of the same comparison but size values were shuffled and randomized (reflected in histogram of slopes). Finally, we compared the mean soil values for our known populations to randomly selected populations.

* [hypothesis_5_water_availability.R](./analyses/hypothesis_5_water_availability.R) is a well commented R script to test the hypothesis is Q. brandegeei is predominantly restricted by water availability? We predicted the closer the individuals are to the river, the larger they would be and the further they were, the smaller they would be. For this, we did linear regressions for each population to see if the trees distance to the river had a relationship with their size.  

  The results are summarized in this google sheet entitled [General_Test_Results_Wanger_REEF_2025](https://docs.google.com/spreadsheets/d/1BemVj7ev1UcTnCs2zXe8JUpJoGbSO78u0YJfAZuehLs/edit?gid=0#gid=0) 

**Unrefined and unused scripts are stored in [Old_Scripts](./analyses/Old_Scripts)**
