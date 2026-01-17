# QUBR_GenGeoEcoDemoCorr

#### Authors: Rebecca Wanger and Ash Hamilton

Looking for correlations in the genetics, geography, ecology, and demographics of *Quercus brandegeei*

*** 

 The [data](./data) folder contains all data files used in all scripts in this project that are cleaned. 

 The [analyses](./analyses) folder contains all the R scripts split by unique analysis type used in this project (all done in R version 2025.05.1+513). The folder contains scripts that are well-commented and cleaned.
The **Unrefined and unused scripts are stored in [Old_Scripts](./analyses/Old_Scripts)**

**[analyzing_moprho_data_cleaned.R](./analyses/analyzing_moprho_data_cleaned.R)** is a well-commented R script which has been modified only slightly from the code Ash wrote for their morphometrics project in the Winter of 2024

**[Data_Processing_Script.R](./analyses/Data_Processing_Script.R)** is a well-commented script for processing/cleaning the
tree spatial, topographic (slope, aspect, and elevation), distance to rivers, and soil metric data. All spatial data/shapefiles were transformed to
UTM 12 N. Shapefiles were created of river buffers (100 m for LM and LC and 70 m for SD), the BCS boundary, and bounding boxes around the individual tree locations
or each population. Field-collected elevation data was corrected for errors. We used INEGI 15 m resolution rasters to generate a finer-resolution elevation column and slope and aspect columns. 
Aspect was recategorized for 4 cardinal directions (N,E,S,W) and 8 directions (N,NE,etc.). We used 250 m resolution soil rasters from [Soil Grids](https://soilgrids.org/) to extract the soil metrics. 
We included 22 soil metrics at 0-5 and 100-200 cm (clay, silt, and sand content, pH, organic carbon density, cation exchange capacity, bulk density, sand/silt field capacity (volume of water at -10 kPa),
clay/loam field capacity (volume of water at -33 kPa), permanent wilting point (volume of water at -1500 kPa), soil organic carbon). We also added four soil metrics: sand available water
and clay/loam available water at 0-5 and 100-200 cm. Available water columns are produced by subtracting field capacity by permanent wilting point. Finally, we processed the spatial, size/shape, and soil data of the 20 known QUBR population locations. We transformed the data to UTM 12 N. We generated a
7,000 km buffer shapefile around the population points in which we cropped the soil metric rasters and extracted the population soil metrics. 

**There are 5 main types of dataframes created in [Data_Processing_Script.R](./analyses/Data_Processing_Script.R) script for the use of the subsequent hypothesis scripts:**

1) *fixed_field_data_processed_sf_trans_coordinates* (the processed spatial/size/shape data for individual trees)
        This dataframe was then filtered for each population, i.e. LM_fixed_field_data_processed
2) *all_points_fixed_field_data_processed_terrain* (the processed spatial, tree size/shape, ELEVATION, SLOPE, and ASPECT for individual trees)
        This dataframe was then filtered for each population, i.e. LM_fixed_field_data_processed_terrain
3) *LM_fixed_field_data_processed_terrain_dist* (the processed spatial, tree size/shape, elevation, slope, aspect, and DISTANCE TO RIVER for individual trees)
        There is only the three dataframes filtered for each population, i.e. LM_fixed_field_data_processed_terrain_dist
4) *LM_fixed_field_data_processed_soils* (the processed spatial, tree size/shape, elevation, slope, aspect, distance to river, and SOIL METRICS for individual trees)
        There is only the three dataframes filtered for each population, i.e. LM_fixed_field_data_processed_soils
5) *all_known_pop_soils* (the processed spatial and soil metrics for the 20 known QUBR population points)


**For the following hypotheses, all of the scripts depend on the data processing script ("Data_Processing_Script.R") being run using the source() function. If the
script was sourced already, it does not need to be run again unless the environment becomes cleared.** 

**Hypothesis 1:** *Quercus brandegeei* seeds are predominantly distributed either by heavy rainfall or gravity, 
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
of *Q. brandegeei* individuals across all sites is impacted by distance to other individuals of the same species either due to 
competition or facilitation. If they are impacted by facilitation, we would expect closer trees would be bigger. If they are impacted 
by competition, we would expect closer trees to be smaller. This hypothesis is explored by four scripts:  



* [hypothesis_2_moran_CLEANED.R](./analyses/hypothesis_2_moran_CLEANED.R) is a well-commented R script to test hypothesis 2 with **Global and Local Moran's I** to determine whether values 
of Short Canopy Axis, Long Canopy Axis, Crown Spread, Canopy Area, and DBH that were closer together were more similar in value or not. The global Moran's I looked for general spatial
autocorrelation, and the local Moran's I looked for areas were values were more similar than other areas. We observed that there is significant global spatial autocorrelation
at all combinations of populations and size and shape metrics, except for La Cobriza and San Dionisio with DBH. We also found varying levels of local spatial autocorrelation with 
every combination of populations and size and shape metrics, except for Las Matancitas with Crown Spread and DBH and San Dionisio with Short Canopy Axis which showed none.

* [hypothesis_2_linear_model_CLEANED.R](./analyses/hypothesis_2_linear_model_CLEANED.R) is well-commented R script to test this hypothesis 2 with **Generalized Linear Regressions (GLRs)** to see if for focal trees (a subset of trees 
that are independent and randomly selected), there is a relationship between the size of the trees and how much competition they face 
(based on the size and distance of the neighbors from the focal trees). A negative relationship would suggest that there is competition between the trees. A positive relationship would suggest there 
is facilitation between the trees. We tested the usefulness of models with different levels of control for spatial autocorrelation (none, exponential, guassian, spherical, linear, rational quadratics) and checked they met the best 
model (based on AIC) met the conditions. If they did, we used a slope test, and if they did not, we used a non-parametric Kendall's Tau Test to check for any relationships. We observed that 
only San Dionisio with Canopy Area and competition (with no control for spatial autocorrelation or outliers) showed a significant negative relationship, demonstrating potential competition. 
We observed that only in San Dionisio's canopy area showed significant competition between trees.

* [hypothesis_2_linear_model_no_outliers_CLEANED.R](./analyses/hypothesis_2_linear_model_no_outliers_CLEANED.R) is well-commented R script to test hypothesis 2 with **Generalized Linear Regressions (GLRs)** that uses the same methodologies
as for [hypothesis_2_linear_model_CLEANED.R](./analyses/hypothesis_2_linear_model_CLEANED.R), except the outliers are removed from the dataframe when the Cook's D
values are greater than three times the mean Cook's D value. Removing the outliers allows the models to better meet the conditions but sacrifices potentially informative trends in the data.
We observed that only San Dionisio's short canopy axis, long canopy axis, canopy area, and crown spread showed significant competition between trees. 

* [hypothesis_2_linear_models_permutations_CLEANED.R](./analyses/hypothesis_2_linear_models_permutations_CLEANED.R) is well-commented R script to test hypothesis 2 with **Generalized Linear Regressions (GLRs)** that 
uses the same methodologies as for [hypothesis_2_linear_model_no_outliers_CLEANED.R](./analyses/hypothesis_2_linear_model_no_outliers_CLEANED.R), except we permutate the randomly selected focal trees (500 times)
to test how robust are the slope test/Kruskal-Wallis test results. We observed that across all of the combinations of size and shape metrics had mean slopes demonstrating competition,
but only San Dionisio's short canopy axis and long canopy axis median p-values showed significant competition.   

**Hypothesis 3:** *Quercus brandegeei* is predominantly restricted by water availability (proxied by distance to river, aspect, slope, and elevation). 
*Q. brandegeei* will exhibit larger sizes with higher soil moisture content. 
This hypothesis is explored by two scripts: 

* [hypothesis_3_size_and_exposure_CLEANED.R](./analyses/hypothesis_3_size_and_water_availability_CLEANED.R) is a well commented R script to test hypothesis 3. To test it, we explored using **Single Linear Regressions (SLRs)**
to see if there was a relationship between size characteristics (SCA, LCA, CS, CA, DBH) and each of the quantiative explanatory variables (elevation, slope, and distance to river)
slope. Finally, we compared the average size values to their aspect (for both N,E,S,W and for N, NW, W, SW, S, SE, E, NE) with ANOVAs/Kruskal-Wallis 
Models. We used 15 m elevation/slope/aspect/distance to river rasters. We observed that for elevation vs. size/shape metrics, Las Matancitas showed significant positive relationships between short canopy axis,
long canopy axis, canopy area, and crown spread and elevation, La Cobriza showed a significant negative relationship between DBH and elevation, San Dionisio showed significant negative relationships between SCA, CA, CS,
and DBH and elevation, and across all populations, there was a significant positive relationship between DBH and elevation. As for slope vs. size/shape metrics, Las Matancitas showed a significant negative relationship between LCA and slope,
San Dionisio showed significant negative relationships between SCA, LCA, CA, and CS and slope, and across all populations, LCA, CA, and CS had significant negative relationships with slope. We observed that distance to river never had a significant relationship with the size and shape of the trees, possibly 
because of how close the trees are to water access already. As for comparing the 8 and 4 categories of 
aspect with the size/shape of trees, in general, we observed the north and northeast facing aspects tending to have larger mean sizes compared to more southern and western facing aspects. 

* [hypothesis_3_GAM_size_and_exposure_CLEANED.R](./analyses/hypothesis_3_GAM_size_and_exposure_CLEANED.R) is a well-commented R script to test hypothesis 3 using **Generalized Additive Models (GAMs)**
as a non-parametric and more flexible method to explore relationships between elevation, slope, and aspect and size characteristics (SCA, LCA, CS, CA, DBH), allowing for non-linear relationships. 

* [hypothesis_3_GAM_water_availability_CLEANED.R](./analyses/hypothesis_3_GAM_water_availability_CLEANED.R) is a well-commented
R script to test hypothesis 5 in which we used **Generalized Additive Models (GAMs)** to find any linear or non-linear relationships
between distance to river, slope, aspect, and elevation (water availability proxies) and the size/shape of the trees. We also explored for 
any interactions. We observed, in general, significant, non-linear relationships between elevation and size/shape metrics for Las Matancias, and significant
non-linear relationships between slopes and size/shape metrics for San Dionisio. La Cobriza tree sizes did not show much significant influences from water availability proxies. 

**Hypothesis 4:** *Quercus brandegeei* persist in a similar, narrow niche breadth of soil characteristics. The 20 known populations have more similar soil characteristics (texture, moisture, Ph, etc.) than 
areas where we know there are no *Q. brandegeei* populations. For this we performed three analyses. The first was a comparison of the mean soil values between our three populations: La Cobriza, 
San Dionisio, and Las Matancitas using ANOVA and non-parametric ANOVA subsitute tests. The second was to compare the 
slopes of size values vs. soil characteristics for each population to the slopes of the same comparison but size values 
were shuffled and randomized (reflected in histogram of slopes). Finally, we compared the mean soil values for our known 
populations to randomly selected populations.
This hypothesis is explored by three scripts:

* [hypothesis_4_soil_between_populations_CLEANED.R](./analyses/hypothesis_4_soil_between_populations_CLEANED.R) is a well-commented R script to 
test hypothesis 4 by exploring whether the three populations (Las Mantacitas, La Cobriza, and San Dionisio) have similar mean soil metric values using **Difference in Means tests**.
We compared the mean soil metrics (soil texture, pH, carbon, water, nitrogen, etc.) between the three populations using parametric ANOVA test (and
the post-hoc Tukey's HSD) or a non-parametric test (Kruskal-Wallis test with the post-hoc Wilcoxon Rank Sum Tests or Welch's ANOVA with post-hoc
Tamhane's test) depending on how well conditions were met. Because of previously observed spatial autocorrelation, we used grid cells to randomly select 
trees from each population for this analysis. We used 250 m resolution soil rasters from Soil Grids to extract the soil metrics. 
When running the non-parametric Kruskal-Wallis Test across all of the population-soil metric combinations, we observed that all of 
the populations showed significant differences in mean soil metrics except for clay/loam field capacity at 100-200 cm. 
These results suggest water availability may explain why the trees are found across these particular locations. 

* [hypothesis_4_soil_characteristics_vs_size_CLEANED.R](./analyses/hypothesis_4_soil_characteristics_vs_size_CLEANED.R) is a well-commented 
R script to test hypothesis 4 by using **Single Linear Regressions (SLRs)** to investigate if there are significant linear relationships between each combination of 
soil metric (soil texture, pH, carbon, water, nitrogen, etc.) and size metrics (Short Canopy Axis, Long Canopy Axis, Canopy Area, Crown Spread, 
and DBH) for each population (Las Mantacitas, La Cobriza, and San Dionisio). To test the validity of our results, we compared the observed SLR slopes
to slopes where the size variables were randomly shuffled across the trees (1000 permutations). We observed that, in general, water availability and soil texture
tended to relate to the size and shapes of the trees most often. 

* [hypothesis_4_soils_between_pops_and_no_pops_CLEANED.R](./analyses/hypothesis_4_soils_between_pops_and_no_pops_CLEANED.R) is a well-commented R script
to test hypothesis 4 by using a **Spatial Null Model Analysis** comparing the mean soil metrics of our 20 known populations to mean soil metrics
if the populations were randomly distributed in a 7000 m buffer zone around the known populations (1000 permutations). We observed that our known populations
had significantly higher clay/loam field capacity (0-5 and 100-200 cm) and sand available water (0-5 cm), suggesting water availability may 
influence the distinct locations of the populations. 

The results are summarized in this google spreadsheet: [General_Test_Results_Wanger_REEF_2025](https://docs.google.com/spreadsheets/d/1BemVj7ev1UcTnCs2zXe8JUpJoGbSO78u0YJfAZuehLs/edit?gid=0#gid=0) 

