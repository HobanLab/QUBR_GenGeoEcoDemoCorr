#### Loading libraries and relevant data####

library(tidyverse)
library(moments) # for calculating the moments of each variable
library(sf) # for plotting spatial objects
library(smatr)
library(ggpmisc)
library(PMCMRplus) # for Dunn test
library(geomtextpath) # for PCA graphing

`%notin%` <- Negate(`%in%`) # Make a function that is the opposite of the %in% function

# Read in the raw data from the Field_datasheets_filled_before_KA_check_copy.csv
field_data_raw <- read.csv("./data/Field_datasheets_filled_before_KA_check_copy.csv", na.strings = c("NA", "")) 

#### Process the data ####
field_data_processed <- field_data_raw %>%
  mutate(Canopy1 = as.numeric(Canopy1), 
         Canopy2 = as.numeric(Canopy2)) %>% # Make the canopy valuables numeric bc they read as characters
  mutate(across(starts_with("DBH"), ~ (.x/100)^2)) %>% # Divide all of the DBHs to be in m rather than cm so they are in the same unit as the canopy measurements -> then square all of the DBH's --> then next row of code
  mutate(DBH_ag = sqrt(rowSums(across(starts_with("DBH")), na.rm = T))) %>% # ->  sum all the the new DBH values -> then take the squareroot the sum --> this is the aggregated DBH (DBH_ag) which is how people in the US tend to calculate a single DBH from multi-stem dbh's (which is what 167 of these trees are)
  mutate(multistemmed = case_when(is.na(DBH2) ~ F, 
                                  !is.na(DBH2) ~ T)) %>% # Create a column that is a logical vector (True/False) that describes if the tree has multiple stems (T) or not (F)
  select(!c(DBH1, DBH2, DBH3, DBH4, DBH5, DBH6)) %>% # Now that we have the aggregated DBH, we can remove all other DBH columns
  filter(DBH_ag != 0) %>% # Filter to keep only individuals with DBH measurements (if DBH_ag = 0, then there were only NAs in all DBH columns to start)
  mutate(Canopy_short = case_when(Canopy1 <= Canopy2 ~ Canopy1, 
                                  Canopy2 < Canopy1 ~ Canopy2,), 
         Canopy_long = case_when(Canopy1 >= Canopy2 ~ Canopy1, 
                                  Canopy2 > Canopy1 ~ Canopy2,)) %>% #Create a Canopy_short column and a Canopy_long column from our current, unranked canopy values (Canopy1 and Canopy2); for the Canopy_short column, use the smaller of the 2 values, for the Canopy_long column, use the bigger of the 2 values
  mutate(Crown_spread = (Canopy_short + Canopy_long)/2) %>% #Create a column with crown spread values -> the avg of the two measures of canopy length
  mutate(a = Canopy_long/2,
         b = Canopy_short/2,
         c = sqrt(a^2 - b^2),
         eccentricity = c/a, 
         Canopy_area = pi*a*b) %>% # Create a bunch of columns which are neccesary to calculate eccentricity and canopy area (assuming the canopy is an ellipse); equations for found here: https://www.andrews.edu/~rwright/Precalculus-RLW/Text/07-03.html#:~:text=Eccentricity%20is%20a%20measure%20of,ellipse%20is%20almost%20a%20line., 
  filter(!is.na(Canopy_short) & !is.na(Canopy_long)) %>% # Remove individuals with no canopy measurements
  mutate(Locality = as.factor(Locality)) %>% # Make sure Locality is a factor, not a character
  mutate(W = paste0("-", W)) %>% # Add a negative to the column with the West values so that this can be treated as longitude
  rename(lat = N, 
         long = W) # Rename N and W columns to be lat and long respectively



# Create a dataframe where all lat values are in decimal degrees instead of some being decimal degrees and others being in in degree decimal minutes
fix_deg_dec_min_lat <- field_data_processed %>%
  filter(grepl(" ", lat)) %>% # Filter so new df only has rows with latitude values that contain a space (because this is true of all entries that are in deg_dec_min)
  mutate(lat = measurements::conv_unit(lat, from = 'deg_dec_min', to = 'dec_deg')) %>% # Convert all late values to dec_deg
  select(c(QUBR_ID, lat)) #keep only the QUBR_ID and latitude columns 

# Create a dataframe where all long values are in decimal degrees instead of some being decimal degrees and others being in in degree decimal minutes
fix_deg_dec_min_long <- field_data_processed %>%
  filter(grepl(" ", long)) %>%
  mutate(long = measurements::conv_unit(long, from = 'deg_dec_min', to = 'dec_deg')) %>%
  select(c(QUBR_ID, long))

# Make a new df with all lat and longs in decimal degrees 
fixed_field_data_processed <- field_data_processed %>% 
  left_join(fix_deg_dec_min_long, by=join_by(QUBR_ID)) %>% # Attach the fix_deg_dec_min_long df to this df by the QUBR_ID
  left_join(fix_deg_dec_min_lat, by=join_by(QUBR_ID)) %>% # Attach the fix_deg_dec_min_lat df to this df by the QUBR_ID
  mutate(lat = case_when(is.na(lat.y) ~ lat.x, 
                         !is.na(lat.y) ~ lat.y), 
         long = case_when(is.na(long.y) ~ long.x, 
                          !is.na(long.y) ~ long.y)) %>% # Make a new lat and long column and use the old lat and long (lat.x and long.x) whenever there is no value from the left_joined dfs and use the new values from the left_joined dfs when they exist 
  mutate(lat = as.numeric(lat), 
         long = as.numeric(long)) %>% # Make sure the lat and the long columns are numeric rather than charaters
  select(!c(lat.y, lat.x, long.y, long.x)) # Get rid of the lat and long columns that aren't all in decimal degree


#downloading fixed_field_data_processed as csv to help import into google maps and google earth
write.csv(fixed_field_data_processed, "./analyses/fixed_field_data_processed.csv")

# Make an sf object with all the associated data for each individual
fixed_field_data_sf <- st_as_sf(fixed_field_data_processed, coords = c("long", "lat"), crs= 4326) # Turn the fixed_field_data_processed df into an sf object with the coordinate reference system (crs) of lat/long coords (4326)

#### Describing univariate dist'ns in the 2023 data####

# Create a df which contains the "classical" univariate dist'n stats of each of the important variables
field_data_summarized <- field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity) %>%  # Keep only the columns we are interested in getting summary values of
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE)) # Create columns which summarize the mean, median, variance, and standard deviation of each of the selected columns --> these will be used on the hisogram plots

# Create a df which contains the statistical moments of each of the important variables
field_data_moments_summarized <- field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity) %>%
  summarise(across(everything(), list(skewness= skewness, kurtosis = kurtosis), na.rm = TRUE)) #create columns which summarize the moments (skewness and kurtosis) of each of the selected columns 


#### Basic histograms of each measure of size ####

# DBH_ag histogram
field_data_processed %>% # Use the field_data_processed df for the plot unless otherwise specified
  ggplot() + # Generate the base plot
  geom_histogram(aes(x = DBH_ag, y = ..density..)) + # Plot the aggregated DBH on the x axis and make it a density plot on the y axis
  geom_density(aes(x = DBH_ag)) + # Add a density line over the plot from the DBH_ag values  
  geom_vline(aes(xintercept = field_data_summarized$DBH_ag_mean), color = "red") + # Add a red vertical line that intersects with the x axis at the DBH_ag_mean (from the field_data_summarized df)
  geom_vline(aes(xintercept = field_data_summarized$DBH_ag_median), color = "blue") + # Add a blue vertical line that intersects with the x axis at the DBH_ag_median (from the field_data_summarized df)
  geom_vline(aes(xintercept = field_data_summarized$DBH_ag_mean + field_data_summarized$DBH_ag_sd ), color = "gray") + # Add a grey vertical line that intersects with the x axis at the value of the sum of the DBH_ag_mean  and the DBH_ag_sd (from the field_data_summarized df)
  geom_vline(aes(xintercept = field_data_summarized$DBH_ag_mean - field_data_summarized$DBH_ag_sd), color = "gray") + # Add a grey vertical line that intersects with the x axis at the value of the difference of the DBH_ag_mean  and the DBH_ag_sd (from the field_data_summarized df)
  theme_classic() # Make the plot pretty using the classic theme

# Short canopy length axis (SCA) histogram, the same as the DBH_ag plot but with the SCA values
ggplot(field_data_processed) +
  geom_histogram(aes(x = Canopy_short, y = ..density..)) +
  geom_density(aes(x = Canopy_short)) +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_short_mean), color = "red") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_short_median), color = "blue") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_short_mean + field_data_summarized$Canopy_short_sd ), color = "gray") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_short_mean - field_data_summarized$Canopy_short_sd), color = "gray") +
  xlab("Short Canopy Axis") +
  theme_classic()

#ggsave("./figs/Canopy_short_hist.png", width = 5, height = 5) # Save the most recently run plot as a png with the name Canopy_short_hist in the figs folder

# Long axis of canopy length (LCA) histogram, the same as the DBH_ag plot but with LCA values
ggplot(field_data_processed) +
  geom_histogram(aes(x = Canopy_long, y = ..density..)) +
  geom_density(aes(x = Canopy_long)) +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_long_mean), color = "red") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_long_median), color = "blue") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_long_mean + field_data_summarized$Canopy_long_sd ), color = "gray") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_long_mean - field_data_summarized$Canopy_long_sd), color = "gray") +
  theme_classic()

# Average crown spread (m) histogram, the same as the DBH_ag plot but with crown spread values
ggplot(field_data_processed) +
  geom_histogram(aes(x = Crown_spread, y = ..density..)) +
  geom_density(aes(x = Crown_spread)) +
  geom_vline(aes(xintercept = field_data_summarized$Crown_spread_mean), color = "red") +
  geom_vline(aes(xintercept = field_data_summarized$Crown_spread_median), color = "blue") +
  geom_vline(aes(xintercept = field_data_summarized$Crown_spread_mean + field_data_summarized$Crown_spread_sd ), color = "gray") +
  geom_vline(aes(xintercept = field_data_summarized$Crown_spread_mean - field_data_summarized$Crown_spread_sd), color = "gray") +
  theme_classic()

# Total canopy area in m^2 histogram, the same as the DBH_ag plot but with canopy area values
ggplot(field_data_processed) +
  geom_histogram(aes(x = Canopy_area, y = ..density..)) +
  geom_density(aes(x = Canopy_area)) +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_area_mean), color = "red") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_area_median), color = "blue") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_area_mean + field_data_summarized$Canopy_area_sd ), color = "gray") +
  geom_vline(aes(xintercept = field_data_summarized$Canopy_area_mean - field_data_summarized$Canopy_area_sd), color = "gray") +
  theme_classic()

# Eccentricity (how close to a circle the canopy is with 1 being a perfect circle) histogram, the same as the DBH_ag plot but with eccentricity values
ggplot(field_data_processed) +
  geom_histogram(aes(x = eccentricity, y = ..density..)) +
  geom_density(aes(x = eccentricity)) +
  geom_vline(aes(xintercept = field_data_summarized$eccentricity_mean), color = "red") +
  geom_vline(aes(xintercept = field_data_summarized$eccentricity_median), color = "blue") +
  geom_vline(aes(xintercept = field_data_summarized$eccentricity_mean + field_data_summarized$eccentricity_sd ), color = "gray") +
  geom_vline(aes(xintercept = field_data_summarized$eccentricity_mean - field_data_summarized$eccentricity_sd), color = "gray") +
  theme_classic()


# Q-Q plots; to visually inspect for normality

# DBH qq plot
ggplot(field_data_processed) + # Plot data from the field_data_processed df
  geom_qq(aes(sample = DBH_ag)) + # Make a qq plot using the DBH_ag values
  theme_classic() # Make the plot pretty using the classic theme

# SCA qq plot
ggplot(field_data_processed) +
  geom_qq(aes(sample = Canopy_short)) +
  theme_classic()
#ggsave("./figs/Canopy_short_qqplot.png", width = 3, height = 3)

# Logged SCA qq plot (to check if normality is better when the variable is logged -> it isn't)
ggplot(field_data_processed) +
  geom_qq(aes(sample = log(Canopy_short))) +
  theme_classic()
#ggsave("./figs/Canopy_short_logged_qqplot.png", width = 3, height = 3)

# LCA qq plot
ggplot(field_data_processed) +
  geom_qq(aes(sample = Canopy_long)) +
  theme_classic()

# Crown spread qq plot
ggplot(field_data_processed) +
  geom_qq(aes(sample = Crown_spread)) +
  theme_classic()

# Canopy area qq plot
ggplot(field_data_processed) +
  geom_qq(aes(sample = Canopy_area)) +
  theme_classic()

# Eccentricity qq plot
ggplot(field_data_processed) +
  geom_qq(aes(sample = eccentricity)) +
  theme_classic()


# Official normality testing using the shapiro.test function from the stats library (base R)
shapiro.test(field_data_processed$DBH_ag)
shapiro.test(field_data_processed$Canopy_short)
shapiro.test(field_data_processed$Canopy_long)
shapiro.test(field_data_processed$Crown_spread)
shapiro.test(field_data_processed$Canopy_area)
shapiro.test(field_data_processed$eccentricity)
# Results: none of them are normal (even when logged, see below) which is really not surprising, so we need to proceed with non-parametric tests 

# Official normality testing on the natural log of each variable
shapiro.test(log(field_data_processed$DBH_ag))
shapiro.test(log(field_data_processed$Canopy_short))
shapiro.test(log(field_data_processed$Canopy_long))
shapiro.test(log(field_data_processed$Crown_spread))
shapiro.test(log(field_data_processed$Canopy_area))
shapiro.test(log(field_data_processed$eccentricity))

#### Looking for differences across sites before size standardization ####

# Create a df which contains the "classical" univariate dist'n stats of each of the important variables FOR EACH SITE
field_data_summarized_by_locality <- field_data_processed %>%
  dplyr::select(DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity, Locality) %>%
  group_by(Locality) %>% # This is the line that modifies summarize such that each stat is found for each of the groups (which here is site aka locality)
  summarise(across(everything(), list(mean = mean, median = median, var = var, sd = sd), na.rm=TRUE))


# Testing for differences in means/medians/"locations" of the groups (aka localities) 
## Since my data aren't normal --> instead of ANOVA, use a Kruskal Wallace test
kruskal.test(field_data_processed$DBH_ag, field_data_processed$Locality) # This command is not tidy compatible I think so --> dependent variable, independent variable
kruskal.test(field_data_processed$Canopy_short, field_data_processed$Locality)
kruskal.test(field_data_processed$Canopy_long, field_data_processed$Locality)
kruskal.test(field_data_processed$Crown_spread, field_data_processed$Locality)
kruskal.test(field_data_processed$Canopy_area, field_data_processed$Locality)
kruskal.test(field_data_processed$eccentricity, field_data_processed$Locality)
# Results: Differences by locality in DBH and marginal differences in short canopy axis


# Post-hoc testing 
## For kruskal.test post-hoc testing is the dunn test

kwAllPairsDunnTest(DBH_ag ~ Locality, data = field_data_processed, method="bh") # Dependent ~ Independent variable 
# Results: LM is sig different from both SD and LC

kwAllPairsDunnTest(Canopy_short ~ Locality, data = field_data_processed, method="bh") 
# Results: LM is sig different from both SD and LC

# Visualization of results via grouped boxplots

# DBH
field_data_processed %>% # Use the df that isn't grouped and summarized field_data_processed for data viz
  ggplot() +
  geom_boxplot(aes(y = DBH_ag, x = Locality, color = Locality)) +
  ylab("DBH (m)") +
  theme_classic()
#ggsave("./Figs/DBH_boxplots_by_locality.png", width = 4, height = 4)

#SCA
field_data_processed %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_short, x = Locality, color = Locality)) +
  ylab("Short Canopy Axis (m)") +
  theme_classic()
#ggsave("./Figs/Short_canopy_boxplots_by_locality.png", width = 4, height = 4)


# Comparing variance/"spread" among groups
## Since my data aren't normal --> Fligner-Killeen test found here [https://www.geeksforgeeks.org/homogeneity-of-variance-test-in-r-programming/]
fligner.test(Canopy_short ~ Locality, data = field_data_processed) #sig
fligner.test(Canopy_long ~ Locality, data = field_data_processed) #sig
fligner.test(Crown_spread ~ Locality, data = field_data_processed) #sig
fligner.test(Canopy_area ~ Locality, data = field_data_processed) #marginally sig
fligner.test(eccentricity ~ Locality, data = field_data_processed) #non sig
fligner.test(DBH_ag ~ Locality, data = field_data_processed) #sig
# Results: all are but one are sig 

# Post-hoc testing
## 3 pairwise Filnger-Killeen tests with a Bonferroni correction to adjust alpha value
### Bonferroni correction alpha value for 3 comparisons =  0.0167 (to be 95% confident)

# SCA
fligner.test(Canopy_short ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "LC"))) # Perform the fligner.test and look for diffs in just LM and LC --> #not sig!!!
fligner.test(Canopy_short ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "SD"))) # Perform the fligner.test and look for diffs in just LM and SD --> #not sig!!!
fligner.test(Canopy_short ~ Locality, data = filter(field_data_processed, Locality %in% c("SD", "LC"))) # Perform the fligner.test and look for diffs in just SD and LC --> #sig!!!
# Results: SD and LC have sig different variances in SCA but neither are sig different from LM

# LCA
fligner.test(Canopy_long ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_long ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Canopy_long ~ Locality, data = filter(field_data_processed, Locality %in% c("SD", "LC"))) #sig!!!
# Results: SD has a sig different variance than both LM and LC which are not sig different from one another

# Crown Spread
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed, Locality %in% c("SD", "LC"))) #sig!!!
# Results: SD has a sig different variance than both LM and LC which are not sig different from one another


# Canopy Area
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "SD"))) #not sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed, Locality %in% c("SD", "LC"))) #not sig!!!
# Results: No sig differences in variance across sites after Bonferonni correction

# DBH
fligner.test(DBH_ag ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(DBH_ag ~ Locality, data = filter(field_data_processed, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(DBH_ag ~ Locality, data = filter(field_data_processed, Locality %in% c("SD", "LC"))) #sig!!!
# Results: SD has a sig different variance than both LM and LC which are not sig different from one another


#Visualizing distributions for sig diffs in variance

# Making the plotting easier by modifying the df to have columns with summary stats
field_data_for_plots <- field_data_processed %>%
  dplyr::select(Locality, DBH_ag, Canopy_short, Canopy_long, Crown_spread, Canopy_area, eccentricity) %>% # Keep only the relevant columns
  group_by(Locality) %>% # Group by site 
  mutate(Canopy_long_mean = mean(Canopy_long, na.rm = T), 
         Crown_spread_mean = mean(Crown_spread, na.rm = T),
         DBH_ag_mean = mean(DBH_ag, na.rm = T), 
         Canopy_long_median = median(Canopy_long, na.rm = T), 
         Crown_spread_median = median(Crown_spread, na.rm = T),
         DBH_ag_median = median(DBH_ag, na.rm = T), 
         Canopy_long_sd = sd(Canopy_long, na.rm = T), 
         Crown_spread_sd = sd(Crown_spread, na.rm = T),
         DBH_ag_sd = sd(DBH_ag, na.rm = T), 
         Canopy_short_mean = mean(Canopy_short, na.rm = T), 
         Canopy_short_median = median(Canopy_short, na.rm = T), 
         Canopy_short_sd = sd(Canopy_short, na.rm = T)) # Instead of using summarize, just using mutate here so that the original data doesn't get condensed to just the summary stats because I need both the summarized and the unsummarized data for my plots below


# SCA by site
ggplot(field_data_for_plots) +
  geom_histogram(aes(x = Canopy_short, y = ..density..)) +
  geom_density(aes(x = Canopy_short)) +
  facet_wrap(~Locality, nrow =3) + # Makes it so a plot will be made for each unique locality, nrow = 3 makes it so these plots stack on top of one another 
  geom_vline(aes(xintercept = Canopy_short_mean, group = Locality), colour = 'red') + # Group = ensures that this vline command is performed uniquely for each group
  geom_vline(aes(xintercept = Canopy_short_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = Canopy_short_mean + Canopy_short_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = Canopy_short_mean - Canopy_short_sd, group = Locality), colour = 'gray') +
  xlab("Short Canopy Axis (m)") +
  theme_classic()
#ggsave("./Figs/Canopy_short_dists_by_locality.png", width = 5, height = 5)


# LCA by site
ggplot(field_data_for_plots) +
  geom_histogram(aes(x = Canopy_long, y = ..density..)) +
  geom_density(aes(x = Canopy_long)) +
  facet_wrap(~Locality, nrow =3) +
  geom_vline(aes(xintercept = Canopy_long_mean, group = Locality), colour = 'red') +
  geom_vline(aes(xintercept = Canopy_long_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = Canopy_long_mean + Canopy_long_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = Canopy_long_mean - Canopy_long_sd, group = Locality), colour = 'gray') +
  theme_classic()

# Average crown spread (m) by site
ggplot(field_data_for_plots) +
  geom_histogram(aes(x = Crown_spread, y = ..density..)) +
  geom_density(aes(x = Crown_spread)) +
  facet_wrap(~Locality, nrow =3) +
  geom_vline(aes(xintercept = Crown_spread_mean, group = Locality), colour = 'red') +
  geom_vline(aes(xintercept = Crown_spread_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = Crown_spread_mean + Crown_spread_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = Crown_spread_mean - Crown_spread_sd, group = Locality), colour = 'gray') +
  theme_classic()

#DBH by site
ggplot(field_data_for_plots) +
  geom_histogram(aes(x = DBH_ag, y = ..density..)) +
  geom_density(aes(x = DBH_ag)) +
  facet_wrap(~Locality, nrow =3) +
  geom_vline(aes(xintercept = DBH_ag_mean, group = Locality), colour = 'red') +
  geom_vline(aes(xintercept = DBH_ag_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = DBH_ag_mean + DBH_ag_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = DBH_ag_mean - DBH_ag_sd, group = Locality), colour = 'gray') +
  xlab("DBH (m)") +
  theme_classic()
#ggsave("./Figs/DBH_dists_by_locality.png", width = 5, height = 5)

#### Size standardization assuming DBH is best measurement of size ####
# Make a new df with all measurements standardized by DBH
field_data_processed_standardized_by_DBH <- field_data_processed %>%
  dplyr::select(Locality, DBH_ag, Canopy_short, Canopy_long) %>% # Keep relevant columns
  mutate(Locality = as.factor(Locality),
         Canopy_short_std = Canopy_short/DBH_ag, 
         Canopy_long_std = Canopy_long/DBH_ag) %>% # Divide SCA and LCA (the other measurements in the dataset) by DBH
  mutate(Crown_spread = (Canopy_short_std + Canopy_long_std)/2,
         a = Canopy_long_std/2,
         b = Canopy_short_std/2,
         c = sqrt(a^2 - b^2),
         eccentricity = c/a, 
         Canopy_area = pi*a*b) # Recalculate the other measurements based on the new LCA and SCA values

# Now I will do the exact same procedure as I did before SS: normality testing -> test for diffs in location --> location post hoc tests --> test for diffs in spread --> spread post hoc tests --> data viz of sig diffs

# Normality testing
shapiro.test(field_data_processed_standardized_by_DBH$Canopy_short_std)
shapiro.test(field_data_processed_standardized_by_DBH$Canopy_long_std)
shapiro.test(field_data_processed_standardized_by_DBH$Crown_spread)
shapiro.test(field_data_processed_standardized_by_DBH$Canopy_area)
shapiro.test(field_data_processed_standardized_by_DBH$eccentricity)
# Results: none are normal --> try logging them (except eccentricity since that is from 0-1)

# Logged normality testing
shapiro.test(log(field_data_processed_standardized_by_DBH$Canopy_short_std))
shapiro.test(log(field_data_processed_standardized_by_DBH$Canopy_long_std))
shapiro.test(log(field_data_processed_standardized_by_DBH$Crown_spread))
shapiro.test(log(field_data_processed_standardized_by_DBH$Canopy_area))
shapiro.test(log(field_data_processed_standardized_by_DBH$eccentricity))
# Results: none of them are normal (even when logged) which is really not surprising, so we need to proceed with non-parametric tests 


# Testing for differences in means/medians/"locations" of the groups (aka localities) 
## Since my data aren't normal --> instead of ANOVA, use a Kruskal Wallace test
kruskal.test(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
kruskal.test(Canopy_long_std ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
kruskal.test(Crown_spread ~ Locality, data = field_data_processed_standardized_by_DBH) #sig
kruskal.test(Canopy_area ~ Locality, data = field_data_processed_standardized_by_DBH) #sig
kruskal.test(eccentricity~ Locality, data = field_data_processed_standardized_by_DBH) # not sig
# Results: there are diffs in shape (after accounting for diffs in DBH as a size proxy)!!!


# Post-hoc testing 
## For kruskal.test post-hoc testing is the dunn test
kwAllPairsDunnTest(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_DBH, method="bh") #LM is sig different from both SD and LC
kwAllPairsDunnTest(Canopy_long_std ~ Locality, data = field_data_processed_standardized_by_DBH, method="bh") #LM is sig different from both SD and LC
kwAllPairsDunnTest(Crown_spread ~ Locality, data = field_data_processed_standardized_by_DBH, method="bh") #LM is sig different from both SD and LC
kwAllPairsDunnTest(Canopy_area ~ Locality, data = field_data_processed_standardized_by_DBH, method="bh")  #LM is sig different from both SD and LC


# Comparing variance/"spread" among groups
## Since my data aren't normal --> Fligner-Killeen test found here [https://www.geeksforgeeks.org/homogeneity-of-variance-test-in-r-programming/]
fligner.test(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
fligner.test(Canopy_long_std ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
fligner.test(Crown_spread ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
fligner.test(Canopy_area ~ Locality, data = field_data_processed_standardized_by_DBH) #sig!
fligner.test(eccentricity ~ Locality, data = field_data_processed_standardized_by_DBH) #non sig

# Post-hoc testing
## 3 pairwise Filnger-Killeen tests with a Bonferroni correction to adjust alpha value
### Bonferroni correction alpha value for 3 comparisons =  0.0167 (to be 95% confident)

# SCA
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "SD"))) # sig!!!
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("SD", "LC"))) #sig!!!
# Results: sig diffs between SD and both other sites which are not diff from each other

# LCA
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("SD", "LC"))) #sig!!!
# Results: sig diffs between SD and both other sites which are not diff from each other

# Crown spread
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("SD", "LC"))) #sig!!!
# Results: sig diffs between SD and both other sites which are not diff from each other

#all site tests for canopy area
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "LC"))) #sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_DBH, Locality %in% c("SD", "LC"))) #sig!!!
# Results: each site is sig diff from each other site

#Group boxplots
field_data_processed_standardized_by_DBH %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_short_std, x = Locality, color = Locality)) +
  ylab("Standardized Short Canopy Axis (m)") +
  theme_classic()
#ggsave("./Figs/Canopy_short_std_boxplot_by_locality.png", width = 4, height = 4)

field_data_processed_standardized_by_DBH %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_long_std, x = Locality, color = Locality)) +
  theme_classic()

field_data_processed_standardized_by_DBH %>% 
  ggplot() +
  geom_boxplot(aes(y = Crown_spread, x = Locality, color = Locality)) +
  theme_classic()

field_data_processed_standardized_by_DBH %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_area, x = Locality, color = Locality)) +
  theme_classic()

# Making the plotting easier by modifying the df to have columsn with summary stats
field_data_DBH_SS_for_plots <- field_data_processed_standardized_by_DBH %>%
  dplyr::select(Locality, DBH_ag, Canopy_short_std, Canopy_long_std, Crown_spread, Canopy_area, eccentricity) %>%
  group_by(Locality) %>%
  mutate(Canopy_long_std_mean = mean(Canopy_long_std, na.rm = T), 
         Crown_spread_mean = mean(Crown_spread, na.rm = T),
         DBH_ag_mean = mean(DBH_ag, na.rm = T), 
         Canopy_long_std_median = median(Canopy_long_std, na.rm = T), 
         Crown_spread_median = median(Crown_spread, na.rm = T),
         DBH_ag_median = median(DBH_ag, na.rm = T), 
         Canopy_long_std_sd = sd(Canopy_long_std, na.rm = T), 
         Crown_spread_sd = sd(Crown_spread, na.rm = T),
         DBH_ag_sd = sd(DBH_ag, na.rm = T), 
         Canopy_short_std_mean = mean(Canopy_short_std, na.rm = T), 
         Canopy_short_std_median = median(Canopy_short_std, na.rm = T), 
         Canopy_short_std_sd = sd(Canopy_short_std, na.rm = T))


# SCA 
ggplot(field_data_DBH_SS_for_plots) +
  geom_histogram(aes(x = Canopy_short_std, y = ..density..), binwidth = 2) +
  geom_density(aes(x = Canopy_short_std)) +
  facet_wrap(~Locality, nrow =3) +
  geom_vline(aes(xintercept = Canopy_short_std_mean, group = Locality), colour = 'red') +
  geom_vline(aes(xintercept = Canopy_short_std_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = Canopy_short_std_mean + Canopy_short_std_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = Canopy_short_std_mean - Canopy_short_std_sd, group = Locality), colour = 'gray') +
  xlab("Short Canopy Axis SS by DBH") +
  theme_classic()
#("./Figs/Canopy_short_std_by_DBH_dists_by_locality.png", width = 5, height = 5)

#### Size standardization using all traditional morphometric measurements ####
# Make a new df with the 3 measurements all size standardized by each other 
field_data_processed_standardized_by_centroid <- field_data_processed %>%
  dplyr::select(Locality, DBH_ag, Canopy_short, Canopy_long) %>%
  filter(!is.na(Canopy_short) & !is.na(Canopy_long)) %>%
  mutate(Locality = as.factor(Locality),
         Centroid_size = sqrt(sum(DBH_ag^2, Canopy_short^2, Canopy_long^2))) %>%
  mutate(DBH_std = DBH_ag/Centroid_size,
         Canopy_short_std = Canopy_short/Centroid_size, 
         Canopy_long_std = Canopy_long/Centroid_size, 
         sanity_check =  sqrt(sum(DBH_std^2, Canopy_short_std^2, Canopy_long_std^2))) %>% 
  mutate(Crown_spread = (Canopy_short_std + Canopy_long_std)/2,
         a = Canopy_long_std/2,
         b = Canopy_short_std/2,
         c = sqrt(a^2 - b^2),
         eccentricity = c/a, 
         Canopy_area = pi*a*b)



# Normality testing
shapiro.test(field_data_processed_standardized_by_centroid$Canopy_short_std)
shapiro.test(field_data_processed_standardized_by_centroid$Canopy_long_std)
shapiro.test(field_data_processed_standardized_by_centroid$Crown_spread)
shapiro.test(field_data_processed_standardized_by_centroid$Canopy_area)
shapiro.test(field_data_processed_standardized_by_centroid$eccentricity)
shapiro.test(field_data_processed_standardized_by_centroid$DBH_std)
# Results: all are non-normal

# Logges normality testing to see if it helps
shapiro.test(log(field_data_processed_standardized_by_centroid$Canopy_short_std))
shapiro.test(log(field_data_processed_standardized_by_centroid$Canopy_long_std))
shapiro.test(log(field_data_processed_standardized_by_centroid$Crown_spread))
shapiro.test(log(field_data_processed_standardized_by_centroid$Canopy_area))
shapiro.test(log(field_data_processed_standardized_by_centroid$DBH_std))
# Results: none of them are normal (even when logged) which is really not surprising, so we need to proceed with non-parametric tests 


# Testing for differences in means/medians/"locations" of the groups (aka localities) 
## Since my data aren't normal --> instead of ANOVA, use a Kruskal Wallace test
kruskal.test(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_centroid) #marginally sig
kruskal.test(Canopy_long_std ~ Locality, data = field_data_processed_standardized_by_centroid) #not sig
kruskal.test(Crown_spread ~ Locality, data = field_data_processed_standardized_by_centroid) #not sig
kruskal.test(Canopy_area ~ Locality, data = field_data_processed_standardized_by_centroid) #not sig
kruskal.test(eccentricity~ Locality, data = field_data_processed_standardized_by_centroid) # not sig
kruskal.test(DBH_std ~ Locality, data = field_data_processed_standardized_by_centroid) #very sig
# Results: there are less diffs in shape than when we assume that DBH the best measure of size

# Comparing variance/"spread" among groups
## Since my data aren't normal --> Fligner-Killeen test found here [https://www.geeksforgeeks.org/homogeneity-of-variance-test-in-r-programming/]
fligner.test(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_centroid) #sig!
fligner.test(Canopy_long_std ~ Locality, data = field_data_processed_standardized_by_centroid) #sig!
fligner.test(Crown_spread ~ Locality, data = field_data_processed_standardized_by_centroid) #sig!
fligner.test(Canopy_area ~ Locality, data = field_data_processed_standardized_by_centroid) #marginally sig!
fligner.test(eccentricity ~ Locality, data = field_data_processed_standardized_by_centroid) #non sig
fligner.test(DBH_std ~ Locality, data = field_data_processed_standardized_by_centroid) #sig!

#Testing which means differ from one another using a Dunn test (from here: https://rcompanion.org/handbook/F_08.html#:~:text=Probably%20the%20most%20popular%20post,making%20a%20type%2DI%20error.)

kwAllPairsDunnTest(Canopy_short_std ~ Locality, data = field_data_processed_standardized_by_centroid, method="bh") #none are sig with corrected p value
kwAllPairsDunnTest(DBH_std ~ Locality, data = field_data_processed_standardized_by_centroid, method="bh") #LM is sig different from both SD and LC


# Post-hoc testing
## 3 pairwise Filnger-Killeen tests with a Bonferroni correction to adjust alpha value
### Bonferroni correction alpha value for 3 comparisons =  0.0167 (to be 95% confident)

# SCA
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "SD"))) #not sig!!!
fligner.test(Canopy_short_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("SD", "LC"))) #sig!!!


# LCA
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Canopy_long_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("SD", "LC"))) #sig!!!


# Crown spread
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(Crown_spread ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("SD", "LC"))) #sig!!!


# Canopy area
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "SD"))) #not sig!!!
fligner.test(Canopy_area ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("SD", "LC"))) #not sig!!!


# DBH
fligner.test(DBH_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "LC"))) #not sig!!!
fligner.test(DBH_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("LM", "SD"))) #sig!!!
fligner.test(DBH_std ~ Locality, data = filter(field_data_processed_standardized_by_centroid, Locality %in% c("SD", "LC"))) #sig!!!


# Data viz

# Group boxplots
field_data_processed_standardized_by_centroid %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_short_std, x = Locality, color = Locality)) +
  ylab("Standardized Short Canopy Axis") +
  theme_classic()
#ggsave("./Figs/Canopy_short_std_centroid_boxplot_by_locality.png", width = 4, height = 4)

field_data_processed_standardized_by_centroid %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_long_std, x = Locality, color = Locality)) +
  theme_classic()

field_data_processed_standardized_by_centroid %>% 
  ggplot() +
  geom_boxplot(aes(y = Crown_spread, x = Locality, color = Locality)) +
  theme_classic()

field_data_processed_standardized_by_centroid %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_area, x = Locality, color = Locality)) +
  theme_classic()

field_data_processed_standardized_by_centroid %>% 
  ggplot() +
  geom_boxplot(aes(y = DBH_std, x = Locality, color = Locality)) +
  ylab("Standardized DBH") +
  theme_classic()
#ggsave("./Figs/DBH_std_centroid_boxplot_by_locality.png", width = 4, height = 4)


# Making the plotting easier by modifying the df to have columsn with summary stats
field_data_centroid_SS_for_plots <- field_data_processed_standardized_by_centroid %>%
  dplyr::select(Locality, DBH_std, Canopy_short_std, Canopy_long_std, Crown_spread, Canopy_area, eccentricity) %>%
  group_by(Locality) %>%
  mutate(Canopy_long_std_mean = mean(Canopy_long_std, na.rm = T), 
         Crown_spread_mean = mean(Crown_spread, na.rm = T),
         DBH_std_mean = mean(DBH_std, na.rm = T), 
         Canopy_long_std_median = median(Canopy_long_std, na.rm = T), 
         Crown_spread_median = median(Crown_spread, na.rm = T),
         DBH_std_median = median(DBH_std, na.rm = T), 
         Canopy_long_std_sd = sd(Canopy_long_std, na.rm = T), 
         Crown_spread_sd = sd(Crown_spread, na.rm = T),
         DBH_std_sd = sd(DBH_std, na.rm = T), 
         Canopy_short_std_mean = mean(Canopy_short_std, na.rm = T), 
         Canopy_short_std_median = median(Canopy_short_std, na.rm = T), 
         Canopy_short_std_sd = sd(Canopy_short_std, na.rm = T))


# SCA
ggplot(field_data_centroid_SS_for_plots) +
  geom_histogram(aes(x = DBH_std, y = ..density..)) +
  geom_density(aes(x = DBH_std)) +
  facet_wrap(~Locality, nrow =3) +
  geom_vline(aes(xintercept = DBH_std_mean, group = Locality), colour = 'red') +
  geom_vline(aes(xintercept = DBH_std_median, group = Locality), colour = 'blue') +
  geom_vline(aes(xintercept = DBH_std_mean + DBH_std_sd, group = Locality), colour = 'gray') +
  geom_vline(aes(xintercept = DBH_std_mean - DBH_std_sd, group = Locality), colour = 'gray') +
  xlab("DBH SS by unit size") +
  theme_classic()
#ggsave("./Figs/DBH_std_by_centroid_dists_by_locality.png", width = 5, height = 5)

#### All sections below this are not in the paper but are interesting ####
#### Checking to see if there is a difference in sizes between multistemmed inds and single stemmed individuals ####
#Comparing means between 2 groups of non-normal data with wilcox tests
wilcox.test(DBH_ag ~ multistemmed, data = field_data_processed)
wilcox.test(Canopy_short ~ multistemmed, data = field_data_processed)
wilcox.test(Canopy_long ~ multistemmed, data = field_data_processed)
wilcox.test(Canopy_area ~ multistemmed, data = field_data_processed) 
#all above are sig
wilcox.test(eccentricity ~ multistemmed, data = field_data_processed) #not sig

#Group boxplots show that multistemmed trees are bigger
field_data_processed %>% 
  ggplot() +
  geom_boxplot(aes(y = DBH_ag, x = multistemmed, color = multistemmed)) +
  theme_classic() 

field_data_processed %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_short, x = multistemmed, color = multistemmed)) +
  theme_classic()

field_data_processed %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_long, x = multistemmed, color = multistemmed)) +
  theme_classic()

field_data_processed %>% 
  ggplot() +
  geom_boxplot(aes(y = Canopy_area, x = multistemmed, color = multistemmed)) +
  theme_classic()

#are there more (bigger) multistemmed trees at any particular site? --> no!
field_data_processed %>%
  group_by(Locality, multistemmed) %>%
  summarise(n=n()) %>%
  mutate(rel.freq = n/sum(n))

#### PCA ####
field_data_pca <- field_data_processed %>% 
  select(c(DBH_ag, Canopy_short, Canopy_long)) %>%
  prcomp(., center=TRUE, tol=NULL, scale=TRUE)

axis1_var_exp <- summary(field_data_pca)$importance[2,1] #calculate percentage of variance explained on PC1
axis2_var_exp <- summary(field_data_pca)$importance[2,2] #calculate percentage of variance explained on PC2

#plot the PCA and where each variable loads

as_tibble(field_data_pca$rotation) %>%
  mutate(var = rownames(field_data_pca$rotation)) %>%
  ggplot(aes(x = 0, y = 0)) +
    geom_point() +
    geom_textsegment(aes(label=var, xend = PC1, yend = PC2), arrow = arrow()) +
    geom_hline(yintercept = 0, color = "grey") +
    geom_vline(xintercept = 0, color = "grey") +
    xlab(paste0("PC1 (", round(axis1_var_exp*100, digits = 2), "%)", sep = "")) +
    ylab(paste0("PC2 (", round(axis2_var_exp*100, digits = 2), "%)", sep= "")) +
    theme_classic()

field_data_processed %>%
  cbind(field_data_pca$x) %>%
  ggplot(aes(x = PC1, y = PC2, color = Locality)) +
  #geom_textsegment(aes(label=var, xend = PC1, yend = PC2), arrow = arrow()) +
  geom_hline(yintercept = 0, color = "grey") +
  geom_vline(xintercept = 0, color = "grey") +
  xlab(paste0("PC1 (", round(axis1_var_exp*100, digits = 2), "%)", sep = "")) +
  ylab(paste0("PC2 (", round(axis2_var_exp*100, digits = 2), "%)", sep= "")) +
  geom_point() +
  #geom_text(aes(label = QUBR_ID))+
  theme_classic()
