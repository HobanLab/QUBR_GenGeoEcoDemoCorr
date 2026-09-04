# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%Looking to see if Q. brandegeei focal populations vary in DBH%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# The purpose of this script is to evaluated whether the size of Quercus brandegeei 
# individuals between the focal populations differs in mean, median, sd, and spread

# The script is broken into sections of
# 1) Loading and processing the packages and processed data for the trees
# 2) Difference in variance tests

# NOTE: Uncomment and run line 36, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 

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
library(raster) #to use crop
library(starsExtra) #to use dist_to_nearest
library(geostatsp) # To successfully use as.im
library(tmaptools)

# loading in the processed tree data 
# NOTE: uncomment and run line 51, sourcing Data_Processing_Script.R, if the line has not yet to be run across any of the scripts/the environment has been cleared 
# source("./analyses/Data_Processing_Script.R")

#ensuring there is a column from latitude and longitude in the populations transformed dataframe because those columns are needed in "hypothesis_1_clumping_CLEANED.R" 

LM_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LM")

LC_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "LC") 

SD_fixed_field_data_processed <- fixed_field_data_processed %>%
  filter(Locality == "SD") 


#### checking for normality of residuals #### 

shapiro_wilk_results <- shapiro.test(fixed_field_data_processed$DBH_ag)
shapiro_wilk_results
hist(fixed_field_data_processed$DBH_ag)

shapiro_wilk_LM <- shapiro.test(LM_fixed_field_data_processed$DBH_ag)
shapiro_wilk_LM
hist(LM_fixed_field_data_processed$DBH_ag)

shapiro_wilk_LC <- shapiro.test(LC_fixed_field_data_processed$DBH_ag)
shapiro_wilk_LC
hist(LC_fixed_field_data_processed$DBH_ag)

shapiro_wilk_SD <- shapiro.test(SD_fixed_field_data_processed$DBH_ag)
shapiro_wilk_SD
hist(SD_fixed_field_data_processed$DBH_ag)


#### Difference in Variance Tests ####

## if residuals were normal, we used Bartlett's Tests 

bartlett.test(DBH_ag ~ Locality, data = fixed_field_data_processed)

# Using a list of vectors
bartlett.test(list(LM, LC, SD))

#add Bonferroni p-value corrections

## if residuals were not normal, we used Fligner-Kileen Test

#Fligner-Killeen, more useful than the other equal variance tests when dealing with non-normal and when outliers present
fligner.test(DBH_ag ~ Locality, data = fixed_field_data_processed)


#post-hoc comparisons
fligner.test(DBH_ag ~ Locality, data = filter(
  fixed_field_data_processed,
  Locality %in% c("LM", "LC")))

fligner.test(DBH_ag ~ Locality, data = filter(
  fixed_field_data_processed,
  Locality %in% c("LM", "SD")))

fligner.test(DBH_ag ~ Locality, data = filter(
  fixed_field_data_processed,
  Locality %in% c("LC", "SD")))

my_p_values <- c(0.00000000005559, 0.2625, 0.0000002502, 0.0000000008461)
p.adjust(my_p_values, method = "fdr")


#### Difference in Means Tests ####

kruskal.test(DBH_ag ~ Locality, data = fixed_field_data_processed)

diff_mean_post_hoc <- pairwise.wilcox.test(fixed_field_data_processed$DBH_ag, fixed_field_data_processed$locality)
diff_mean_post_hoc

#add Bonferroni p-value corrections

diff_mean_post_hoc_corrected <- pairwise.wilcox.test(fixed_field_data_processed$DBH_ag, fixed_field_data_processed$locality, 
                                                     p.adjust.method = "fdr")
diff_mean_post_hoc_corrected

#### Difference in Median Tests ####

library(RVAideMemoire)
#chi-square test

mood.medtest(DBH_ag ~ locality, data  = fixed_field_data_processed, exact = FALSE,
             method = "fdr")

mood.medtest(DBH_ag ~ locality, data  = fixed_field_data_processed, exact = FALSE)

#pairwise post-hoc tests
mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LM", "LC")), 
                                               exact = FALSE)
mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LM", "LC")), 
             exact = FALSE, method = "fdr")

mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LC", "SD")), 
             exact = FALSE)
mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LC", "SD")), 
             exact = FALSE, method = "fdr")

mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LM", "SD")), 
             exact = FALSE)
mood.medtest(DBH_ag ~ locality, data  = subset(fixed_field_data_processed, 
                                               locality %in% c("LM", "SD")), 
             exact = FALSE, method = "fdr")


my_p_values <- c(0.0004118, 0.0004448, 0.007181, 0.7689)
p.adjust(my_p_values, method = "fdr")

#### Figures for Paper ####


plot_colors <- c("#875E3BFF", "#E37F3EFF", "#B3C6ADFF")

DBH_summary <- fixed_field_data_processed %>%
  group_by(locality) %>%
  summarise(mean_DBH = mean(DBH_ag, na.rm = TRUE),
            median_DBH = median(DBH_ag, na.rm = TRUE),
            sd_DBH = sd(DBH_ag, na.rm = TRUE))

#esnruing the order displayed is LM then LC and then SD
fixed_field_data_processed$locality <- factor(fixed_field_data_processed$locality, levels = c("LM", "LC", "SD"))

DBH_summary$locality <- factor(DBH_summary$locality, levels = c("LM", "LC", "SD"))

histograms_populations <- ggplot(data = fixed_field_data_processed, aes(x = DBH_ag))+
  geom_histogram(aes(fill = locality)) +
  geom_vline(data = DBH_summary, aes(xintercept = mean_DBH), color = "red", linewidth = 0.5) +
  geom_vline(data = DBH_summary, aes(xintercept = median_DBH), color = "blue", linewidth = 0.5) +
  geom_text(data = DBH_summary, aes(x = mean_DBH, y = Inf, label = paste0("Mean = ", round(mean_DBH, 3))),
    color = "red", size = 5,  vjust = 1.5, hjust = -0.2) +
  geom_text(data = DBH_summary, aes(x = median_DBH, y = Inf, label = paste0("Median = ", round(median_DBH, 3))),
            color = "blue", size = 5, vjust = 1.5, hjust = 1.25) +
  geom_text(data = DBH_summary, x = 1.5, y = 20, aes(label = paste0("SD = ", round(sd_DBH, 3))),
            color = "black", size = 5, vjust = 1.5, hjust = 1.3) +
  # geom_vline(x = mean(fixed_field_data_processed$DBH_ag, na.rm = T), color = "blue") + 
  # geom_vline(x = median(fixed_field_data_processed$DBH_ag, na.rm = T), color = "red") + 
  facet_wrap(~ locality, nrow = 3, ncol = 1) +
  labs(x = "DBH (m)",
       y = "Frequency",
       fill = "Locality") +
  scale_fill_manual(values = plot_colors) +
  theme_minimal() +
  theme(title=element_text(size=18), 
        axis.text=element_text(size=18),  axis.title.x =element_text(size= 18),
        axis.title.y =element_text(size= 18),
        label =element_text(size= 18, family = "serif"),
        text = element_text(family = "serif"),
        legend.text = element_text(size= 18),
        strip.text = element_text(size = 18, family = "serif", face = "bold"))

histograms_populations



#"skyblue"
#"#5C462B"
# shallower_soil_color <- c("#C4B07B")
# deeper_soil_color <- c("#5C462B")

elevation.plot <- ggplot()+
  geom_boxplot(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED, fill = Locality, color = Locality), alpha = 0.2, , outliers=F)+ #fill = plot_color,
  geom_jitter(data = means_fixed_field_data_processed_terrain_dist, aes(Locality, Elevation..m.FIXED, width = 0.2, fill = Locality, color = Locality), alpha = 0.3)+
  ylab("Elevation (m)")+
  xlab("")+
  scale_fill_manual(values = plot_colors) +
  scale_color_manual(values = plot_colors) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(round(min(means_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED)), 
                                  round(max(means_fixed_field_data_processed_terrain_dist$Elevation..m.FIXED)), 
                                  by = 30)) +
  theme(title=element_text(size=18), 
        axis.text=element_text(size=18),  axis.title.x =element_text(size= 18),
        axis.title.y =element_text(size= 18),
        label =element_text(size= 18, family = "serif"),
        text = element_text(family = "serif"),
        legend.text = element_text(size= 18))
elevation.plot
#extracting the plot information and building a dataframe with the significance letters
elevation.plot.build <- ggplot_build(elevation.plot)
elevation.plot.build.df <- elevation.plot.build$data[[1]] #storing the 0-5 cm boxplots
significance_letters_elevation_df <- data.frame(
  Locality = c("LM", "LC", "SD"), #storing the localities
  y = c(elevation.plot.build.df$ymin[1] - 5,  #adding 5 units to the whisker y-maxs
        elevation.plot.build.df$ymin[2] - 5, 
        elevation.plot.build.df$ymin[3] - 5), 
  label = c("a", "b", "c") #adding the significance labels
)
#adding the significance letters to the plot
elevation.plot <- elevation.plot +
  geom_text(data = significance_letters_elevation_df, 
            aes(x = Locality, y=y, 
                label = label),
            color = plot_colors,
            fontface = "bold",
            size = 7)
elevation.plot




