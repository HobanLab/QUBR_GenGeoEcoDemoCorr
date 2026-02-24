###Loading relevant libraries + functions + data####
library(tidyverse)
library(sf)

`%notin%` <- Negate(`%in%`) # Make a function that is the opposite of the %in% function

# Read in the raw data from 2023 sampling trip
field_data_raw <- read.csv("./data/Field_datasheets_filled_before_KA_check_copy.csv", na.strings = c("NA", "")) 

# Read in the csv with Daniels pt data
daniel_data_raw <- read.csv("./data/2021Data_populations_metal_labels.csv", na.strings = "") 

# Read in the raw data from 2025 sampling trip
raw_2025_data <- read_csv("./data/Nov_2025_GPS_data/Revised_QUBR_locations.csv")

# Read in the the bad data taken in the wrong layer (now called mother_points_ruined) during the 2025 sampling trip
LC_weird_data <- read_csv("./data/Nov_2025_GPS_data/LC_data.csv") %>%
  select(!c(lat, long_))

###Adding Daniel's pts into Ash's 2023 data####
daniel_data_in_ash_data <- daniel_data_raw %>%
  mutate(W = paste0("-", W)) %>%
  rename(lat = N, 
         long = W,
         DBH = DAP.Total..m., 
         tree_num = Núm., 
         Daniel_tag = Num..Metal.label, 
         region = Region, 
         locality = Localidad, 
         altitude = Altitud..msnm., 
         fruited_2020 = Con.Bellota..2020., 
         collected_2020 = Colecta..2020., 
         height = Altura..m., 
         DRC = D..Basal..m., 
         Canopy1 = Cob.1..m., 
         Canopy2 = Cob.2..m., 
         notes = Observaciones) %>%
  mutate(across(c(lat, long), ~ gsub('° ', ' ', ., fixed = T)), 
         across(c(lat, long), ~ gsub('. ', '.', ., , fixed = T)), 
         across(c(lat, long), ~ measurements::conv_unit(., from = 'deg_dec_min', to = 'dec_deg')), 
         across(c(lat, long), ~ as.numeric(.))) %>%
  mutate(across(starts_with("DBH"), ~ (.x/pi*100))) %>% #converting their circumference measurements in m to DBH in cm
  mutate(DBH1 = case_when(is.na(DBH1) ~ DBH, 
                          !is.na(DBH1) ~ DBH1)) %>% #making it so all DBH's are in the DBH1 column
  filter(Daniel_tag %in% field_data_raw$Daniel_tag) %>% #filter to only keep rows (based on the number on the silver tags placed by Daniel) which are in Ash's field data
  select(Daniel_tag, DBH1, DBH2, DBH3, DBH4, Canopy1, Canopy2, lat, long)

#need to figure out best way to pull Daniel data into ash data.... thought about left_join but resulted in a shit ton of duplicated columns that I didn't need --> perhaps best to reduce both datasets down to the bare necessities, rbind and then left join that back onto the big ash df?
field_data_processing <- field_data_raw %>%
  left_join(., daniel_data_in_ash_data, by = "Daniel_tag") %>%
  mutate(DBH1.x = ifelse(!is.na(Daniel_tag), DBH1.y, DBH1.x), 
         DBH2.x = ifelse(!is.na(Daniel_tag), DBH2.y, DBH2.x), 
         DBH3.x = ifelse(!is.na(Daniel_tag), DBH3.y, DBH3.x), 
         DBH4.x = ifelse(!is.na(Daniel_tag), DBH4.y, DBH4.x), 
         Canopy1.x = ifelse(!is.na(Daniel_tag), Canopy1.y, Canopy1.x), 
         Canopy2.x = ifelse(!is.na(Daniel_tag), Canopy2.y, Canopy2.x)) %>%
  rename(Daniel_lat = lat, 
         Daniel_long = long) %>%
  select(!ends_with(".y")) %>% #get rid of duplicated cols from Daniel's data
  rename_with( ~ str_remove(., ".x")) #get rid of the .x remaning from the original cols


field_data_processed <- field_data_processing %>%
  mutate(Canopy1 = as.numeric(Canopy1), 
         Canopy2 = as.numeric(Canopy2)) %>%
  mutate(across(starts_with("DBH"), ~ (.x/100)^2)) %>% #Squaring all of the DBH's and then summing them and then sqrting the sum is how people in the US tend to use multi-stem dbh's (which is what 167 of my trees are)
  mutate(DBH_ag = sqrt(rowSums(across(starts_with("DBH")), na.rm = T))) %>%
  mutate(multistemmed = case_when(is.na(DBH2) ~ F, 
                                  !is.na(DBH2) ~ T)) %>% # adding a column that is a logical vector that describes if the tree has multiple stems or not
  select(!c(DBH1, DBH2, DBH3, DBH4, DBH5, DBH6)) %>%
  #filter(DBH_ag != 0) %>% #remove individuals with no DBH measurements
  
  mutate(Canopy_short = case_when(Canopy1 <= Canopy2 ~ Canopy1, 
                                  Canopy2 < Canopy1 ~ Canopy2,), 
         Canopy_long = case_when(Canopy1 >= Canopy2 ~ Canopy1, 
                                 Canopy2 > Canopy1 ~ Canopy2,),
         Crown_spread = (Canopy_short + Canopy_long)/2,
         a = Canopy_long/2,
         b = Canopy_short/2,
         c = sqrt(a^2 - b^2),
         eccentricity = c/a, 
         Canopy_area = pi*a*b) %>% #eccentricity is a measure of how circular the canopy is, equations for this found here: https://www.andrews.edu/~rwright/Precalculus-RLW/Text/07-03.html#:~:text=Eccentricity%20is%20a%20measure%20of,ellipse%20is%20almost%20a%20line., crown spread is the avg of our two measures of canopy length
  filter(QUBR_ID != "LM_156") %>% #I used to remove individuals with no canopy measurements so this *used to* REMOVE POINTS WHICH DANIEL MEASURED AND WE DID NOT REMEASURE (probs like 20 inds) which I didn't want and was why I made this script --> then it was removing trees that simply had canopy axes values either of 0/NA or with a character in them which I also didn't want --> after playing around, it was easiest to just remove this single QUBR_ID that has no lat/long
  mutate(Locality = as.factor(Locality)) %>% # make sure Locality is a factor
  mutate(W = paste0("-", W)) %>%
  rename(lat = N, 
         long = W) #fixing lat long names + data

#checking that only the trees I want to drop are dropped when creating the processed dataset
filter(field_data_processing, QUBR_ID %notin% field_data_processed$QUBR_ID) %>%
  select(c(QUBR_ID, Canopy1, Canopy2, N, W, Daniel_lat, Daniel_long)) 
  #indeed! only LM_156 is dropped

# Fixing lat/long points to be in dec deg for mapping#
fix_deg_dec_min_lat <- field_data_processed %>%
  filter(grepl(" ", lat)) %>%
  mutate(lat = measurements::conv_unit(lat, from = 'deg_dec_min', to = 'dec_deg')) %>%
  select(c(QUBR_ID, lat))

fix_deg_dec_min_long <- field_data_processed %>%
  filter(grepl(" ", long)) %>%
  mutate(long = measurements::conv_unit(long, from = 'deg_dec_min', to = 'dec_deg')) %>%
  select(c(QUBR_ID, long))


fixed_field_data_processed <- field_data_processed %>% 
  left_join(fix_deg_dec_min_long, by=join_by(QUBR_ID)) %>%
  left_join(fix_deg_dec_min_lat, by=join_by(QUBR_ID)) %>%
  mutate(lat = case_when(is.na(lat.y) ~ lat.x, 
                         !is.na(lat.y) ~ lat.y), 
         long = case_when(is.na(long.y) ~ long.x, 
                          !is.na(long.y) ~ long.y)) %>%
  mutate(lat = as.numeric(lat), 
         long = as.numeric(long)) %>%
  select(!c(lat.y, lat.x, long.y, long.x))
#all of the lat and long values are now in dec deg!!

fixed_field_data_sf <- st_as_sf(fixed_field_data_processed, coords = c("long", "lat"), crs= 4326) #making a easy to plot on a map dataset even though I won't use it in this script



#NOTE that this might now be incompatible with Rebecca's scripts because there are individual's with NAs for measurements


###Cleaning Ash's 2025 data####

####Cleaning + formatting the majority of the data in prep for merging####
cleaned_2025_data <- raw_2025_data %>%
  filter(notes != "Practice" | is.na(notes)) %>%
  rename_with(~str_remove(., "esrignss_")) %>%
  rename(horiz_accuracy_m = h_rms, 
         vert_accuracy_m = v_rms) %>%
  mutate(positionsourcetype = case_when(positionsourcetype == 1 ~ "user", 
                                        positionsourcetype == 3 ~ "GPS"), 
         fixtype = case_when(fixtype == 1 ~ "GPS", 
                             fixtype == 2 ~ "DGPS",
                             fixtype == 4 ~ "RTK")) %>%
  dplyr::select(bt_id, dbh_1, notes, height, leaves_coll, locality, date_coll, qubr_id_2, fruiting, Lat, Long, altitude,  horiz_accuracy_m, vert_accuracy_m, positionsourcetype, fixtype, correctionage, numsats) %>%
  rename(Metal_ID = bt_id, 
         QUBR_ID = qubr_id_2, 
         Date = date_coll, 
         DBH1 = dbh_1) %>%
  mutate(DBH2 = str_extract(notes, "(?<=DBH\\s?2: )\\d*\\.?\\d"), 
         DBH3 = str_extract(notes, "(?<=DBH\\s?3: )\\d*\\.?\\d"), 
         Metal_ID_new = str_extract(notes, "(?<=ew BT: )\\d*"), 
         Metal_ID = str_remove(Metal_ID, "^0")) %>%
  filter((str_detect(notes, "AJ", negate = T) & str_detect(notes, "Branch", negate = T)) | is.na(notes)) %>%
  filter(str_detect(Metal_ID, "AJ", negate = T) | is.na(Metal_ID)) %>%
  mutate(QUBR_ID = case_when(str_detect(QUBR_ID, "-") ~ str_replace(QUBR_ID, "-", "_"), 
                             .default = QUBR_ID)) %>%
  mutate(QUBR_ID = case_when(str_detect(notes, "_") ~ str_extract(notes, "[A-Z]*_[0-9]*"), 
                             .default = QUBR_ID))  %>%
  mutate(across(starts_with("DBH"), ~ (as.numeric(.x)/100)^2)) %>% #Squaring all of the DBH's and then summing them and then sqrting the sum is how people in the US tend to use multi-stem dbh's
  mutate(DBH_ag = sqrt(rowSums(across(starts_with("DBH")), na.rm = T))) %>%
  mutate(multistemmed = case_when(is.na(DBH2) ~ F, 
                                  !is.na(DBH2) ~ T), 
         DBH_ag = case_when(DBH_ag != 0 ~ DBH_ag, 
                            .default = NA)) %>% # adding a column that is a logical vector that describes if the tree has multiple stems or not
  select(!c(DBH1, DBH2, DBH3))

####Cleaning + formatting the bad data from LC for merging####
cleaned_LC_weird_data <- LC_weird_data %>%
  filter(str_detect(notes, "ranch", negate = T)) %>%
  select(!c(mother_tree_num, num_branches)) %>%
  separate_wider_delim(notes, ", ", names = c("height", "fruiting", "leaves_coll", "DBH1", "notes"), too_few = "align_start", too_many = "merge") %>% # the notes column contains info split by commas that I am spliting into up to 5 unique columns
  mutate(fruiting = case_when(is.na(qubr_id) & leaves_coll == "Y" ~ "Y", 
                              .default = fruiting), 
         leaves_coll = case_when(is.na(qubr_id) & leaves_coll == "Y" ~ "N", 
                                 .default = leaves_coll)) %>% # cleaning up errored inputs where I accidentally put a N for fruiting and Y for leaves collected but didn't note a QUBR_ID (so almost certainly didn't actually collect leaves) 
  rename(Metal_ID = metal_id, 
         QUBR_ID = qubr_id) %>%
  mutate(DBH2 = str_extract(notes, "(?<=DBH\\s?2: )\\d*\\.?\\d"), 
         DBH3 = str_extract(notes, "(?<=DBH\\s?3: )\\d*\\.?\\d"), 
         Metal_ID_new = str_extract(notes, "(?<=ew BT: )\\d*")) %>% # pulling out DBH2 and 3 and new metal tag ID's from the notes info
  mutate(across(starts_with("DBH"), ~ (as.numeric(.x)/100)^2)) %>% # squaring all of the DBH's and then summing them and then sqrting the sum is how people in the US tend to use multi-stem dbh's
  mutate(DBH_ag = sqrt(rowSums(across(starts_with("DBH")), na.rm = T))) %>%
  mutate(multistemmed = case_when(is.na(DBH2) ~ F, 
                                  !is.na(DBH2) ~ T), 
         DBH_ag = case_when(DBH_ag != 0 ~ DBH_ag, 
                            .default = NA)) %>% # adding a column that is a logical vector that describes if the tree has multiple stems or not
  select(!c(DBH1, DBH2, DBH3)) %>%
  mutate(Metal_ID = as.character(Metal_ID), 
         height = as.numeric(height),
         Locality = "LC") # add locality info since all data from this layer is from LC 

#check if we have any of the same trees in the weird LC data and the better data from 2025
tmp <- LC_weird_data_cleaned %>%
  filter(Metal_ID %in%  cleaned_2025_data$Metal_ID)
#we do not!


####Merging the two 2025 datasets####
cleaned_all_2025_data <- cleaned_2025_data %>%
  rename(Locality = locality) %>%
  full_join(., cleaned_LC_weird_data)

####Checking for individuals which we collected data for twice in 2025#### 
inds_w_typos <- cleaned_all_2025_data %>%
  mutate(rownum = row_number()) %>%
  group_by(Metal_ID) %>%
  filter(n() > 1) %>%
  filter(!is.na(Metal_ID)) %>%
  arrange(Metal_ID)

#There are 4 trees in SD that we collected the same data from 2x because I didn't have the data from Dana's phone downloaded onto my phone so I didn't realize we already had that data --> I want to do some comparing the sets of points to get an idea of the level of accuracy of location + height data but generally, I will consider the points taken with Dana better than the ones taken without Dana (for no particular reason)  

#Quick and dirty comparison of the trees that we sampled twice in 2025
library(geodist)
inds_w_typos_latlong <- inds_w_typos %>%
  select(Long, Lat) %>%
  rename(longitude = Long, 
         latitude = Lat)
#BT 175
geodist(inds_w_typos_latlong[1:2,], sequential = T) #fairly large distance --> uncertain as to why --> will want to take the point with more satelites probably

#BT 177
geodist(inds_w_typos_latlong[3:4,], sequential = T)
#BT 178
geodist(inds_w_typos_latlong[5:6,], sequential = T)
#BT 181
geodist(inds_w_typos_latlong[7:8,], sequential = T)

#other 3 are within 1.5 meters which is generally my resolution

#interestingly the fixtype on 11/27 is the "better" type (RTK instead of DGPS) and better resulting resolution BUT usually more satellites on 11/23 --> since it's close I'm going to choose the RKS/resolution rather than satellite number to decide which point to keep for downstream analyses 

#checking on variability in height estimations just for vibes

#BT 177
abs(inds_w_typos$height[3] -inds_w_typos$height[4])
#BT 178
abs(inds_w_typos$height[5] -inds_w_typos$height[6])
#BT 181
abs(inds_w_typos$height[7] -inds_w_typos$height[8])

#all are within 2 meters... which isn't actually that bad??? worth noting though

#Getting rid of the worse resolution duplicates (those which have the fixtype that ISN'T RTK)
cleaned_all_2025_data_final <- cleaned_all_2025_data %>%
  mutate(rownum = row_number()) %>%
  filter(rownum %notin% filter(inds_w_typos, fixtype != "RTK")$rownum)

####Merging with the data from Daniel/2023####

cleaned_2023_data_final <- fixed_field_data_processed %>%
  select(!c(Stick, Tree_pic, Env_pic, Recorder, Date, Page, Transcriber, Checked.)) %>%
  mutate(Metal_ID = as.character(Metal_ID))

cleaned_all_data <- cleaned_all_2025_data_final %>%
  full_join(., cleaned_2023_data_final, by =join_by(Metal_ID))
#originally this join didn't work bc the old data had 2 inds w/ metal tags of 181 and the new data had 2 inds with tags of 503 --> went through photos from both datasets and found that 1 was meant to be 191 and the other 504 --> fixed in respective datasets (then checked for other non-unique metal tags and found nothing)

####Checking the merged data for errors/mismatches####
#checking to see if there are any mismatched QUBR ID's between datasets
tmp <- cleaned_all_data %>%
  filter(!is.na(QUBR_ID.x)) %>%
  filter(QUBR_ID.x != QUBR_ID.y) %>%
  select(QUBR_ID.x, QUBR_ID.y, Metal_ID)
  # There are 3 here now which are inds that we collected in 2025 thinking they were new but when I went and re-examined the photos, they were clearly the same tree as an already collected individual in the 2023 data --> we only need to maintain the old ID's if there are 2


#checking to see if there are any mismatched site IDs between datasets
tmp <- cleaned_all_data %>%
  filter(Locality.x != Locality.y) %>%
  select(Metal_ID, Locality.x, Locality.y)
  #nope!


#checking to see if we have any individuals which have DBH's measured in both datasets to see if we are getting similar dbh measures
tmp <- filter(cleaned_all_data, !is.na(DBH_ag.x) & !is.na(DBH_ag.y)) %>%
  select(Metal_ID, DBH_ag.x, DBH_ag.y)
  #we do! we have 6 --> 2 = very close, one wasn't close bc I wasn't allowing DBH2 and 3 to be pulled from notes if there wasn't a space in between DBH and 2 (DBH 2 was pulled but DBH2 wasn't) --> fixed! --> another I had to check by hand --> turns out I made a typo and entered height in DBH --> fixed in datasheets --> the other 2 are individuals which we accidentally re-recorded as "new individuals" in the 2025 data --> both are pretty challenging bc 1 is a tree that's leaning over a lot so hard to know where to measure and the other was counted as a single clone in 2023 where we split it into 2 different individuals in 2025 --> in both cases we should keep the 2025 data over the 2023 data


####Cleaning + formatting the full merged dataset####

#I want am overwriting old DBHs, multistemmed, locality, lat/long BUT I don't want to overwrite old QUBR_ID's --> if we have 2 I want to maintain the older one (which is the right dataset in the merge so the .y cols with identical names)
cleaned_all_data_final <- cleaned_all_data %>% 
  mutate(QUBR_ID = case_when(is.na(QUBR_ID.x) ~ QUBR_ID.y,
                             !is.na(QUBR_ID.x) & !is.na(QUBR_ID.y) ~ QUBR_ID.y,
                             .default = QUBR_ID.x), 
         DBH_ag = case_when(is.na(DBH_ag.x) ~ DBH_ag.y,
                            .default = DBH_ag.x), 
         multistemmed = case_when(is.na(multistemmed.x) ~ multistemmed.y,
                                  .default = multistemmed.x), 
         locality = case_when(is.na(Locality.x) ~ Locality.y,
                              .default = Locality.x), 
         lat = case_when(is.na(Lat) ~ lat,
                         .default = Lat),  
         long = case_when(is.na(Long) ~ long,
                          .default = Long)) %>%
  rename(notes_new_data = notes,
         notes_old_data = Comments) %>%
  select(c(Metal_ID, QUBR_ID, locality, lat, long, altitude, fruiting, DBH_ag, multistemmed, height, Canopy_short, Canopy_long,  notes_new_data, notes_old_data, horiz_accuracy_m, vert_accuracy_m, Crown_spread, eccentricity, Canopy_area, positionsourcetype, fixtype, numsats))
write_csv(cleaned_all_data_final, "./data/all_data_merged.csv")

###Final sanity checks####
#checking for inds with no QUBR_ID --> this should be essentially no trees but instead there were 10 --> hand checked these to figure out + assign their ID's
tmp <- filter(cleaned_all_data_final, is.na(QUBR_ID)) %>%
  select(Metal_ID, QUBR_ID, DBH_ag, height, notes_new_data, notes_old_data)  
  #2 inds that have notes about not reaching leaves (so obviously no QUBR_ID)
  #1 ind which was ID'ed from photos as BT 287 but during post-processing, I realized we found the real 287 and there's no indication of what other BT this could have been so it has no BT and no QUBR_ID
  #Originally 550, 124, 384 bc these were removed from 2023 data earlier in processing bc there are missing Canopy axes/have characters in thier canopy axes --> revised that code to retain these inds
  #There were 5 other new inds where I accidentally didn't record the QUBR_ID but I could determine the QUBR_ID by looking at the date/time that the point was taken and find the IDs of the new pts taken right before and after and use the ID in between those (since we assigned QUBR_IDs in numerical order)


#figuring out how many new individuals we sampled in 2025 (that we didn't in 2023)
tmp <- cleaned_all_data_final %>%
  filter(QUBR_ID %notin% fixed_field_data_processed$QUBR_ID)
  #115 individuals!


#figuring out how many individuals weren't in the new data that were in the old data....
inds_not_found_in_2025 <- fixed_field_data_processed %>%
  filter(Metal_ID %notin% cleaned_all_2025_data_final$Metal_ID) %>%
  select(Metal_ID, QUBR_ID, lat, long)
#write.csv(inds_not_found_in_2025, "/Users/Ashley/Desktop/Hoban/Quercus_brandegeei/exploring_data/Nov_2025_GPS_data/inds_not_found_in_2025.csv")

#uploaded this csv to arcgis online so I could examine these individuals 1 by 1 next to the 2025 and 2023 data to try to identify the reason behind us missing them in 2025 --> 11 were the result of typos where the Metal tag was simply typed wrong and since we have photos of the metal tags, this was easy to identify --> 20 remain:
  #BT 645 and 647 were not sampled bc of bees 
  #BT 561 there's a good chance it's the same tree as 2308 BUT I can't be sure from the photos so assuming they are different for now
  #BT 559 there's a good chance it's the same tree as 2314 BUT I can't be sure from the photos so assuming they are different for now
  #BT 534 I think we just missed since we were in the brush
  #BTs 505, 509, 510, 511, 524 we missed I think just bc this is where we started the process + I had to stop and we switched recorders so we just didn't have our system down yet
  #BTs 167, 170, 171, 174, 176 were all behind a new fence so we couldn't resample
  #BT 120 I think we just missed
  #BT 352 just straight up didn't exist anymore --> we looked everywhere it should have been and it was just gone --> possible swept away by a storm
  #BT 233 I think we just missed
  #BT 263 wasn't in the old data that we had mapped (probably because it didn't have any value for DBH? but I'm not 100% sure) so we could resample points --> notes say it's on a high cliff so actually not surprised we didn't end up happening to resample it
  #BT 400 just straight up didn't exist anymore --> we looked everywhere it should have been and it was just gone --> possible swept away by a storm
