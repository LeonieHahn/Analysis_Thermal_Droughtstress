# Steps:
# - Data cleaning and preprocessing (Dendrometer):
#     - Bring all data to the same length
#     - Remove obvious artefacts and measurement errors manually
#     - Interpolate datagaps up to 6 hours with spline interpolation 
#     - Run jump-correction with dendRoAnalyst (see Aryal et al. 2020) for 
#       all dendrometers 
#     - Interpolate larger datagaps with network interpolation from 
#       dendRoAnalyst by using data from dendrometers of the same tree species 
#       and treatment and rerun jump-correction
#     - Compare data cleaning results visually per tree and create plots
#     - Add stem diameter from the end of experiment to have absolute tree
#       stem diameter values
#     - Temperature correction according to manuafturer (https://ecomatik.de/site/assets/files/13854/poster_-catch_the_drift-_testing_thermal_expansion_of_dendrometer_measurements.pdf)
#     - calculation of mean and sd TSD at beginning of experiment

# - Calculate different measures like TWD and GRO with dendRoAnalyst 
# - Calculate logarithmic TWD to compensate for tree size effects
# - Normalize log(TWD) to max. log(TWD) per tree species 
# - Combine TWD with WP

library(readr)
library(readxl)
library(tidyr)
library(dplyr)
library(ggplot2)
library(zoo)
library(treenetproc)
library(imputeTS)
library(dendRoAnalyst)

### Load functions
source("./Code/Functions.R")

### Read in data
# dendrometer
dendro <- read_csv("./data/dendro_avg_GH1+GH2+Treatment.csv", 
                   col_select = -1)

dendro$TIMESTAMP <- paste(dendro$TIMESTAMP, " CEST")
dendro$TIMESTAMP <- ifelse(nchar(dendro$TIMESTAMP) == 16,  
                           paste0(substr(dendro$TIMESTAMP, 1, 10),
                                  " 00:00:00 CEST"),  
                           dendro$TIMESTAMP)

dendro$TIMESTAMP <- as.POSIXct(dendro$TIMESTAMP, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")
# tree stem diameter at the end of experiment
sd <- read.csv("./data/Tree_Metainfos_Dendrometer.csv") 

# plant water potentials
wp <- read.csv("./data/Scholander_WP_data.csv")

climate <- read_excel("./data/GH2023_Climate_combined.xlsx", sheet = 1)

climate <- climate[, c(1, 21)]

# Convert to CEST
climate$DateTime <- sub(" UTC$", " CEST", climate$DateTime)

# fill in "00:00:00" for the midnight hours, since they are missing
climate$DateTime <- ifelse(nchar(climate$DateTime) == 10,
                          paste(climate$DateTime, "00:00:00"),
                          climate$DateTime)
climate$DateTime <- as.POSIXct(climate$DateTime, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")
colnames(climate)[1] <- "Datetime"


### Datacleaning and preprocessing 

## Plant water potentials
colnames(wp)[3]<- "Date"
wp <- wp[, c(3:ncol(wp))]
wp$Date <- as.POSIXct(strptime(wp$Date, format = "%d.%m.%Y"))
wp$Tree.ID <- paste0(wp$Tree.ID, "_", wp$Treatment)
wp <- pivot_wider(wp[, c(1, 2, 5, 22)], names_from = pd.md, 
                   names_prefix = "WP_", values_from = WP_average..MPa.)

## Dendrometer
dendro <- dendro[1:13511,] # remove last 2 rows with missing data

dendro_long <- pivot_longer(dendro, cols=2:41, names_to = "Tree.ID")
dendro_long$Treatment <- sub(".*(W|D).*", "\\1",  dendro_long$Tree.ID)

dendro_long <- dendro_long %>%
  mutate(Tree.Species = case_when(
    startsWith(Tree.ID, "B") ~ "Fagus sylvatica",
    startsWith(Tree.ID, "D") ~ "Pseudotsuga menziesii",
    startsWith(Tree.ID, "E") ~ "Quercus robur",
    startsWith(Tree.ID, "T") ~ "Abies alba")) %>%
  filter(Tree.Species != "Abies alba") # exclude A. alba

# Plot dendrometer data per treatment and tree species
ggplot(dendro_long)+geom_line(aes(x=TIMESTAMP, y=value, color=Tree.ID))+
  facet_grid(Tree.Species ~ Treatment)

# Bring data all on the same length

# since some dendrometer don't have data at the end, interpolate these missing
# values with the last recording to improve the recalculation of the start 
# diameter ("real" tree stem diameter only measured at the end, when 
# dendrometers were taken off)

dendro <- dendro %>%
  mutate(B3_W = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(B3_W), 
                       zoo::na.locf(B3_W), B3_W),
         B27_D = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(B27_D), 
                        zoo::na.locf(B27_D), B27_D),
         B44_W = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(B44_W), 
                        zoo::na.locf(B44_W), B44_W),
         D43_D = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(D43_D), 
                        zoo::na.locf(D43_D), D43_D),
         D7_D = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(D7_D), 
                       zoo::na.locf(D7_D), D7_D),
         E22_W = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(E22_W), 
                        zoo::na.locf(E22_W), E22_W),
         T5_D = ifelse(TIMESTAMP > "2023-09-06 16:20:00" & is.na(T5_D), 
                       zoo::na.locf(T5_D), T5_D)) 



###############Prüfen##########
# use only data from 2023-06-21 10:00:00 onwards so less NAs exist and until
# 2023-09-07 00:00:00
# which will make the data cleaning and processing more easy
dendro <- dendro %>% filter(TIMESTAMP>="2023-06-21 10:00:00" &
                              TIMESTAMP<="2023-09-07 00:00:00")


# Eliminate obviuos measurement errors and artefacts
# some dendrometers need some corrections, where the potentiometer of the 
# dendrometer was at the edge and needed to be readjusted
# D33 D also needs some corrections -> has extreme drop of values

# ggplot(dendro)+geom_line(aes(x=TIMESTAMP, y=E12_D)) +
#   xlim(as.POSIXct("2023-07-05"), as.POSIXct("2023-07-30"))
dendro$D15_W[dendro$TIMESTAMP > "2023-08-06" & 
               dendro$TIMESTAMP < "2023-08-18 12:00:00"] <- NA
dendro$D25_W[dendro$TIMESTAMP > "2023-08-11" & 
               dendro$TIMESTAMP < "2023-08-18 12:00:00"] <- NA
dendro$D39_W[dendro$TIMESTAMP > "2023-08-28" & 
               dendro$TIMESTAMP < "2023-08-31 12:00:00"] <- NA
dendro$E22_W[dendro$TIMESTAMP > "2023-08-28" & 
               dendro$TIMESTAMP < "2023-08-31 12:00:00"] <- NA
dendro$E34_W[dendro$TIMESTAMP > "2023-08-21 22:30" &
               dendro$TIMESTAMP < "2023-08-22 9:10" |
               dendro$TIMESTAMP > "2023-08-23 01:30" &
               dendro$TIMESTAMP < "2023-08-23 08:50" |
               dendro$TIMESTAMP > "2023-08-23 23:50" &
               dendro$TIMESTAMP < "2023-08-24 9:00" |
               dendro$TIMESTAMP > "2023-08-24 13:50" & 
               dendro$TIMESTAMP < "2023-08-31 12:00:00"] <- NA
dendro$E3_W[dendro$TIMESTAMP > "2023-08-24 23:00" & 
              dendro$TIMESTAMP < "2023-08-25 10:40:00" |
              dendro$E3_W >= 
              (max(dendro$E3_W, na.rm =TRUE))-4] <- NA
dendro$D33_D[dendro$TIMESTAMP > "2023-08-23 10:30:00" & 
               dendro$TIMESTAMP < "2023-08-23 15:00:00"] <- NA 
dendro$B27_D[dendro$TIMESTAMP >= "2023-07-10 10:40:00" &
               dendro$TIMESTAMP < "2023-07-15 5:00:00"] <- NA

# Spline interpolation per column for all gaps that are max 6 hours long
interpolated <- as.data.frame(na.spline(dendro[2:41], maxgap = 36))
interpolated <- cbind(dendro$TIMESTAMP, interpolated)

# Remove outliers and correct jumps with DendRoAnalyst (Aryal et al. 2020)

# rename datecolumn
colnames(interpolated)[1] <- "ts" 

# jump-correct with dendRoAnalyst
interpolated_j <- jump.locator(interpolated, v = 30)
# check how many missing values still exist in total
sum(is.na(interpolated_j)) 

# count how many missing values still exist per column
dendro_NA_count <- as.data.frame(apply(interpolated_j, 2, 
                                       function(x) sum(is.na(x))))

colnames(dendro_NA_count)[1] <- "Amount_Missing_Values"

# filter dendrometer with missing data
dendro_misdat_filt <- dendro_NA_count%>%
  filter(dendro_NA_count$Amount_Missing_Values!= 0)
dendro_misdat_names <- row.names(dendro_misdat_filt)

# Plot NA distribution per dendrometer
for(d in 1:nrow(dendro_misdat_filt )){
  dname <- rownames(dendro_misdat_filt )[d]
  NA_stats <- statsNA(interpolated[,dname ])
  NA_Plot <- ggplot_na_distribution(interpolated[, dname], title = dname)
  print(NA_Plot)
}

# To interpolate the missing dendrometer data with data from the same species
# and treatment with DendRoAnalyst (Aryal et al. 2020), we need to divide 
# the dendrometer dataframe into several dataframes according to tree species
# and treatment

df_list <- list(c("BW", "B", "W"), c("BD", "B", "D"), 
                c("EW", "E", "W"), c("ED", "E", "D"),
                c("DW", "D", "W"), c("DD", "D", "D"))
df_sorted <- lapply(df_list, function(params) 
  tst_df_creation(interpolated_j, params[1], params[2], params[3]))

# create dataframe per tree species and treatment with dendrometers that have
# missing data
df_NA_dendros_select <- lapply(df_sorted, select_dendro_with_NAs)

# create dataframe per tree species and treatment without dendrometers that have
# missing data
df_sorted_cleaned <- lapply(df_sorted, function(x) 
  remove_dendro_with_NAs(x, dendro_misdat_names))


# network interpolate and jump correct dendrometer data with data gaps and 
# combine all dendrometer data back together in one dataframe
corrected_dendro2 <- lapply(seq_along(df_NA_dendros_select), function(i) {
  process_single_dendro(df_NA_dendros_select[[i]], df_sorted_cleaned[[i]])
})

# bind the dataframes per tree species and treatment into one dataframe 
corrected_dendro_all_combined <- bind_cols(corrected_dendro2)
colnames(corrected_dendro_all_combined)[1] <- "Datetime"
corrected_dendro_all_combined <- corrected_dendro_all_combined %>% 
  dplyr::select(-contains("TIME."))

# Compare the correction of dendrometers visually

# Extract Tree.IDs without suffixes
tree_id_list <- gsub("_int$|_jcor$", "", colnames(corrected_dendro_all_combined)
                    [!grepl("_int$|_jcor$|Datetime",
                            names(corrected_dendro_all_combined))])
  
# Initialize results list
results <- list()
  
# Iterate through each tree ID
for (j in tree_id_list) {
# Define column names based on tree ID
    interpoliert_col <- paste0(j, "_int")
    jumpkorrigiert_col <- paste0(j, "_jcor")
    
    # Plot with comparison of original, interpolated, and jump-corrected data
    results[[j]] <- plot_dendrometer_comparison(
      data = corrected_dendro_all_combined,
      tree_col = j,
      png_path = "./graphics/Dendrometer_Processing/Corrections_Singletree/"
    )
}

# select dendRoAnalyst network interpolated and jump-corrected data, if no
# corrections were performed with dendRoAnalyst, keep the treentproc corrected
# data
selected_columns <- c("Datetime")  

for (id in tree_id_list) {
  jcor_col <- paste0(id, "_jcor")
  int_col <- paste0(id, "_int")
  
  if (jcor_col %in% names(corrected_dendro_all_combined)) {
    # add jcor-column if available
    selected_columns <- c(selected_columns, jcor_col)  
  } else if (id %in% names(corrected_dendro_all_combined)) {
    # add column without dendRoAnalyst corrections, 
    # if jcor-column is not available
    selected_columns <- c(selected_columns, id)  
  }
}

# Select columns from corrected_dendro_all_combined and rename columns
corrected_dendro_selected <- corrected_dendro_all_combined %>%
  dplyr::select(all_of(selected_columns))
colnames(corrected_dendro_selected) <- gsub("_jcor", "", 
                                            colnames(corrected_dendro_selected))

# Recalculate the actual diameter of the tree stems of every dendrometer 
# recording, since we only measured the actual diameter for the end of the 
# experiment (listed in sd)

sd$Tree.ID_orig <- sd$Tree.ID
sd$Tree.ID <- paste0(sd$Tree.ID, "_", sd$Treatment)
sd$sd_microm <- sd$Diameter_07092023*1000 # recalculate mm to micrometer

dendro_diam_adj <- corrected_dendro_selected

for (tree_id in sd$Tree.ID) {
  # select end diameter per Tree.ID listed in sd
  end_diameter <- sd$sd_microm[sd$Tree.ID == tree_id]
  
  # select last dendrometer recording per Tree.ID
  last_measured_value <- tail(corrected_dendro_selected[[tree_id]], n = 1)
  
  # calculate adjustment factor as the difference between the end diameter
  # from the dendrometer-recording and the actual measured stem diameter
  adjustment_factor <- end_diameter - last_measured_value
  
  # add the adjustment factor to the timeseries of the selected Tree.ID
  dendro_diam_adj[[tree_id]] <- corrected_dendro_selected[[tree_id]] + 
    adjustment_factor
}

# Plot corrected tree stem diameters per tree species and treatment
dendro_diam_adj_long <- pivot_longer(dendro_diam_adj, 
                                     cols = -Datetime, 
                                     names_to = "Tree.ID", 
                                     values_to = "TSD")

dendro_diam_adj_long$Tree.Species <- substr(dendro_diam_adj_long$Tree.ID, 1, 1)

dendro_diam_adj_long$Treatment <- sub(".*(W|D).*", "\\1", 
                                      dendro_diam_adj_long$Tree.ID)


ggplot(dendro_diam_adj_long)+
  geom_line(aes(x=Datetime, y=TSD, color = Tree.ID))+
  facet_grid(Tree.Species~Treatment)+theme_bw()+
  labs(x = "Date", y = "Tree stem diameter [µm]", color = "Tree ID", 
       title = "Tree stem diameter") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(path= 
         "./graphics/Dendrometer_Processing/Overview_Cleaned_Dendrometer/",
       filename = "Overview_TSD_cleaned.png",
       width=2560, height=1489, units = "px")


# Temperature correction of the dendrometer data 
# see also https://ecomatik.de/site/assets/files/13854/poster_-catch_the_drift-_testing_thermal_expansion_of_dendrometer_measurements.pdf
# climate data

error_per_degree <- 0.02  # µm per °C, see manufacturer information

# calculate temperature differences between actual and prev. temperature value
climate <- climate %>%
  mutate(temperature_change = Tair_min - lag(Tair_min),
         temp_error = temperature_change * error_per_degree)

dendro_diam_adj_climate <- left_join(dendro_diam_adj, climate, by = "Datetime")

dendro_diam_adj_climate %>%
  filter(is.na(temp_error))

# correct dendrometer data for temperature effects
dendro_temp_corrected  <-  dendro_diam_adj_climate %>%
#dendro_temp_corrected  <-  dendro_diam_adj_climate_20230724_1400_filled %>%
  mutate(across(.cols = -c(Datetime, temp_error, Tair_min, temperature_change), 
                ~ . - temp_error))  

# inspect when maximum temperature error was reached
ggplot() +
  geom_line(data = dendro_diam_adj_climate, 
            aes(x = Datetime, y = B2_W), color = "blue") +
  geom_line(data = dendro_temp_corrected, aes(x = Datetime, y = B2_W),
            linetype = "dashed") +
  xlim(as.POSIXct("2023-07-13 15:00"), as.POSIXct("2023-07-13 16:00")) +
  ylim(9182, 9187)

dendro_cor <- dendro_temp_corrected[, c(1:31)]
write.csv(dendro_cor, file = "./data/processed/GH_TSD_temp_cor.csv")

# calculate means and standard deviations per treatment and tree species before
# the experiment started, use first recording for this, calculate the values in
# mm
dendro_cor_long <- dendro_cor[1 ,] %>%
  pivot_longer(cols = 2:31, names_to = "Tree.ID", values_to = "TSD")
dendro_cor_long$Treatment <- sub(".*(W|D).*", "\\1",  dendro_cor_long$Tree.ID)
dendro_cor_long <- dendro_cor_long %>%
  mutate(Tree.Species = case_when(
    startsWith(Tree.ID, "B") ~ "Fagus sylvatica",
    startsWith(Tree.ID, "D") ~ "Pseudotsuga menziesii",
    startsWith(Tree.ID, "E") ~ "Quercus robur"))

dendro_start_dm_stats <- dendro_cor_long %>%
  group_by(Treatment, Tree.Species) %>%
  summarise(mean_tsd = round(mean(TSD)/1000, 1),
            sd_tsd = round(sd(TSD)/1000, 2))

write.csv(dendro_start_dm_stats, 
          file = "./data/processed/GH_TSD_start_diam_stats.csv")

# define start and end dates of dataframe and recalculate into DOY
# extract first date of dataframe and recalculate into DOY
first_row <- head(dendro_cor, 1)
first_date <- as.Date(first_row$Datetime)
first_day_of_year <- as.numeric(format(first_date, "%j"))

# extract last date of dataframe and recalculate into DOY
last_row <- tail(dendro_cor, 1)
last_date <- as.Date(last_row$Datetime)
last_day_of_year <- as.numeric(format(last_date, "%j"))

# calculate metrics like TWD and GRO from dendRoAnalyst and create overview plot 
# for every tree 

# create empty dataframes to store results in
ZG_cycle_DF <- data.frame()
ZG_phase_DF <- data.frame()
daily_data_DF <- data.frame()

for(i in 2:length(dendro_cor)){
  calc_plot_DA(i, dendro_cor, 2023, 
               "./graphics/Dendrometer_Processing/TWD/DA-TWD_Singletree/") }

# replace Inf values with NA
daily_data_DF <- as.data.frame(lapply(daily_data_DF, replace_inf))

write.csv(daily_data_DF, 
          "./data/processed/daily_data_DF.csv")
write.csv(ZG_phase_DF, 
          "./data/processed/ZG_phase_DF.csv")
write.csv(ZG_cycle_DF, 
          "./data/processed/ZG_cycle_DF.csv")  

# ZG_phase_DF <- read_csv( 
#                    "./../data/2_DA_results/ZG_phase_DF.csv", col_select = -1)

# change tree species and treatment labels
ZG_phase_DF$Tree.Species <- factor(ZG_phase_DF$Tree.Species,
                                   levels = c("B", "E", "D", "T"),
                                   labels = c("B" = "Fagus sylvatica",
                                              "E" = "Quercus robur",
                                              "D" = "Pseudotsuga menziesii",
                                              "T" = "Abies alba"))

ZG_phase_DF$Treatment <- factor(ZG_phase_DF$Treatment, 
                                labels = c("D" = "Droughted", 
                                           "W" = "Watered"))

daily_data_DF$Tree.Species <- factor(daily_data_DF$Tree.Species,
                                     levels = c("B", "E", "D", "T"),
                                     labels = c("B" = "Fagus sylvatica",
                                                "E" = "Quercus robur",
                                                "D" = "Pseudotsuga menziesii",
                                                "T" = "Abies alba"))

daily_data_DF$Treatment <- factor(daily_data_DF$Treatment, 
                                  labels = c("D" = "Droughted", 
                                             "W" = "Watered"))

# normalize GRO to 0 at the first measurement for comparison
ZG_phase_DF <- ZG_phase_DF %>%
  group_by(Tree.ID) %>%                     
  mutate(GRO_norm = GRO - first(GRO)) %>%         
  ungroup() 

# create overview plots per treatment and tree species without abies alba, 
# since they are further excluded from the analysis due to lacking
# strong drought stress reactions

# calculate min, max and mean of dm, TWD and GRO per day, treatment and 
# tree  species
ZG_phase_DF_stats <- ZG_phase_DF %>%
  group_by(TIME, Treatment, Tree.Species) %>%
  summarise(Min_TSD = min(dm),
            Max_TSD = max(dm),
            Mean_TSD = mean(dm),
            Min_TWD = min(TWD),
            Max_TWD = max(TWD),
            Mean_TWD = mean(TWD),
            Min_GRO = min(GRO_norm),
            Max_GRO = max(GRO_norm),
            Mean_GRO = mean(GRO_norm))

# Plots for TSD per tree
ggplot(ZG_phase_DF) +
  geom_line(aes(x=TIME, y=dm, color=Tree.ID))+
  facet_grid(Tree.Species ~ Treatment) + theme_bw()+
  xlab("Date") + ylab("TSD [µm]")
ggsave(path = "./graphics/Dendrometer_Processing/TSD/",
       filename = "Overview_TSD_Singletree.png",
       width = 2560, height = 1550, units = "px")

# plot daily mean TWD and its ranges
ggplot(ZG_phase_DF_stats) +
  geom_ribbon(aes(x = TIME, ymin = Min_TSD, ymax = Max_TSD, 
                  fill = Treatment), alpha = 0.3)+
  geom_line(aes(x = TIME,
                y = Mean_TSD, 
                color = Treatment))+
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue")) +
  facet_grid(Tree.Species ~ .) + theme_bw()+
  xlab("Date") + ylab("TSD [µm]") +
  ggtitle("Mean tree stem diameter (TSD) and its ranges") +
  theme(plot.title = element_text(hjust = 0.5))
ggsave(
  path = "./graphics/Dendrometer_Processing/TSD/",
  filename = "Overview_MeanMinMax_TSD.png",
  width = 2560, height = 1550, units = "px")

# Plots for TWD per tree
ggplot(ZG_phase_DF) +
  geom_line(aes(x=TIME, y=TWD, color=Tree.ID)) +
  facet_grid(Tree.Species ~ Treatment) + theme_bw() +
  xlab("Date") + ylab("TWD [µm]")
ggsave(path = "./graphics/Dendrometer_Processing/TWD/DA-TWD_Overview/",
       filename = "Overview_TWD_Singletree.png",
       width = 2560, height = 1550, units = "px")

# calculate logarithmic TWD compensate for tree size effects
ZG_phase_DF$TWD_log <- log((ZG_phase_DF$TWD+1), base=10)

# calculate mean, min and max log(TWD) per treatment and tree species
ZG_phase_log_stats <- ZG_phase_DF %>%
  group_by(Date = as.Date(TIME, tz = "Europe/Berlin"), 
           Tree.Species, Treatment) %>%
  summarise(Min_TWD_log = min(TWD_log, na.rm=TRUE),
            Max_TWD_log = max(TWD_log, na.rm=TRUE),
            Mean_TWD_log = mean(TWD_log, na.rm=TRUE))

# plot mean, min and max log(TWD) per treatment and tree species
ggplot(ZG_phase_log_stats) +
  geom_ribbon(aes(x = Date, ymin = Min_TWD_log, ymax = Max_TWD_log, 
                  fill = Treatment), alpha = 0.3)+
  geom_line(aes(x = Date,
                y = Mean_TWD_log, 
                color = Treatment))+
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  facet_grid(Tree.Species ~ .) + theme_bw()+
  xlab("Date") + ylab((expression(Log[10](TWD)))) +
  ggtitle(expression(Mean~log[10](TWD)~and~ranges)) +
  theme(plot.title = element_text(hjust = 0.5))

# select max log(TWD) per tree species and store results in a dataframe
species_list <- c("Fagus sylvatica", "Quercus robur", "Pseudotsuga menziesii")

max_data_log <- bind_rows(lapply(species_list, get_max_TWD_log_info))

# Calculate relative log(TWD) by normalizing the log(TWD) to the max log(TWD)
# observed per species 
TWD_rel_df <- bind_rows(lapply(species_list, normalize_TWD))
write.csv(TWD_rel_df, "./data/processed/ZG_phase_rel_logTWD.csv")

# Plot rel. log(TWD) per treatment and tree species
ggplot(TWD_rel_df) +
  geom_line(aes(x = TIME, y = rel_TWD_log, color = Tree.ID)) +
  facet_grid(Tree.Species~Treatment) + theme_bw()+
  xlab("Date") + ylab("Rel. TWD")
ggsave(path=
         "./graphics/Dendrometer_Processing/TWD/DA-TWD_Overview/",
       filename = "Overview_Rel_logTWD_Singletree.png",
       width = 2560, height = 1550, units = "px")

# calculate mins, maxs and means of rel. log(TWD)
TWD_rel_stats <- TWD_rel_df %>%
  group_by(Date = as.Date(TIME, tz = "Europe/Berlin"), 
           Treatment, Tree.Species) %>%
  summarise(Min_rel_TWD = min(rel_TWD_log, na.rm=TRUE),
            Max_rel_TWD = max(rel_TWD_log, na.rm=TRUE),
            Mean_rel_TWD = mean(rel_TWD_log, na.rm=TRUE))

# Plot mean rel. log(TWD) and ranges per treatment and tree species
ggplot(TWD_rel_stats)+
  geom_ribbon(aes(x = as.POSIXct(Date, tz = "Europe/Berlin"), 
                  ymin = Min_rel_TWD, ymax = Max_rel_TWD, 
                  fill = Treatment), alpha = 0.3)+
  geom_line(aes(x = as.POSIXct(Date, tz = "Europe/Berlin"), 
                y = Mean_rel_TWD, color = Treatment)) +
  geom_vline(aes(xintercept = as.POSIXct("2023-06-26", tz = "Europe/Berlin")),
             color = "chartreuse4", linewidth = 1) +
  facet_grid(Tree.Species ~ .) + theme_bw() +
  scale_fill_manual(values = c("red", "blue")) +
  scale_color_manual(values = c("red", "blue")) +
  scale_x_datetime(date_labels = "%d %b") +
  xlab("2023") + ylab("Rel. TWD") +
  theme(plot.title = element_text(hjust = 0.5),
        legend.position = "bottom")
ggsave(path=
         "./graphics/Dendrometer_Processing/TWD/DA-TWD_Overview/",
       filename = "Overview_MeanMinMax_Rel_logTWD.png",
       width = 2560, height = 1950, units = "px")

# select min. and max. rel. log(TWD) per tree
min_max_twd_per_day_and_tree <- TWD_rel_df %>%
  group_by(Date = as.Date(TIME, tz = "Europe/Berlin"), Tree.ID) %>%
  summarise(Min_rel_TWD = min(rel_TWD_log),
            Max_rel_TWD = max(rel_TWD_log))

# combine TWD with WP
combined<- left_join(min_max_twd_per_day_and_tree, wp,
                      by = c("Date", "Tree.ID"))
# add treatment and tree species
combined$Treatment <- gsub(".*([DW]).*", "\\1", combined$Tree.ID)
combined$Tree.Species <- str_sub(combined$Tree.ID, 1, 1)

write.csv(combined, "./data/processed/TWD_WP_combined.csv")


#################################################
# =============================================================================
# - Calculate irrigation thresholds ITlow and ITup per tree in the greenhouse 
# - Visualize the thresholds in combination with environmental measurements 
#   per tree and as overview heatmaps for all trees
# - Select max. irrigation demand for when the trees still recovered back to 
#   no irrigation demand (for the droughted trees only)
# =============================================================================

# library(readr)
# library(readxl)
# library(dplyr)
# library(purrr)
# library(ggplot2)
# library(patchwork)
# library(hrbrthemes)
# library(tidyr)
# library(ggtext)
# library(zoo)
# library(purrr)

# Load functions
#source("./code/Functions.R")

# model_results <- read_csv(
#   "./data/GH/results/log_model_results/Log_model_results.csv", 
#   col_select = -1)
# 
# threshold_perc <- read_csv(
#   "./data/GH/results/log_model_results/p12thres_percentages.csv", 
#   col_select = -1)

# table with output from DendRoAnalyst (phase.zg) and rel. log TwD
rel_TWD_df<-read_csv(
  "./data/processed/ZG_phase_rel_logTWD.csv",
  col_select = c(2:12))

rel_TWD_df$TIME <- sub(" UTC$", " CEST", rel_TWD_df$TIME)

# fill in "00:00:00" for the midnight hours, since they are missing
rel_TWD_df$TIME <- ifelse(nchar(rel_TWD_df$TIME) == 10,
                          paste(rel_TWD_df$TIME, "00:00:00"),
                          rel_TWD_df$TIME)
rel_TWD_df$TIME <- as.POSIXct(rel_TWD_df$TIME, format="%Y-%m-%d %H:%M:%S", 
                              tz = "Europe/Berlin")
rel_TWD_df$Date <- as.Date(rel_TWD_df$TIME, tz = "Europe/Berlin")

# table with output from DendRoAnalyst (daily.data)
# daily_df <- read_csv("./data/GH/input_data/daily_data_DF.csv", 
#                      col_select = -1)
# daily_df$DATE <- as.Date(daily_df$DATE, tz = "Europe/Berlin")
# colnames(daily_df)[1] <- "Date"

# climate data
climate <- read_excel(
  "./data/GH2023_Climate_combined.xlsx",
  sheet = 2)
colnames(climate) <- c("Date", "VPD_adj_dailyavg_kPa", "VPD_adj_dailymax_kPa",
                       "VPD_adj_dailymin_kPa", "n_dailyavg", 
                       "VPD_meas_dailyavg_kPa", "VPD_meas_dailymax_kPa",
                       "VPD_meas_dailymin_kPa")
climate$Date <- as.Date(climate$Date, tz = "Europe/Berlin")

# vol. soil water content 
vwc <- read_csv(
  "./data/GH_VWC_dendrotrees.csv", 
  col_select = -1)
vwc$Tree.ID <- paste0(vwc$Tree.ID, "_", vwc$Treatment)
vwc$Date <- as.Date(vwc$Date, tz = "Europe/Berlin")

# table with TWD and canopy water potentials
combined <- read_csv("./data/processed/TWD_WP_combined.csv", 
                     col_select = -1)
combined$Date <- as.Date(combined$Date, tz = "Europe/Berlin")

# # select best model for predicting TWD according to R²
# best_model <- model_results %>%
#   group_by(Tree.Species) %>%
#   slice_max(order_by = R2_value, n = 1) %>%
#   ungroup()
# 
# # select rel. TWD percentages (= irrigation demand) according to best R²
# threshold_perc_sel <- threshold_perc %>%
#   group_by(Tree.Species) %>%
#   slice_max(order_by = R2_value, n = 1) %>%
#   ungroup()
# 
# # select p12 thresholds according to tree species and best R²
# thres_B <- filter(best_model, Tree.Species == "B")$p12_thres
# thres_E <- filter(best_model, Tree.Species == "E")$p12_thres
# thres_D <- filter(best_model, Tree.Species == "D")$p12_thres

# combine dataframes
combined <- left_join(rel_TWD_df, climate[, 1:4], by = "Date")
combined <- left_join(combined, vwc[, c(1:3)], by = c("Date", "Tree.ID")) 

# add min. rel TWD per day and Tree.ID 
combined  <- combined  %>%
  group_by(Tree.ID, Date) %>%
  mutate(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) %>%
  ungroup()


# # calculate when min TWD > 0
MinTWD_not0 <- combined %>%
  group_by(Date = Date, Tree.ID, Treatment, Tree.Species) %>%
  summarise(Min_TWD = min(TWD)) %>%
  filter(Min_TWD != 0)
MinTWD_not0$MinTWDnot0_Date <- as.Date(MinTWD_not0$Date)

# ITup - calculate when 12% Xylem conductivity loss corresponding min. TWD value 
# threshold is crossed
# p12_thres_crossed <- combined %>%
#   mutate(Date_rel_minTWD_thres_crossed = case_when(
#     Tree.Species == "Fagus sylvatica" & min_rel_TWD >= thres_B ~ as.Date(Date),
#     Tree.Species == "Quercus robur" & min_rel_TWD >= thres_E ~ as.Date(Date),
#     Tree.Species == "Pseudotsuga menziesii" & 
#       min_rel_TWD >= thres_D ~ as.Date(Date),
#     TRUE ~ as.Date(NA)
#   )) %>%
#   select(Tree.ID, Date, Date_rel_minTWD_thres_crossed) %>%
#   unique()

#  calculate last day of GRO and the day after last GRO
combined_sorted <- combined %>%
  arrange(Tree.ID, Date)
lastGRO <- combined_sorted %>%
  group_by(Tree.ID) %>%
  filter(Phases == "2") %>%
  slice_tail(n = 1)
lastGRO$FirstDayAfterLastGRO <- as.Date(lastGRO$Date) + 1
lastGRO$LastGRODate <- as.Date(lastGRO$Date)

FirstDayAfterLastGRO_df <- lastGRO[, c(6, 18, 18)]
colnames(FirstDayAfterLastGRO_df)[3] <- "Date"

# set "2023-09-07 and following dates to NA since this Date does not exist in 
# dataframe (data only exists until "2023-09-06")
FirstDayAfterLastGRO_df <- FirstDayAfterLastGRO_df %>%
  mutate(FirstDayAfterLastGRO = if_else(
    FirstDayAfterLastGRO >= as.Date("2023-09-07"), 
    as.Date(NA), 
    FirstDayAfterLastGRO
  ))

# merge results in one dataframe
results <- reduce(list(combined, MinTWD_not0[, c(1, 2, 6)], 
                       FirstDayAfterLastGRO_df),
                  left_join, by = c("Date", "Tree.ID")) 

# fill in FirstDayAfterLastGRO all FirstDayAfterLastGRO-dates after 
# that occurred
results <- results %>%
  group_by(Tree.ID) %>%
  mutate(FirstDayAfterLastGRO = if_else(is.na(FirstDayAfterLastGRO), 
                                        lag(FirstDayAfterLastGRO), 
                                        FirstDayAfterLastGRO)) %>%
  fill(FirstDayAfterLastGRO, .direction = "down") %>%
  ungroup()

results <- results %>%
  dplyr::select(Tree.ID, Treatment, Tree.Species, Date, VPD_adj_dailyavg_kPa,
                VPD_adj_dailymax_kPa,	VPD_adj_dailymin_kPa,	VWC_mean,
                min_rel_TWD, FirstDayAfterLastGRO) %>%
  unique()


write.csv(results, "./data/processed/GH_Dendro_2023.csv", row.names = FALSE)

# # Round threshold_perc_sel$rounded_rel_TWD and change
# # add color to percentage value
# threshold_perc_sel <- threshold_perc_sel %>%
#   group_by(Tree.Species) %>%
#   mutate(
#     rounded_rel_TWD = round(p12_thres_scaled, 8) 
#   )

# # tree species and corresponding tree species abbreviation in df_base_perc
# species_mapping <- tibble(
#   results_species = c("Fagus sylvatica", "Quercus robur", 
#                       "Pseudotsuga menziesii"),
#   threshold_species = c("B", "E", "D")
# )
# 
# # apply calc_thres_perc on all tree species and combine results
# results <- species_mapping %>%
#   pmap_dfr(function(results_species, threshold_species) {
#     calc_thres_perc(
#       df_results = results, 
#       df_base_perc = threshold_perc_sel, 
#       treespecies_l = results_species, 
#       treespecies_s = threshold_species
#     )
#   })

# results$Tree.Species.y <- NULL
# colnames(results)[8] <- "Tree.Species"

# # create new percentage column, where 0 = no threshold crossed,
# # 1-99 threshold ITlow is crossed, 100 = threshold ITup is reached and 
# # 101 = threshold ITup is crossed
# results <- results %>%
#   mutate(
#     Percentage_w_cros = case_when(
#       Percentage == 100 &
#         !is.na("Date_rel_minTWD_thres_crossed") ~ 101, 
#       TRUE ~ Percentage                               
#     )
#   )
# 
# results <- results %>%
#   mutate(
#     Percentage_color = case_when(
#       Percentage_w_cros == 0 ~ "lightblue",      
#       Percentage_w_cros == 101 ~ "darkmagenta",
#       Percentage_w_cros > 0 & Percentage_w_cros <= 100 ~ scales::col_numeric(
#         palette = c("goldenrod2", "mediumorchid3"),
#         domain = c(1, 100)
#       )(Percentage_w_cros),  # Continous color scheme for values from 1-99
#       TRUE ~ "black"                  
#     )
#   )
# 
# # Create Demo color palette
# palette_df <- data.frame(Percentage_w_cros = -10:110) %>%
#   mutate(
#     Percentage_color = case_when(
#       Percentage_w_cros <= 0 ~ "lightblue",      
#       Percentage_w_cros >= 101 ~ "darkmagenta",
#       Percentage_w_cros > 0 & Percentage_w_cros <= 100 ~ scales::col_numeric(
#         palette = c("goldenrod2", "mediumorchid3"),
#         domain = c(1, 100)
#       )(Percentage_w_cros),  # Continous color scheme for values from 1-99
#       TRUE ~ "black"                  
#     )
#   )
# 
# ggplot(palette_df, aes(x = Percentage_w_cros, fill = Percentage_color)) +
#   geom_tile(aes(y = 1)) +
#   scale_fill_identity() +
#   scale_x_continuous(name = "%", 
#                      breaks = c(-5, 1, 25, 50, 75, 100, 110),
#                      labels = c( "0", "1", "25", "50", 
#                                  "75", "100", ">100")) +
#   theme_minimal() +
#   theme(
#     axis.ticks.y = element_blank(),       
#     axis.text.y = element_blank(),       
#     axis.title.y = element_blank(),       
#     axis.text.x = element_text(size = 22), 
#     axis.title.x = element_text(size = 24)  
#   )
# 
# ggsave(path = "./graphics/GH/Irrigation_Thresholds/",
#        filename = "Color_Pallete_Irrigation_Threshold_Percentages.png",
#        width = 2560, height = 1489, units = "px")
# 
# 
# # reverse plot
# ggplot(palette_df, aes(x = Percentage_w_cros, fill = Percentage_color)) +
#   geom_tile(aes(y = 1)) +
#   scale_fill_identity() +
#   scale_x_reverse( 
#     name = "%", 
#     breaks = c(110, 100, 75, 50, 25, 1, -5), 
#     labels = c(">100", "100", "75", "50", "25", "1", "0")  
#   ) +
#   theme_minimal() +
#   theme(
#     axis.ticks.y = element_blank(),       
#     axis.text.y = element_blank(),       
#     axis.title.y = element_blank(),       
#     axis.text.x = element_text(size = 22),  
#     axis.title.x = element_text(size = 24) 
#   )
# 
# ggsave(path =
#          "./graphics/GH/Irrigation_Thresholds/",
#        filename = "Color_Pallete_Irrigation_Threshold_Percentages_rev.png",
#        width = 2560, height = 1489, units = "px")
# 
# # plot irrigation thresholds per tree
# trees <- as.list(unique(results$Tree.ID))
# 
# for (i in trees){
#   Tree.ID_Nr <- i
#   print(Tree.ID_Nr)
#   selected_tree_species <- substr(Tree.ID_Nr, 1, 1)
#   rel_TWD_max_thres <- ifelse(selected_tree_species == "E", thres_E,
#                               ifelse(selected_tree_species == "B", thres_B,
#                                      thres_D))
#   filtered_data <- filter(results, Tree.ID == Tree.ID_Nr)
#   filtered_data$Date <- as.POSIXct(filtered_data$Date)
#   plot_irrigation_thresholds(filtered_data, 
#                              "./graphics/GH/Irrigation_Thresholds/Irrigation_Thresholds_Singletree/")
# }
# 
# # create stacked area charts with ggplot to visualize the dates when the 
# # irrigation thresholds occurred
# thresholds_daily <- results[, c(6:8, 12, 13:21, 23:24)] %>%
#   group_by(Date, Tree.Species, Treatment)
# thresholds_daily_unique <- thresholds_daily %>% distinct()

# calculate mean SWC, standard deviation and rel. SWC per treatment 
# and tree species for plotting later on
# vwc_stats <- results[, c(7, 8, 12, 16)] %>%
#   group_by(Date, Treatment, Tree.Species) %>%
#   summarise(SWC_mean_treat = mean(VWC_mean, na.rm = TRUE),
#             SWC_sd_treat = sd(VWC_mean, na.rm = TRUE))
# 
# # create new column with information that no threshold was reached 
# # (but we would have had data to differentiate NAs from thresholds not reached)
# thresholds_daily_unique <- thresholds_daily_unique %>%
#   mutate(Date_nothres = case_when(
#     is.na(MinTWDnot0_Date) & 
#       is.na(FirstDayAfterLastGRO) & 
#       is.na(Date_rel_minTWD_thres_crossed) ~  as.Date(Date),
#     FALSE ~ as.Date(NA)))
# 
# # create heatmaps with irrigation thresholds per tree
# heatmap_data <- thresholds_daily_unique %>%
#   pivot_longer(
#     cols = c(MinTWDnot0_Date,  
#              Date_rel_minTWD_thres_crossed, 
#              Date_nothres),
#     names_to = "Threshold",
#     values_to = "ThresholdDate"
#   ) %>%
#   filter(!is.na(ThresholdDate)) %>%  # Filter dates where threshold is crossed
#   mutate(ThresholdDate = as.Date(ThresholdDate),  
#          Date = as.Date(Date))                    
# 
# # Prioritize irrigation thresholds and keep only threshold with highest priority
# heatmap_data_filtered <- heatmap_data %>%
#   mutate(Threshold_Priority = case_when(
#     Threshold == "Date_rel_minTWD_thres_crossed" ~ 2,
#     Threshold == "MinTWDnot0_Date" ~ 1,
#     Threshold == "Date_nothres" ~ 0
#   )) %>%
#   group_by(Tree.ID, Date) %>%
#   slice_max(order_by = Threshold_Priority, with_ties = FALSE) %>%  
#   ungroup()
# 
# heatmap_data_filtered$Date <- as.POSIXct(heatmap_data_filtered$Date, 
#                                          format = "%Y-%m-%d")
# 
# # sort Tree.IDs descending after irrigation threshold percentages
# threshold_summary <- heatmap_data_filtered %>%
#   group_by(Tree.ID) %>%
#   summarise(
#     High_Percentage_Count = sum(Percentage_w_cros == "101", na.rm = TRUE),
#     Mean_Percentage = mean(as.numeric(Percentage_w_cros), na.rm = TRUE)
#   ) %>%
#   arrange(desc(High_Percentage_Count), desc(Mean_Percentage))
# 
# heatmap_data_sorted <- heatmap_data_filtered %>%
#   mutate(Tree.ID = factor(Tree.ID, levels = threshold_summary$Tree.ID)) %>%
#   arrange(Tree.ID, Date)
# 
# FirstDayAfterLastGRO_df2 <- FirstDayAfterLastGRO_df %>%
#   mutate(
#     Tree.Species = case_when(
#       grepl("B", Tree.ID) ~ "Fagus sylvatica",
#       grepl("E", Tree.ID) ~ "Quercus robur",
#       TRUE ~ "Pseudotsuga menziesii"),
#     Treatment = ifelse(grepl("W", Tree.ID), "Watered", "Droughted")
#   )
# 
# x_axis_limits <- c(as.POSIXct("2023-06-15"), as.POSIXct("2023-09-06"))
# x_axis_settings <- scale_x_datetime(
#   date_labels = "%d %b" ,
#   limits = x_axis_limits
# )
# 
# # Scaling factor for distributing SWC evenly on the 5 trees in the plot
# swc_scaling_factor <- 0.1
# 
# heatmap_data_sorted <- heatmap_data_sorted %>%
#   mutate(Tree.ID_num = as.numeric(Tree.ID)) %>%
#   group_by(Tree.Species, Treatment) %>%
#   mutate(Tree.ID_num_cat = as.numeric(factor(Tree.ID))) %>%
#   ungroup()  
# 
# FirstDayAfterLastGRO_df2 <- left_join(FirstDayAfterLastGRO_df2,
#                                       heatmap_data_sorted[, c(1, 13:18)], 
#                                       by = "Tree.ID")
# 
# FirstDayAfterLastGRO_df2 <- unique(FirstDayAfterLastGRO_df2)
# 
# 
# # Create Heatmaps with irrigation thresholds per tree and additional plots with
# # the environmental conditions
# VPD_plot_1 <- ggplot(filter(combined, Date >= "2023-06-19")) +
#   geom_ribbon(aes(x = as.POSIXct(Date), 
#                   ymin = VPD_adj_dailymin_kPa,
#                   ymax = VPD_adj_dailymax_kPa), 
#               fill ="black", alpha = 0.2) +
#   geom_line(aes(x = as.POSIXct(Date), y = VPD_adj_dailyavg_kPa), 
#             color = "grey35", linewidth = 1, linetype = "dashed") +
#   x_axis_settings +
#   scale_y_continuous(
#     name = "VPD [kPa]",
#     breaks = c(1, 3, 5),
#     labels = c("1", "3", "5"),
#     sec.axis = sec_axis(
#       ~ .,  
#       name = "")) +
#   theme_bw() +
#   theme(
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.title.y.left = element_text(size = 16),  
#     axis.title.x = element_blank(),
#     axis.text.y.left = element_text(size = 14),
#     axis.title.y.right = element_blank(),
#     axis.text.y.right = element_blank(),
#     axis.ticks.y.right = element_blank()
#   )
# 
# 
# VPD_plot_2 <- ggplot(filter(combined, Date >= "2023-06-19")) +
#   geom_ribbon(aes(x = as.POSIXct(Date), 
#                   ymin = VPD_adj_dailymin_kPa,
#                   ymax = VPD_adj_dailymax_kPa), 
#               fill ="black", alpha = 0.2) +
#   geom_line(aes(x = as.POSIXct(Date), y = VPD_adj_dailyavg_kPa), 
#             color = "grey35", linewidth = 1, linetype = "dashed") +
#   x_axis_settings +
#   scale_y_continuous(
#     name = "VPD [kPa]",
#     breaks = c(1, 3, 5),
#     labels = c("1", "3", "5"),
#     sec.axis = sec_axis(
#       ~ ., 
#       name = "")) +
#   theme_bw() +
#   theme(
#     strip.text.x = element_blank(),
#     strip.text.y = element_blank(),
#     axis.text.x = element_text(size = 14),
#     axis.title.y.left = element_blank(), 
#     axis.title.x = element_blank(),
#     axis.text.y.left = element_blank(),
#     axis.title.y.right = element_blank(),
#     axis.text.y.right = element_blank(),
#     axis.ticks.y = element_blank()
#   )
# 
# heatmap_plot_D <- ggplot() +
#   geom_tile(
#     data = filter(heatmap_data_sorted, Treatment == "Droughted"), 
#     aes(x = as.POSIXct(Date), y = Tree.ID_num_cat, fill = Percentage_color), 
#     color = "white") +
#   scale_fill_identity() +  
#   # add scaled SWC
#   geom_ribbon(data = filter(vwc_stats, Treatment == "Droughted"),
#               aes(x = as.POSIXct(Date),
#                   ymin =(SWC_mean_treat - SWC_sd_treat) * swc_scaling_factor +1,
#                   ymax = (SWC_mean_treat + SWC_sd_treat) * swc_scaling_factor +1),
#               alpha = 0.5) +
#   geom_line(
#     data = filter(vwc_stats, Treatment == "Droughted"), 
#     aes(
#       x = as.POSIXct(Date), 
#       y = SWC_mean_treat * swc_scaling_factor + 1,  
#       group = 1
#     ),  
#     color = "black",
#     linewidth = 1
#   ) +
#   # add FirstDayAfterLastGRO
#   geom_point(
#     data = filter(FirstDayAfterLastGRO_df2, Treatment == "Droughted"),
#     aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
#         y = Tree.ID_num_cat),
#     shape = 19,
#     color = "#00008b",
#     size = 2
#   ) +
#   scale_y_continuous(
#     breaks =  unique(filter(heatmap_data_sorted, 
#                             Treatment == "Droughted")$Tree.ID_num_cat),
#     labels = 
#       unique(filter(heatmap_data_sorted, 
#                     Treatment == "Droughted")$Tree.ID_num_cat),
#     sec.axis = sec_axis(
#       ~ (. - 1) / swc_scaling_factor,  
#       name = "SWC [%]",
#       breaks = c(0, 20, 40, 60),
#       labels = c(0,20, 40, 60)
#     )
#   ) +
#   geom_label(
#     data = filter(heatmap_data_sorted, 
#                   Treatment == "Droughted"),
#     aes(x = min(as.POSIXct(Date)) - (86400 * 5), y = Tree.ID_num_cat,
#         label = gsub("_", "", Tree.ID)),
#     color = "black",
#     fill = "white",  
#     size = 5,
#     fontface = "plain",
#     label.padding = unit(0.2, "lines"),  
#     label.r = unit(0.2, "lines") 
#   ) +
#   facet_grid(Tree.Species ~ Treatment, scales = "free", space = "free",
#              switch = 'y',
#              labeller = labeller(Tree.Species = as_labeller(c(
#                "Fagus sylvatica" = "F. sylvatica",
#                "Pseudotsuga menziesii" = "P. menziesii",
#                "Quercus robur" = "Q. robur"
#              )))) +
#   theme_bw() +
#   x_axis_settings +
#   theme(
#     legend.position = "none",
#     strip.text.y = element_text(size = 16),
#     strip.text.x = element_text(size = 16),
#     axis.text.y = element_blank(),
#     axis.text.y.right = element_blank(),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(), 
#     axis.title.y.left = element_blank(), 
#     axis.ticks.y.left = element_blank(),
#     axis.ticks.y.right = element_blank(),  
#     axis.ticks.x = element_blank(), 
#     axis.text.y.left = element_blank(),  
#     axis.title.y.right = element_blank()
#   )
# 
# print(heatmap_plot_D)
# 
# heatmap_plot_W <- ggplot() +
#   geom_tile(
#     data = filter(heatmap_data_sorted, Treatment == "Watered"), 
#     aes(x = as.POSIXct(Date), y = Tree.ID_num_cat, fill = Percentage_color), 
#     color = "white"
#   ) +
#   scale_fill_identity() +  
#   # add scaled SWC
#   geom_ribbon(data = filter(vwc_stats, Treatment == "Watered"),
#               aes(x = as.POSIXct(Date),
#                   ymin =(SWC_mean_treat - SWC_sd_treat) * swc_scaling_factor + 1,
#                   ymax = (SWC_mean_treat + SWC_sd_treat) * swc_scaling_factor + 1),
#               alpha = 0.5) +
#   geom_line(
#     data = filter(vwc_stats, Treatment == "Watered"), 
#     aes(
#       x = as.POSIXct(Date), 
#       y = SWC_mean_treat * swc_scaling_factor + 1,  
#       group = 1
#     ),  
#     color = "black",
#     linewidth = 1
#   ) +
#   # add FirstDayAfterLastGRO
#   geom_point(
#     data = filter(FirstDayAfterLastGRO_df2, Treatment == "Watered"),
#     aes(x = as.POSIXct(FirstDayAfterLastGRO, format = "%Y-%m-%d"), 
#         y = Tree.ID_num_cat),
#     shape = 19,
#     color = "#00008b",
#     size = 2
#   ) +
#   scale_y_continuous(
#     breaks =  unique(filter(heatmap_data_sorted, 
#                             Treatment == "Watered")$Tree.ID_num_cat),
#     labels = 
#       unique(filter(heatmap_data_sorted, 
#                     Treatment == "Watered")$Tree.ID_num_cat),
#     sec.axis = sec_axis(
#       ~ (. - 1) / swc_scaling_factor,  
#       name = "SWC [%]",
#       breaks = c(0, 20, 40, 60),
#       labels = c(0, 20, 40, 60)
#     )
#   ) +
#   geom_label(
#     data = filter(heatmap_data_sorted, Treatment == "Watered"),
#     aes(x = min(Date, na.rm = TRUE) - (86400 * 5), y = Tree.ID_num_cat,
#         label = gsub("_", "", Tree.ID)),
#     color = "black",
#     fill = "white",  
#     size = 5,
#     fontface = "plain",
#     label.padding = unit(0.2, "lines"),  
#     label.r = unit(0.2, "lines")  
#     
#   ) +
#   facet_grid(Tree.Species ~ Treatment, scales = "free", space = "free",
#              switch = 'y',
#              labeller = labeller(Tree.Species = as_labeller(c(
#                "Fagus sylvatica" = "F. sylvatica",
#                "Pseudotsuga menziesii" = "P. menziesii",
#                "Quercus robur" = "Q. robur"
#              )))) +
#   theme_bw() +
#   x_axis_settings +
#   theme(
#     legend.position = "none",
#     strip.text.y = element_blank(),
#     strip.text.x = element_text(size = 16),
#     axis.text.y = element_text(size = 14),
#     axis.text.y.right = element_text(color = "black"),
#     axis.title.x = element_blank(),
#     axis.text.x = element_blank(),
#     axis.title.y.left = element_blank(),  
#     axis.ticks.y.left = element_blank(),  
#     axis.ticks.x = element_blank(),  
#     axis.text.y.left = element_blank(),  
#     axis.title.y.right = element_text(size = 16, angle = 270, color = "black")
#   )
# 
# print(heatmap_plot_W)
# 
# combined_plot_all_perc <- heatmap_plot_D /  heatmap_plot_W / 
#   VPD_plot_1 / VPD_plot_2 +
#   plot_layout(ncol = 2, nrow = 2, heights = c(1, 0.15), widths = c(1, 1))
# combined_plot_all_perc

write.csv(heatmap_data_sorted, 
          "./data/GH/results/Irrigation_thresholds/Irrigation_thres_GH_2023.csv")

# find max. percentage_w_cros where each tree did still recover to an irrigation
# demand of 0
result_max_still_recovery <- heatmap_data_sorted %>%
  filter(Date >= "2023-06-15" & Date <= "2023-09-06",
         Treatment == "Droughted") %>%
  group_by(Tree.ID) %>%
  group_split() %>%
  lapply(find_threshold_change) %>%
  bind_rows()

