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
# - Combine TWD with waterpotential measurements (WP)

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
                        zoo::na.locf(E22_W), E22_W)) 

# use only data from 2023-06-21 10:00:00 onwards so less NAs exist and until
# 2023-09-07 00:00:00
# which will make the data cleaning and processing more easy
dendro <- dendro %>% filter(TIMESTAMP>="2023-06-21 10:00:00" &
                              TIMESTAMP<="2023-09-07 00:00:00")


# Eliminate obvious measurement errors and artefacts
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
# and treatment with network interpolation from 
# DendRoAnalyst (Aryal et al. 2020), we need to divide 
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

# select DendRoAnalyst network interpolated and jump-corrected data, if no
# corrections were performed with dendRoAnalyst, keep the original data
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
          file = "./data/processed/GH_TSD_start_diam_stats_new.csv")

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

# create empty dataframe to store results in
ZG_phase_DF <- data.frame()

for(i in 2:length(dendro_cor)){
  calc_plot_DA(i, dendro_cor, 2023, 
               "./graphics/Dendrometer_Processing/TWD/DA-TWD_Singletree/") }

write.csv(ZG_phase_DF, 
          "./data/processed/ZG_phase_DF_new.csv")

# change tree species and treatment labels
ZG_phase_DF$Tree.Species <- factor(ZG_phase_DF$Tree.Species,
                                   levels = c("B", "E", "D"),
                                   labels = c("B" = "Fagus sylvatica",
                                              "E" = "Quercus robur",
                                              "D" = "Pseudotsuga menziesii"))

ZG_phase_DF$Treatment <- factor(ZG_phase_DF$Treatment, 
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

TWD_rel_df$Date <- as.Date(TWD_rel_df$TIME, tz = "Europe/Berlin")

TWD_rel_df  <- TWD_rel_df  %>%
  group_by(Tree.ID, Date) %>%
  mutate(min_rel_TWD = min(rel_TWD_log, na.rm = TRUE)) 
  
write.csv(TWD_rel_df, "./data/processed/TWD_rel_df.csv")

Min_TWD_df <- TWD_rel_df %>%
  dplyr::select(Tree.ID, Date, min_rel_TWD) %>%
  unique()

write.csv(Min_TWD_df, "./data/processed/Min_rel_TWD_df.csv")
