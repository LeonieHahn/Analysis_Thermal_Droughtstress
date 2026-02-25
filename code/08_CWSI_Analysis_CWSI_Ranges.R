# - calculate CWSI per date and Tree.ID using the median values after the 
#   different percentile clipping versions
# - visualize CWSI per tree species, treatment and date
# - combine CWSI values with environmental and ecophysiological data
# - calculate time differences between CWSI acquisition, min., max. and rel. TWD
# - create overview with timestamps of thermal image recordings and 
# - get earliest and latest image time per date
# - add solar noon per date
# - explore how CWSI is affected by image acquisition time / 
#   distance tosolar noon 



library(readxl)
library(stringr)
library(readr)
library(tidyr)
library(purrr)
library(ggplot2)
library(lubridate)
library(scales)
library(suncalc)
library(ppcor)
library(dplyr)


segment_info <- read_csv("./data/thermal/TreeSpecies_Dates_TreeIDs.csv")
segment_info$Date <- as.Date(as.character(segment_info$Date), format = "%Y%m%d")
segment_info <- segment_info %>%
  mutate(Tree.ID = str_extract(Tree_ID, "[A-Z][0-9]+_[DW]"))

temp_stats_df <- read_csv("./data/thermal/results/percentile_and_allpixels_stats.csv")

merged_stats <- left_join(temp_stats_df[, -1], segment_info[, -c(9, 11)],
                          by = c("Tree.ID", "Date"))

merged_stats <- merged_stats %>%
  filter(Action != "Exclude" | is.na(Action))

# keep only those observations where data for the same position exists for
# every object. In case there is data available for both conditions 
# (sunny/shady), use what is listed in column "Action"
complete_per_position <- merged_stats %>%
  group_by(Date, Tree.ID, Position) %>%
  summarise(n_objects = n_distinct(object), .groups = "drop") %>%
  filter(n_objects == 3)

complete_with_action <- complete_per_position %>%
  left_join(
    merged_stats %>%
      dplyr::select(Date, Tree.ID, Action) %>%
      distinct(),
    by = c("Date", "Tree.ID")
  )

preferred_position <- complete_with_action %>%
  group_by(Date, Tree.ID) %>%
  mutate(
    keep = case_when(
      Action == "use Tree_Shade" & Position == "Shade" ~ TRUE,
      Action == "use Tree_Sun"   & Position == "Sun"   ~ TRUE,
      is.na(Action) ~ TRUE,  # no action --> allow both
      TRUE ~ FALSE
    )
  ) %>%
  filter(keep) %>%
  dplyr::select(Date, Tree.ID, Position) %>%
  ungroup()

cleaned_df <- merged_stats %>%
  inner_join(preferred_position, by = c("Date", "Tree.ID", "Position"))

cleaned_df <- cleaned_df[, 1:16] %>%
  filter(Tree.Species != "Abies alba") # remove A. alba from analysis
names(cleaned_df)[8:15] <- paste0("temp_", names(cleaned_df)[8:15])

# check for duplicates
duplicates_position <- cleaned_df %>%
  distinct(Date, Tree.ID, Position) %>%
  group_by(Date, Tree.ID) %>%
  filter(n() > 1) %>%
  arrange(Date, Tree.ID)

duplicates_position

# calculate CWSI for median_objectarea

# but: keep rangetype (from quantile clipping)
CWSI_df <- cleaned_df %>%
  pivot_wider(
    id_cols = c(Date, Tree.ID, Treatment, Tree.Species, RangeType),
    names_from = object,
    values_from = c(
      temp_mean_objectarea, temp_median_objectarea, temp_min_objectarea,
      temp_max_objectarea, temp_sd_objectarea, Position, temp_pixelamount_used,
      temp_mean_median_dif, temp_range
    ),
    names_sep = "."
  ) %>%
  mutate(
    consistent = Position.Tree == Position.Green_Dry & 
      Position.Tree == Position.Green_Wet,
    CWSI = ifelse(
      consistent,
      (temp_median_objectarea.Tree - temp_median_objectarea.Green_Wet) /
        (temp_median_objectarea.Green_Dry - temp_median_objectarea.Green_Wet),
      NA_real_
    )
  ) %>%
  dplyr::select(
    Date,
    RangeType,  # <- keep rangetype
    Position = Position.Tree,
    Tree.ID,
    Treatment,
    Tree.Species,
    temp_mean_objectarea.Tree,
    temp_median_objectarea.Tree,
    temp_min_objectarea.Tree,
    temp_max_objectarea.Tree,
    temp_sd_objectarea.Tree,
    temp_pixelamount_used.Tree,
    temp_mean_median_dif.Tree,
    temp_range.Tree,
    CWSI
  )


CWSI_df <- left_join(CWSI_df, cleaned_df[, c(1, 3, 4, 7, 16)],
                     by = c("Date", "Position", "Tree.ID", "RangeType")) %>%
  mutate(DateTime = str_c(str_extract(Tree_ID, "\\d{4}-\\d{2}-\\d{2}"),  
                          str_extract(Tree_ID, "(?<=_)\\d{2}-\\d{2}-\\d{2}(?=\\+)"), 
                          sep = " ") ) %>%
  mutate(DateTime_POSIX = as.POSIXct(DateTime, format = "%Y-%m-%d %H-%M-%S",
                                     tz = "Europe/Berlin"),
         DateTime_rounded = round_date(DateTime_POSIX, unit = "10 minutes")) %>%
  dplyr::select(-DateTime_POSIX) %>%
  mutate(
    Treatment = case_when(
      Treatment == "W" ~ "Watered",
      Treatment == "D" ~ "Droughted",
      TRUE ~ Treatment
    )
  )

CWSI_df <- CWSI_df%>%
  distinct()

# read in climate data for adding RH, Temp and VPD at thermal measurement
climate <- read_excel(
  "./data/GH2023_Climate_combined.xlsx",
  sheet = 1) 

climate <- climate[, c(1, 24, 21, 26)] # version with Tair_min corrected VPD+RH

# climate <- climate[, c(1, 17, 16, 18)] # version with uncorrected VPD

# climate <- climate[, c(1, 36, 33, 38)] # version with Tair_median corrected VPD+RH

climate$DateTime <- force_tz(climate$DateTime, tzone = "Europe/Berlin")

colnames(climate) <- c("DateTime", "RH", "Air_temp", "VPD_kPa")

# calculate daily mean, min and max VPD
climate <- climate %>%
  mutate(DateTime = force_tz(DateTime, tzone = "Europe/Berlin")) %>%
  mutate(Date = as_date(DateTime)) %>%
  group_by(Date) %>%
  mutate(
    VPD_adj_dailyavg_kPa = mean(VPD_kPa, na.rm = TRUE),
    VPD_adj_dailymax_kPa = max(VPD_kPa, na.rm = TRUE),
    VPD_adj_dailymin_kPa = min(VPD_kPa, na.rm = TRUE)
  ) %>%
  ungroup() %>%
  dplyr::select(-Date)

# read in irrigation demand, twd, environmental data
irrig <- read_csv("./data/processed/Irrigation_thres_GH_2023.csv",
                  col_select = -c(1, 6, 7, 8))

# read in min. and max. TWD and WP data (dendrometer trees)
TWD_WP <- read_csv("./data/processed/TWD_WP_combined.csv",
                   col_select = -1)

# read in WP measurements of all trees 
wp <- read_csv("./data/Scholander_WP_data.csv",
               col_select = -c(1, 2))

colnames(wp)[1]<- "Date"
wp$Date <- as.Date(strptime(wp$Date, format = "%d.%m.%Y"))
wp$Tree.ID <- paste0(wp$Tree.ID, "_", wp$Treatment)
wp <- pivot_wider(wp[, c(1, 2, 5, 22)], names_from = pd.md, 
                  names_prefix = "WP_", values_from = WP_average..MPa.)

# read in TWD from all dendrometers in 10 minute resolution
TWD_10 <- read_csv("./data/processed/ZG_phase_rel_logTWD.csv",
                   col_select = -1)
TWD_10$TIME <- sub(" UTC$", " CEST", TWD_10$TIME)
# fill in "00:00:00" for the midnight hours, since they are missing
TWD_10$TIME <- ifelse(nchar(TWD_10$TIME) == 10,
                      paste(TWD_10$TIME, "00:00:00"),
                      TWD_10$TIME)
TWD_10$TIME <- as.POSIXct(TWD_10$TIME, format="%Y-%m-%d %H:%M:%S", 
                          tz = "Europe/Berlin")

# create info that it is a dendrometer tree
dendro_trees <- data.frame(Tree.ID = unique(TWD_10$Tree.ID)) %>%
  mutate(Dendro_info = "DendroTree")

# read in soil moisture of all trees
swc <- read_csv(
  "./data/GH_VWC_TDR_means_cleaned+daily_interp_womaxVWC_alltrees.csv",
  col_select = -1)

swc$Tree.ID <- paste0(swc$Tree.ID, "_", swc$Treatment)

# join data to CWSI dataframe
irrig_CWSI <- full_join(CWSI_df, irrig,
                        by = c("Tree.ID", "Treatment", "Tree.Species", "Date"))

irrig_CWSI <- full_join(irrig_CWSI, TWD_WP[, 1:4], 
                        by = c("Tree.ID", "Date"))

irrig_CWSI <- left_join(irrig_CWSI, wp, by = c("Date", "Tree.ID"))


irrig_CWSI <- left_join(
  irrig_CWSI,
  TWD_10,
  by = c("DateTime_rounded" = "TIME", "Tree.ID", "Treatment", "Tree.Species")
)

irrig_CWSI$VWC_mean <- NULL
irrig_CWSI <- full_join(irrig_CWSI, swc[, 1:3], 
                        by = c("Date", "Tree.ID"))

irrig_CWSI <- left_join(irrig_CWSI, dendro_trees, by = "Tree.ID")
irrig_CWSI <- left_join(irrig_CWSI, climate, 
                        by = c("DateTime_rounded" = "DateTime"))

irrig_CWSI <- irrig_CWSI %>%
  dplyr::select(-c(Tree_ID, DateTime, min_rel_TWD, FirstDayAfterLastGRO, 
            Percentage_color, Phases, ThresholdDate, Tree.ID_num, 
            Tree.ID_num_cat)) %>%
  mutate(Dendro_info = ifelse(grepl("DendroTree", Dendro_info), "DendroTree",
                              "NonDendroTree"),
         Treatment = ifelse(grepl("_D$", Tree.ID), "D",
                            ifelse(grepl("_W$", Tree.ID), "W", NA)),
         Tree.Species = ifelse(grepl("^E", Tree.ID), "Q. robur",
                               ifelse(grepl("^D", Tree.ID), "P. menziesii",
                                      ifelse(grepl("^T", Tree.ID), "A. alba",
                                             ifelse(grepl("^B", Tree.ID), "F. sylvatica", NA)))))


# join all 10 min data with CWSI
combined_df <- left_join(TWD_10, climate, by = c("TIME" = "DateTime"))
combined_df$Date <- as.POSIXct(
  format(combined_df$TIME, "%Y-%m-%d"),  
  tz = "Europe/Berlin"
)

combined_df <- left_join(combined_df, swc[, 1:3], 
                         by = c("Date", "Tree.ID"))


CWSI_df_joinsel <- CWSI_df %>%
  filter(RangeType == "Q10_Q90")

combined_df <- left_join(combined_df, TWD_WP[, 1:6], by = c("Date", "Tree.ID"))
combined_df <- left_join(combined_df, CWSI_df_joinsel[, c(1, 4, 15)],
                         by = c("Date", "Tree.ID"))

# exclude CWSI values below 0 and above 1 since they are unrealistic ->
# measurement errors

# count how often CWSI values are below 0 or above 1
count_CWSI_unreal <-  irrig_CWSI %>%
  filter(CWSI < 0 | CWSI > 1) %>%
  filter(RangeType == "Q10_Q90")

# count how many CWSI values there are in total
count_CWSI <-  irrig_CWSI %>%
  filter(RangeType == "Q10_Q90")

# cout per tree species
count_CWSI_unreal_spec <- count_CWSI_unreal %>%
  group_by(Tree.Species) %>%
  count()

irrig_CWSI <- irrig_CWSI %>%
  mutate(CWSI = ifelse(CWSI < 0 | CWSI > 1, "NA", CWSI)) %>%
  mutate(CWSI = as.numeric(CWSI))

combined_df <- combined_df%>%
  mutate(CWSI = ifelse(CWSI < 0 | CWSI > 1, "NA", CWSI)) %>%
  mutate(CWSI = as.numeric(CWSI))

# remove A. alba
combined_df <- combined_df %>%
  filter(Tree.Species != "Abies alba")
sum(duplicated(combined_df))

write.csv(combined_df, "./data/results/Env_Eco_CWSI_combined_10min.csv",
          row.names = FALSE)

# write.csv(combined_df, "./data/results/Env_Eco_CWSI_combined_10min_VPD_uncor.csv",
#           row.names = FALSE)
# 
# write.csv(combined_df, "./data/results/Env_Eco_CWSI_combined_10min_VPD_median_cor.csv",
#           row.names = FALSE)

# combine CWSI_df and TWD info to calculate the time differences between
# CWSI aqcuisition and min., max. and rel. TWD

# calculate time of min. and max. TWD

TWD_10_min_max_df <-  TWD_10 %>%
  filter(Tree.Species != "Abies alba") %>%

  # change timezone
  mutate(
    TIME = with_tz(TIME, tzone = "Europe/Berlin"),
    Date = as_date(TIME),
    # extract time
    Time_only = format(TIME, "%H:%M:%S")
  ) %>%
  group_by(Tree.ID, Date) %>%
  summarize(
    # select daily min. TWD before noon
    min_rel_TWD_log = min(rel_TWD_log[Time_only < "12:00:00"], na.rm = TRUE),
    
    # select latest time at this morning
    min_time = max(
      TIME[
        rel_TWD_log == min(rel_TWD_log[Time_only < "12:00:00"], na.rm = TRUE) &
          Time_only < "12:00:00"
      ]
    ),
    
    # select daily max. TWD (earliest time)
    max_rel_TWD_log = max(rel_TWD_log, na.rm = TRUE),
    max_time = TIME[which.max(rel_TWD_log)],
    
    .groups = "drop"
  )

CWSI_TWD_timecompare <- left_join(TWD_10_min_max_df, 
                                  CWSI_df[, c(1, 2, 4:6, 15, 17, 18)] %>%
                                    filter(RangeType == "Q10_Q90") %>%
                                     mutate(CWSI = ifelse(CWSI < 0 | CWSI > 1, 
                                                          "NA", CWSI)) %>%
                                             mutate(CWSI = as.numeric(CWSI)), 
                                  by = c("Tree.ID", "Date")) %>%
  dplyr::select(-c("RangeType", "CWSI")) 

colnames(CWSI_TWD_timecompare)[9] <- "CWSI_acq_time"


CWSI_TWD_timecompare <- CWSI_TWD_timecompare %>%
  
  # Converte CWSI_acq_time in POSIXct (Europe/Berlin)
  mutate(
    CWSI_acq_time_clean = gsub("-", ":", CWSI_acq_time),      
    CWSI_acq_time_clean = sub(":([0-9]{2}):([0-9]{2})$", "-\\1-\\2", 
                              CWSI_acq_time_clean), 
  ) %>%
  
  mutate(
    CWSI_acq_time_POSIX = parse_date_time(
      CWSI_acq_time_clean,
      orders = "Ymd HMS",
      tz = "Europe/Berlin"
    )
  ) %>%
  
  mutate(
    DateTime_rounded_POSIX = with_tz(DateTime_rounded, "Europe/Berlin"),
    
    # Calculate time difference between CWSI aqcuisition and rounded time 
    # (= time of rel.TWD used) in seconds
    rel_TWD_time_diff_sec = as.numeric(
      difftime(DateTime_rounded_POSIX, CWSI_acq_time_POSIX, units = "secs")
    ),
    # Calculate time difference between CWSI aqcuisition and min. TWD time
    min_TWD_time_diff_sec = as.numeric(
      difftime(min_time, CWSI_acq_time_POSIX, units = "secs")
    ),
    # Calculate time difference between CWSI aqcuisition and max. TWD time
    max_TWD_time_diff_sec = as.numeric(
      difftime(max_time, CWSI_acq_time_POSIX, units = "secs")
    )
  )

# calculate average and min-max ranges for the time differences
CWSI_TWD_timecompare_stats <- CWSI_TWD_timecompare %>%
  summarize(mean_rel_TWD_time_diff_sec = mean(rel_TWD_time_diff_sec, na.rm = TRUE),
            min_rel_TWD_time_diff_sec = min(rel_TWD_time_diff_sec, na.rm = TRUE),
            max_rel_TWD_time_diff_sec = max(rel_TWD_time_diff_sec, na.rm = TRUE),
            mean_min_TWD_time_diff_sec = mean(min_TWD_time_diff_sec, na.rm = TRUE),
            min_min_TWD_time_diff_sec = min(min_TWD_time_diff_sec, na.rm = TRUE),
            max_min_TWD_time_diff_sec = max(min_TWD_time_diff_sec, na.rm = TRUE),
            mean_max_TWD_time_diff_sec = mean(max_TWD_time_diff_sec, na.rm = TRUE),
            min_max_TWD_time_diff_sec = min(max_TWD_time_diff_sec, na.rm = TRUE),
            max_max_TWD_time_diff_sec = max(max_TWD_time_diff_sec, na.rm = TRUE))

# recalculate values from seconds into minutes and hours
CWSI_TWD_timecompare_stats <- CWSI_TWD_timecompare_stats %>%
  mutate(
    across(ends_with("_sec"), 
           list(
             min  = ~ .x / 60,       
             hour = ~ .x / 3600       
           ),
           .names = "{.col}_{.fn}")
  )


irrig_CWSI <- irrig_CWSI %>%
  filter(!is.na(CWSI)) %>%
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),  
    Date_label = factor(Date_label, levels = unique(Date_label))  
  )

# fill VPD data for all trees (so far only for dendro trees, but this data is
# the same for all trees)           
vpd_per_day <- irrig_CWSI %>%
  filter(!is.na(VPD_adj_dailyavg_kPa)) %>%
  group_by(Date) %>%
  summarise(
    VPD_adj_dailyavg_kPa = first(VPD_adj_dailyavg_kPa),
    VPD_adj_dailymax_kPa = first(VPD_adj_dailymax_kPa),
    VPD_adj_dailymin_kPa = first(VPD_adj_dailymin_kPa),
    .groups = "drop"
  )

irrig_CWSI <- irrig_CWSI %>%
  dplyr::select(-VPD_adj_dailyavg_kPa, -VPD_adj_dailymax_kPa, -VPD_adj_dailymin_kPa) %>%
  left_join(vpd_per_day, by = "Date")

irrig_CWSI_allranges <- irrig_CWSI %>%
  rename(Mean_VPD_kPa = VPD_adj_dailyavg_kPa,
         Min_VPD_kPa = VPD_adj_dailymin_kPa,
         Max_VPD_kPa = VPD_adj_dailymax_kPa,
         Irrigation_demand = Percentage_w_cros,
         SWC = VWC_mean)

write.csv(irrig_CWSI_allranges, "./data/results/Irrig_CWSI_allranges.csv", 
          row.names = FALSE)

# write.csv(irrig_CWSI_allranges, "./data/results/Irrig_CWSI_allranges_VPD_uncor.csv",
#           row.names = FALSE)
# 
# write.csv(irrig_CWSI_allranges, "./data/results/Irrig_CWSI_allranges_VPD_median_cor.csv",
#           row.names = FALSE)


# compare CWSI values when using different percentiles for outlier thresholding

ggplot(irrig_CWSI_allranges, 
       aes(x = Date, y = CWSI, color = RangeType, group = RangeType)) +
  stat_summary(fun = mean, geom = "line", linewidth = 1) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2) +
  theme_minimal() +
  facet_grid(Tree.Species ~ Treatment) +
  labs(title = "Mean CWSIper Rangetype and error margins")


# Calculate differences between CWSI after percentile clipping and CWSI 
# using all values in long format

cwsi_diffs_long <- irrig_CWSI_allranges %>%
  dplyr::select(Date, Tree.Species, Treatment, Tree.ID, RangeType, CWSI) %>%

  pivot_wider(
    names_from = RangeType,
    values_from = CWSI,
    values_fn = mean  
  ) %>%
  
  pivot_longer(
    cols = c(Q10_Q90, Q5_Q95, Q1_Q99, Q15_Q85, IQR_based),
    names_to = "RangeType",
    values_to = "CWSI"
  ) %>%
  
  mutate(diff_to_allpix = CWSI - All_Pixels)


# Boxplot for the CWSI differences calculated from different median values 
# before and after the percentile clipping

# Plot for S1
ggplot(cwsi_diffs_long %>%
         mutate( Date_label = format(Date, "%d.%m."),
                 Date_label = factor(Date_label, 
                                     levels = unique(format(sort(unique(Date)), 
                                                            "%d.%m."))),
           RangeType = dplyr::recode(RangeType,
                                   "Q1_Q99"    = "1.–99.% Quantile",
                                   "Q5_Q95"    = "5.–95.% Quantile",
                                   "Q10_Q90"   = "10.–90.% Quantile",
                                   "Q15_Q85"   = "15.–85.% Quantile",
                                   "IQR_based" = "IQR±1.5×IQR"
         ), 
       RangeType = factor(
         RangeType,
         levels = c(
           "IQR±1.5×IQR",
           "1.–99.% Quantile",
           "5.–95.% Quantile",
           "10.–90.% Quantile",
           "15.–85.% Quantile"
         ))),
       aes(x = RangeType, y = diff_to_allpix, fill = RangeType)) +
  geom_boxplot(alpha = 0.7) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
  facet_grid(Tree.Species ~ Date_label ) +
  theme_minimal() +
  labs(
    y = expression("CWSI difference (" * CWSI[cleaned] - 
                     CWSI[all] * ")"),
    x = "Quantile ranges"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(size = 14, angle = 90, hjust = 1, vjust = 0.5),
    axis.text.y = element_text(size = 14),
    axis.title.x = element_text(size = 16,  margin = margin(t = 10)),
    axis.title.y = element_text(size = 16,  margin = margin(r = 10)),
    strip.text = element_text(size = 14, face = "bold"),  
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    panel.grid.minor = element_blank(),
    panel.grid.major = element_line(linewidth = 0.3, color = "grey80")
  )

ggsave("./graphics/results/Overview_CWSI_difs_percentils_all_pixels_Date_Species_BED.png", 
       width = 36, 
       height = 25,
       units = "cm",
       dpi = 600)

irrig_CWSI <- irrig_CWSI_allranges %>%
  filter(RangeType == "Q10_Q90") 
irrig_CWSI$RangeType <- NULL

write.csv(irrig_CWSI, "./data/results/Irrig_CWSI.csv", row.names = FALSE)
# write.csv(irrig_CWSI, "./data/results/Irrig_CWSI_VPD_uncor.csv", row.names = FALSE)
# write.csv(irrig_CWSI, "./data/results/Irrig_CWSI_VPD_median_cor.csv", row.names = FALSE)


# create overview plots
# median leaf temperatures all trees
pvals_median_temp <- irrig_CWSI %>%
  filter(!is.na(temp_median_objectarea.Tree)) %>%
  group_by(Tree.Species, Date, Position) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(temp_median_objectarea.Tree ~ Treatment)$p.value, 
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)     ~ "",
      p_value < 0.001    ~ "***",
      p_value < 0.01     ~ "**",
      p_value < 0.05     ~ "*",
      TRUE               ~ ""  
    )
  ) %>% 
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),
    Date_label = factor(Date_label, levels = unique(Date_label))
  )


# Plot with median temperatures
ggplot(irrig_CWSI %>%
         filter(!is.na(temp_median_objectarea.Tree))) +
  geom_boxplot(aes(x = Treatment, y = temp_median_objectarea.Tree, 
                   color = Treatment)) +
  facet_grid(Tree.Species ~ Position ~ Date_label) +
  scale_color_manual(values = c("D" = "red", "W" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  ylab("Median leaf temperature [°C]") +
  # add significance level
  geom_text(
    data = pvals_median_temp,
    aes(
      x = 1.5,  
      y = Inf, 
      label = signif_label
    ),
    vjust = 1.5, inherit.aes = FALSE, size = 5
  )

ggsave("./graphics/results/Overview_median_leaf_temp_BED.png", 
       width = 36, 
       height = 25,
       units = "cm",
       dpi = 600)

# median leaf temperatures DendroTrees
pvals_median_temp_dendro <- irrig_CWSI %>%
  filter(Dendro_info == "DendroTree" & !is.na(temp_median_objectarea.Tree)) %>%
  group_by(Tree.Species, Date, Position) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(temp_median_objectarea.Tree ~ Treatment)$p.value, 
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)     ~ "",
      p_value < 0.001    ~ "***",
      p_value < 0.01     ~ "**",
      p_value < 0.05     ~ "*",
      TRUE               ~ ""  
    )
  ) %>% 
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),
    Date_label = factor(Date_label, levels = unique(Date_label))
  )


# Plot with median temperatures
ggplot(irrig_CWSI %>%
         filter(Dendro_info == "DendroTree" & 
                  !is.na(temp_median_objectarea.Tree))) +
  geom_boxplot(aes(x = Treatment, y = temp_median_objectarea.Tree, 
                   color = Treatment)) +
  facet_grid(Tree.Species ~ Position ~ Date_label) +
  scale_color_manual(values = c("D" = "red", "W" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  ylab("Median leaf temperature [°C]") +
  # add significance level
  geom_text(
    data = pvals_median_temp_dendro,
    aes(
      x = 1.5,  
      y = Inf, 
      label = signif_label
    ),
    vjust = 1.5, inherit.aes = FALSE, size = 5
  )

ggsave("./graphics/results/Overview_median_leaf_temp_DendroTrees_BED.png", 
       width = 36, 
       height = 25,
       units = "cm",
       dpi = 600)               

# Mean leaf temperatures
pvals_mean_temp <- irrig_CWSI %>%
  filter(!is.na(temp_mean_objectarea.Tree)) %>%
  group_by(Tree.Species, Date, Position) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(CWSI ~ Treatment)$p.value,
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)     ~ "",
      p_value < 0.001    ~ "***",
      p_value < 0.01     ~ "**",
      p_value < 0.05     ~ "*",
      TRUE               ~ ""  
    )
  ) %>% 
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),
    Date_label = factor(Date_label, levels = unique(Date_label))
  )


ggplot(irrig_CWSI %>%
         filter(Dendro_info == "DendroTree" & 
                  !is.na(temp_mean_objectarea.Tree))) +
  geom_boxplot(aes(x = Treatment, y = temp_mean_objectarea.Tree, 
                   color = Treatment)) +
  facet_grid(Tree.Species ~ Position ~ Date_label) +
  scale_color_manual(values = c("D" = "red", "W" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  ylab("Mean leaf temperature [°C]") +
  # add significance level
  geom_text(
    data = pvals_mean_temp,
    aes(
      x = 1.5,  
      y = Inf, 
      label = signif_label
    ),
    vjust = 1.5, inherit.aes = FALSE, size = 5
  )

ggsave("./graphics/results/Overview_mean_leaf_temp_DendroTrees.png", 
       width = 36, 
       height = 25,
       units = "cm",
       dpi = 600)

# CWSI

# Plot with CWSI all trees
pvals_CWSI <- irrig_CWSI %>%
  filter(!is.na(CWSI)) %>%
  group_by(Tree.Species, Date) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(CWSI ~ Treatment)$p.value, 
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)     ~ "",
      p_value < 0.001    ~ "***",
      p_value < 0.01     ~ "**",
      p_value < 0.05     ~ "*",
      TRUE               ~ ""  
    )
  ) %>% 
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),
    Date_label = factor(Date_label, levels = unique(Date_label))
  )


# Plot with CWSI all trees
ggplot(irrig_CWSI %>%
         filter(!is.na(CWSI))) +
  geom_boxplot(aes(x = Treatment, y = CWSI, 
                   color = Treatment)) +
  facet_grid(Tree.Species ~ Date_label) +
  scale_color_manual(values = c("D" = "red", "W" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  ylab("CWSI") +
  # add significance level
  geom_text(
    data = pvals_CWSI,
    aes(
      x = 1.5,  
      y = Inf, 
      label = signif_label
    ),
    vjust = 1.5, inherit.aes = FALSE, size = 5
  )

ggsave("./graphics/results/Overview_CWSI_BED.png", 
       width = 36, 
       height = 20,
       units = "cm",
       dpi = 600)

# CWSI DendroTrees
pvals_CWSI_dendro <- irrig_CWSI %>%
  filter(Dendro_info == "DendroTree" & !is.na(CWSI)) %>%
  group_by(Tree.Species, Date) %>%
  summarise(
    p_value = tryCatch(
      wilcox.test(CWSI ~ Treatment)$p.value, 
      error = function(e) NA
    ),
    .groups = "drop"
  ) %>%
  mutate(
    signif_label = case_when(
      is.na(p_value)     ~ "",
      p_value < 0.001    ~ "***",
      p_value < 0.01     ~ "**",
      p_value < 0.05     ~ "*",
      TRUE               ~ ""  
    )
  ) %>% 
  arrange(Date) %>%
  mutate(
    Date_label = format(Date, "%d.%m."),
    Date_label = factor(Date_label, levels = unique(Date_label))
  )


# Plot with CWSI for DendroTrees
ggplot(irrig_CWSI %>%
         filter(Dendro_info == "DendroTree" & 
                  !is.na(CWSI))) +
  geom_boxplot(aes(x = Treatment, y = CWSI, 
                   color = Treatment)) +
  facet_grid(Tree.Species ~ Date_label) +
  scale_color_manual(values = c("D" = "red", "W" = "blue")) +
  theme_bw() +
  theme(axis.text.x = element_text(hjust = 1),
        strip.text = element_text(size = 14), 
        axis.text = element_text(size = 12),
        axis.title = element_text(size = 16),
        legend.title = element_text(size = 16),
        legend.text = element_text(size = 12)) +
  ylab("CWSI") +
  # add significance level
  geom_text(
    data = pvals_CWSI_dendro,
    aes(
      x = 1.5,  
      y = Inf,
      label = signif_label
    ),
    vjust = 1.5, inherit.aes = FALSE,
    size = 5
  )

ggsave("./graphics/results/Overview_CWSI_DendroTrees_BED.png", 
       width = 36, 
       height = 20,
       units = "cm",
       dpi = 600)    


# create overview with timestamps of thermal image recordings and 

temp_files <- lapply(list.files(
  "./data/Thermal/tempvalue_tables/", full.names = TRUE), 
  function(f) read_csv(f, col_select = -1))

temp_df <- Reduce(rbind, temp_files)
temp_df <- temp_df %>%
  # exclude Abies alba as it is not used in this study
  filter(Tree.Species != "Abies Alba") %>%
  mutate(
    date_str = str_extract(imgname, "\\d{4}-\\d{2}-\\d{2}"),
    time_str = str_extract(imgname, "(?<=_)\\d{2}-\\d{2}-\\d{2}") %>%
      str_replace_all("-", ":"),
    datetime = ymd_hms(paste(date_str, time_str), tz = "Europe/Berlin"),
    Time = format(datetime, "%H:%M:%S")
  )

time_range_per_day <- temp_df %>%
  group_by(Date) %>%
  summarise(
    earliest_time = format(min(datetime, na.rm = TRUE), "%H:%M:%S"),
    latest_time   = format(max(datetime, na.rm = TRUE), "%H:%M:%S"),
    .groups = "drop"
  )

# get earliest and latest image time per date
time_range_per_day <- temp_df %>%
  group_by(Date) %>%
  summarise(
    earliest_time = min(Time, na.rm = TRUE),
    latest_time   = max(Time, na.rm = TRUE),
    .groups = "drop"
  )

# add solar noon per date
solar_noon_df <- getSunlightTimes(
  date = time_range_per_day$Date,
  lat  = 48.399400,    
  lon  = 11.717798,    
  tz   = "Europe/Berlin",
  keep = "solarNoon"
)

time_range_with_noon <- time_range_per_day %>%
  bind_cols(solar_noon_df %>% dplyr::select(solarNoon)) %>%
  mutate(
    solarNoonTime = format(solarNoon, "%H:%M:%S") 
  )

# calculate time differences between solar noon and earliest and lates image
time_range_with_noon <- time_range_with_noon %>%
  mutate(
    earliest_POSIX = as.POSIXct(paste(Date, earliest_time), 
                                format = "%Y-%m-%d %H:%M:%S", 
                                tz = "Europe/Berlin"),
    latest_POSIX   = as.POSIXct(paste(Date, latest_time),   
                                format = "%Y-%m-%d %H:%M:%S", 
                                tz = "Europe/Berlin"),
    solarNoon_POSIX = solarNoon,  
    Diff_Earliest  = as.numeric(difftime(solarNoon_POSIX, 
                                         earliest_POSIX, units = "mins")),
    Diff_Latest    = as.numeric(difftime(latest_POSIX, 
                                         solarNoon_POSIX, units = "mins")),
    Diff_Earliest_h  = as.numeric(difftime(solarNoon_POSIX, 
                                           earliest_POSIX, units = "hours")),
    Diff_Latest_h    = as.numeric(difftime(latest_POSIX, 
                                           solarNoon_POSIX, units = "hours"))
  ) %>%
  dplyr::select(Date, earliest_time, latest_time, 
         solarNoonTime, Diff_Earliest, Diff_Latest, 
         Diff_Earliest_h, Diff_Latest_h) 


ggplot(temp_df %>%
         filter(Position  %in% 
                  c("Tree_Shade", "Green_Dry_Shade", "Green_Wet_Shade",
                    "Tree_Sun", "Green_Dry_Sun", "Green_Wet_Sun",
                    "Geen_Dry_Shade", "GReen_Dry_Shade", "Green_wet_Sun"))) +
  geom_point(aes(x =as.POSIXct(paste(Time), 
                               format = "%H:%M:%S", 
                               tz = "Europe/Berlin"), 
                 y = Temperature, color = Treatment)) +
  facet_grid(Position ~ Tree.Species) +
  xlab("Time") +
  theme_bw()

# create table with image acquisition times and time difference to solar noon
imagetimes_solarnoon <- left_join(temp_df, solar_noon_df, 
                                  by = c("Date" = "date"))
imagetimes_solarnoon <- imagetimes_solarnoon %>%
  mutate(TimeDiff_SolarNoon_h  = as.numeric(
                                difftime(solarNoon, datetime, 
                                         units = "hours")),
         Tree.ID = paste0(Tree.ID, "_", Treatment)) %>%
  dplyr::select("Date", "Tree.ID", "Treatment", "Time", "solarNoon",
         "TimeDiff_SolarNoon_h") %>%
  unique()

# compare how CWSI is affected by image acquisition time and time difference to
# solar noon
cwsi_solarnoon <- left_join(irrig_CWSI, imagetimes_solarnoon,
                            by = c("Tree.ID", "Treatment", "Date"))


# calculate R² for the relationship CWSI and TimeDiff_SolarNoon_min
r2_df <- cwsi_solarnoon %>%
  summarise(
    model = list(lm(CWSI ~ TimeDiff_SolarNoon_h, 
                    data = pick(everything())))
  ) %>%
  mutate(
    R2 = summary(model[[1]])$r.squared,
    p_value = summary(model[[1]])$coefficients[2, 4],
    signif_code = case_when(
      p_value < 0.001 ~ "***",
      p_value < 0.01  ~ "**",
      p_value < 0.05  ~ "*",
      TRUE ~ ""
    ),
    R2_label = paste0("R² = ", round(R2, 2), signif_code)
  ) %>%
  select(R2, p_value, R2_label)

# Calculate partial correlation between CWSI and TimeDiff_SolarNoon_min
# without the effect of VPD_kPa
df_test <- cwsi_solarnoon %>%
  dplyr::select(CWSI, TimeDiff_SolarNoon_h, VPD_kPa) %>%
  na.omit()

pcor_res <- pcor.test(
  df_test$CWSI,
  df_test$TimeDiff_SolarNoon_h,
  df_test$VPD_kPa,
  method = "spearman"
)

pcor_df <- tibble(
  pcor_r = pcor_res$estimate,
  p_value = pcor_res$p.value,
  label = paste0(
    "partial r = ", round(pcor_res$estimate, 2),
    ifelse(pcor_res$p.value < 0.001, " ***",
           ifelse(pcor_res$p.value < 0.01, " **",
                  ifelse(pcor_res$p.value < 0.05, " *", "")))
  )
)

pcor_df


# Plot with R² and partial correlation results
ggplot(cwsi_solarnoon, aes(x = TimeDiff_SolarNoon_h, y = CWSI, 
                           color = Tree.Species)) +
  geom_point(alpha = 0.6) +
  geom_smooth(aes(group = 1), method = "lm", se = TRUE, color = "black") +
 # geom_smooth(aes(color = Tree.Species), method = "lm", se = FALSE) +
  scale_color_manual(values = c(
    "F. sylvatica" = "#298c8c",
    "P. menziesii" = "orange",
    "Q. robur" = "darkmagenta"
  )) +
  
  # R² label 
  geom_text(
    data = r2_df,
    aes(
      x = 4,
      y = max(cwsi_solarnoon$CWSI, na.rm = TRUE),
      label = R2_label
    ),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    color = "black"
  ) +
  
  # Partial correlation label
  geom_text(
    data = pcor_df,
    aes(
      x = 4,
      y = max(cwsi_solarnoon$CWSI, na.rm = TRUE) - 0.05, 
      label = label
    ),
    inherit.aes = FALSE,
    hjust = 1, vjust = 1,
    color = "black"
  ) +
  
  labs(
    x = "Time difference from solar noon [hours]",
    y = "CWSI",
    color = "Tree species"
  ) +
  xlim(-3.3, 4) +
  theme_minimal()

ggsave("./graphics/results/correlations/correlation_CWSI_Time_SolarNoon.png")

# calculate how many cwsi measurements were taken between +- 2 hours around 
# solar noon
cwsi_near_noon <- cwsi_solarnoon %>%
  filter(TimeDiff_SolarNoon_h >= -2 & TimeDiff_SolarNoon_h <= 2)

n_values <- nrow(cwsi_near_noon)

perc_cwsi_near_noon <- n_values/nrow(cwsi_solarnoon)

