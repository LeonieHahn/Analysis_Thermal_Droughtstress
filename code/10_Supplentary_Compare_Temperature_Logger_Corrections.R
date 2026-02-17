# Compare the effects of the logger temperature correction for VPD calculation:
# - Plot Tair from the different sensors for a selected period of time
# - Compare correlation results between CWSI and the different Tairs, RHs 
#   and VPDs: using the Tmin corrected values, Tmedian instead of Tmin, and the
#   actual measurements from the sensor with potential overheating effects
# - Compare GAM performances when using different Tairs

    
library(ggplot2)
library(tidyr)
library(dplyr)
library(lubridate)
library(readxl)
library(stringr)
library(readr)

# read in climate data for adding RH, Temp and VPD at thermal measurement
climate <- read_excel(
  "./data/GH2023_Climate_combined.xlsx",
  sheet = 1,
  col_types = c("date", rep("text", 38))  
)

climate <- climate %>%
  mutate(
    DateTime = force_tz(DateTime, tzone = "Europe/Berlin"),
    across(
      -DateTime,
      ~ as.numeric(str_replace(as.character(.x), ",", "."))
    )
  )

climate <- climate %>%
  mutate(across(
    c(
      `T_Air1_LAT_Avg(GH1-1)_Deg C`,
      `T_Air2_LAT_Avg(GH1-3)_Deg C`,
      `T_Air3_LAT_Avg(GH2-1)_Deg C`,
      `T_Air4_LAT_Avg(GH2-3)_Deg C`
    ),
    as.numeric
  ))


# Plot with airtemperatures measured with different sensors

climate_long <- climate %>%
  select(DateTime,
         Temperature_Hobo1, Temperature_Hobo2,
         Temperature_Hobo3, Temperature_Hobo4,
         `T_Air1_LAT_Avg(GH1-1)_Deg C`,
         `T_Air2_LAT_Avg(GH1-3)_Deg C`,
         `T_Air3_LAT_Avg(GH2-1)_Deg C`,
         `T_Air4_LAT_Avg(GH2-3)_Deg C`,
         Temperature_Testo) %>%
  pivot_longer(-DateTime,
               names_to = "Sensor",
               values_to = "Tair")

ggplot(climate_long, aes(x = DateTime, y = Tair, color = Sensor)) +
  geom_line(alpha = 0.7) +
  xlim(as.POSIXct("2023-07-05"), as.POSIXct("2023-07-20")) +
  xlab("2023") +
  ylab("Tair [°C]") +
  theme_bw() +
  theme(legend.position = "bottom") +
  scale_color_manual(values = c(
    "Temperature_Hobo1" = "darkgreen",
    "Temperature_Hobo2" = "blue",
    "Temperature_Hobo3" = "chartreuse3",
    "Temperature_Hobo4" = "aquamarine1",
    "T_Air1_LAT_Avg(GH1-1)_Deg C" = "orange",
    "T_Air2_LAT_Avg(GH1-3)_Deg C" = "red",
    "T_Air3_LAT_Avg(GH2-1)_Deg C" = "magenta",
    "T_Air4_LAT_Avg(GH2-3)_Deg C" = "purple",
    "Temperature_Testo" = "black"
  ))


ggsave("graphics/results/Tair_sensor_differences.png", 
       width = 10, height = 6, dpi = 300)

# Compare Tmin, Tmedian and T from the potential overheated sensor
climate_difs <- climate %>%
  mutate(Dif_T_Tair_min = Temperature_Testo - Tair_min,
         Dif_T_Tair_median = Temperature_Testo - Tair_median)


ggplot(climate_difs) +
  geom_line(aes(x = DateTime, y = Dif_T_Tair_min)) +
  geom_line(aes(x = DateTime, y = Dif_T_Tair_median), 
             color = "red") +
  xlim(as.POSIXct("2023-06-20"), as.POSIXct("2023-08-20")) +
  xlab("2023") +
  ylab(expression(
    "Temperatur difference [°C]: " ~ T - T[min] ~ " / " ~ T - T[median])) +
  theme_bw()

# Compare VPD from Tmin, Tmedian and T from the potential overheated sensor
climate_difs <- climate %>%
  mutate(Dif_VPD_VPD_min = VPD_kPA - `VPD_adj (kPa)...26`,
         Dif_VPD_VPD_median = VPD_kPA - `VPD_adj (kPa)...38`)


ggplot(climate_difs) +
  geom_line(aes(x = DateTime, y = Dif_VPD_VPD_min, color = "VPD - VPDmin")) +
  geom_line(aes(x = DateTime, y = Dif_VPD_VPD_median, 
                color = "VPD - VPDmedian")) +
  xlim(as.POSIXct("2023-06-20"), as.POSIXct("2023-08-20")) +
  xlab("2023") +
  ylab(expression("VPD difference [kPa]")) +
  scale_color_manual(name = "", 
                     values = c("VPD - VPDmin" = "red", 
                                "VPD - VPDmedian" = "black")) +
  theme_bw() +
  theme(legend.position = "bottom")

ggsave("graphics/results/VPD_differences.png")


# compare correlations with CWSI
cor_min <- read_csv("data/results/statistical_analysis/Correlations_BDE_formatted.csv")
cor_median <- read_csv("data/results/statistical_analysis_VPD_median_cor/Correlations_BDE_formatted_VPD_median_cor.csv")
cor_orig <- read_csv("data/results/statistical_analysis_VPD_uncor/Correlations_BDE_formatted_VPD_uncor.csv")


cor_min$Input <- "Tmin"
cor_median$Input <- "Tmedian"
cor_orig$Input <- "Torig"

cor <- Reduce(function(x, y) full_join(x, y, 
                                       by = c("variable", "Overall", 
                                              "F. sylvatica", "P. menziesii",
                                              "Q. robur", "Input")), 
       list(cor_orig, cor_min, cor_median)) %>%
  filter(variable %in% c("Mean_VPD_kPa", "Min_VPD_kPa", "Max_VPD_kPa", 
                         "VPD_kPa", "RH", "Air_temp"))

write.csv(cor, "data/results/statistical_analysis/correlation_comparison.csv")


# compare GAMs
gam_res_min <- read_csv("data/results/statistical_analysis/GAM_results_wide.csv",
                        col_select = -c(3, 4, 6))
gam_res_median <- read_csv("data/results/statistical_analysis_VPD_median_cor/GAM_results_wide.csv",
                           col_select = -c(3, 4, 6))
gam_res_orig <- read_csv("data/results/statistical_analysis_VPD_uncor/GAM_results_wide.csv",
                         col_select = -c(3, 4, 6))


gam_res_min$Input <- "Tmin"
gam_res_median$Input <- "Tmedian"
gam_res_orig$Input <- "Torig"


gams <- Reduce(function(x, y) 
  full_join(x, y, by = c("Metric", "CWSI+SWC+VPD", "SWC+VPD", "CWSI+VPD" ,
                         "VPD", "Input")), 
              list(gam_res_orig, gam_res_min, gam_res_median))
write.csv(gams, "data/results/statistical_analysis/GAM_comparison.csv")
