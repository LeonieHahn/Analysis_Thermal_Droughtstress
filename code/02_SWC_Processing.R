library(ggplot2)
library(dplyr)
library(tidyr)
library(lubridate)
library(zoo)
library(ggrepel)
library(RColorBrewer)


file_vwc<- read.csv("./data/GH_VWC_TDR.csv")
file_vwc$Date <- as.Date(file_vwc$Date, format="%m/%d/%Y")
file_vwc$VWC <- file_vwc$VWC*100

# file with metainfos
metainfos <- read.csv("./data/Tree_Metainfos.csv")

# check if for every date values of every tree.id exist
date_u <- (unique(file_vwc$Date))

check_result <- file_vwc %>%
  filter(Date %in% date_u) %>%
  group_by(Date, Tree.ID) %>%
  summarise(Count = n()) %>%
  ungroup() %>%
  complete(Date = date_u, Tree.ID) %>%
  group_by(Date) %>%
  mutate(Complete = if_else(is.na(Count), FALSE, TRUE)) %>%
  ungroup()

check_result_f <- check_result %>% filter(Complete==FALSE)

# check missing values
count_NAs <- function(x){
  sum(is.na(x))
}

#overview of NA amount per Tree-ID
vwc_NA_count <- as.data.frame(apply(file_vwc, 2, FUN = count_NAs)) 

# remove measurements that have been taken with a different cable length than 
# 10.5
file_vwc$VWC[file_vwc$Remarks == "cable length 11"] <- NA


# remove VWC values > max. possible VWC since they do not seem to be reasonable,
# maximum vwc ist 57.03873 % (derived from extra measurements)

file_vwc$VWC[file_vwc$VWC > 57.03873 ] <- NA

vwc_NA_count2 <- as.data.frame(apply(file_vwc, 2, FUN = count_NAs)) 

# check for Missing values per Tree.ID
vwc_NA <- as.data.frame(is.na(file_vwc[4])) 
vwc_NA_2 <- cbind(file_vwc[,c(2,3,5:12)], vwc_NA)
vwc_NA_2_f <- vwc_NA_2 %>% filter(VWC=="TRUE")
vwc_NA_2_f <- arrange(vwc_NA_2_f, Tree.ID)

ggplot(vwc_NA_2_f)+geom_point(aes(x=as.Date(Date), y=Tree.ID, color=Treatment))+
  facet_grid(Tree.Species~., scales = "free_y", space = "free_y",
             switch = "y") +theme_bw() +
  scale_x_date(date_breaks = "1 week")+labs(x = "Date", 
                                            title = "Missing values ")+
  theme(strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(hjust=0.5))


#calculate mean for every pot and date
pot_means <- file_vwc %>%
  group_by(Date, Tree.ID) %>%
  summarize(VWC_mean = mean(VWC, na.rm = TRUE))

pot_means$Tree.ID_Date <- paste0(pot_means$Tree.ID, pot_means$Date)

#create dataset to merge vwc-means to 
vwc_1 <- file_vwc%>%
  filter(Measurement.Nr==1)
vwc_1$Tree.ID_Date <- paste0(vwc_1$Tree.ID, vwc_1$Date)

vwc_means <- right_join(pot_means, vwc_1, by="Tree.ID_Date")

# clean vwc_means dataset
vwc <- vwc_means[,c(1:3,9, 11:16)]
colnames(vwc)[1] <- "Date"
colnames(vwc)[2] <- "Tree_ID"

# use only data from 26.06.2023 on, 
# 26.06.2023 was last day of irrigation of all trees
vwc <- filter(vwc, as.Date(Date) >= as.Date("2023-06-26"))
vwc$Date <- as.Date(vwc$Date)
vwc$Date2 <- format(vwc$Date, "%d.%m.")

# check for Missing values per Tree.ID (Means)
vwc_NA_mean <- as.data.frame(is.na(vwc$VWC_mean)) 
vwc_NA_mean <- cbind(vwc[,-3], vwc_NA_mean)
colnames(vwc_NA_mean)[11] <- "Missing_Value"
vwc_NA_mean_f <- vwc_NA_mean %>% filter(Missing_Value=="TRUE")
vwc_NA_mean_f <- arrange(vwc_NA_mean_f, Tree_ID)

ggplot(vwc_NA_mean_f)+geom_point(aes(x=as.Date(Date), y=Tree_ID, 
                                     color=Treatment))+
  facet_grid(Tree.Species~., scales = "free_y", space = "free_y",
             switch = "y") +theme_bw() +
  scale_x_date(date_breaks = "1 week")+labs(x = "Date", 
                                      title = "Missing values (Mean vol. SWC)")+
  theme(strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(hjust=0.5))


# Linear Interpolation for values before 31.08.2023, afterwards, only Douglas
# and silver fir have been measured

vwc_interpolated <- vwc%>%
  group_by(Tree_ID) %>%
  # interpolate missing values linearly
  mutate(VWC_mean = na.approx(VWC_mean, na.rm=FALSE))
  

# check for Missing values after interpolation
vwc_interpolated_NA <- sum(is.na(vwc_interpolated$VWC_mean)) 

# check for missing values per Tree.ID (means) after interpolation
vwc_NA_mean <- as.data.frame(is.na(vwc_interpolated$VWC_mean)) 
vwc_NA_mean <- cbind(vwc_interpolated[,-3], vwc_NA_mean)
colnames(vwc_NA_mean)[11] <- "Missing_Value"
vwc_NA_mean_f <- vwc_NA_mean %>% filter(Missing_Value=="TRUE")
vwc_NA_mean_f <- arrange(vwc_NA_mean_f, Tree_ID)

ggplot(vwc_NA_mean_f)+geom_point(aes(x=as.Date(Date), y=Tree_ID, 
                                     color=Treatment))+
  facet_grid(Tree.Species~., scales = "free_y", space = "free_y",
             switch = "y") +theme_bw() +
  scale_x_date(date_breaks = "1 week")+labs(x = "Date", 
                                            title = "Missing values (Mean vol. SWC) after linear interpolation")+
  theme(strip.text.y.left = element_text(angle = 0),
        plot.title = element_text(hjust=0.5))


# filter only "measure" trees and do not further use weighing trees
vwc_cleaned <- vwc_interpolated
colnames(vwc_cleaned)[2] <- "Tree.ID"
vwc_cleaned <- left_join(vwc_cleaned, metainfos[c(2,8)], by = "Tree.ID")
vwc_cleaned <- vwc_cleaned %>% filter(Type=="Measure")

# calculate mean, min, max and standard deviation per Treatment and Species 
# per Date
stat <- vwc_cleaned %>%
  group_by(as.Date(Date), Tree.Species, Treatment) %>%
  summarise(mean_VWC = mean(VWC_mean, na.rm=TRUE),
            min_VWC = min(VWC_mean, na.rm = TRUE),
            max_VWC = max(VWC_mean, na.rm = TRUE),
            sd_VWC = sd(VWC_mean, na.rm =TRUE),
            quan05_VWC = quantile(VWC_mean, 0.05, na.rm = TRUE),
            quan95_VWC = quantile(VWC_mean, 0.95, na.rm = TRUE))
colnames(stat)[1] <- "Date"
stat$Date <- as.Date(stat$Date, format="%m/%d/%Y")
stat$min_VWC[stat$min_VWC== "Inf"] <- NA
stat$max_VWC[stat$max_VWC== "-Inf"] <- NA
stat$mean_VWC[stat$mean_VWC== "NaN"] <- NA

# check normal distribution
ggplot(vwc_cleaned, aes(x = VWC_mean)) +
  geom_histogram(binwidth=1)+
  facet_grid(Tree.Species ~ Treatment) +  
  labs(x = "VWC_mean", y = "Frequency", title = "Histogram of VWC_mean")

shapiro_test <- vwc_cleaned %>%
  group_by(Tree.Species, Treatment) %>%
  summarize(shapiro_p_value = shapiro.test(VWC_mean)$p.value)


# Lineplot per Tree.ID + Means per Treatment and Tree.Species
ggplot() +
  # geom_boxplot() +
  geom_line(data=vwc_cleaned, aes(x = as.Date(Date), y = VWC_mean, 
                                  color = Tree.ID)) +
  geom_line(data=stat, aes(x = as.Date(Date), y = mean_VWC), 
                           color = "black")+
  geom_ribbon(data = stat, aes(x = as.Date(Date), ymin = quan05_VWC, 
                               ymax = quan95_VWC, 
                               fill = Treatment), alpha = 0.3)  +
  facet_grid(Tree.Species~Treatment) +
  labs(x = "Date", y = "Vol. SWC in %", title = "Volumetric Soil Water Content") +

  ylim(0, 70)+theme_bw()+ 
  theme(
    axis.title.x = element_text(size = 14, hjust = 0.5),  
    # axis.text.x = element_blank(), 
    strip.text = element_text(size = 9),  
    plot.title = element_text(hjust = 0.5, size = 16) 
  ) 


# Pointplot with Means and 5th and 95th quantils per Treatment and Tree.Species
ggplot() +
  geom_ribbon(data = stat, aes(x = as.Date(Date), ymin = quan05_VWC, 
                               ymax = quan95_VWC, 
                               fill = Treatment), alpha = 0.2)  +
  geom_point(data=vwc_cleaned, aes(x = as.Date(Date), y = VWC_mean, 
                                  color = Treatment, shape=Treatment)) +
  geom_line(data=stat, aes(x = as.Date(Date), y = mean_VWC, 
            color = Treatment), linewidth=1.2)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  facet_wrap(.~Tree.Species) +
  labs(x = "Date", y = "Vol. SWC in %", title = "Volumetric Soil Water Content") +
  ylim(0, 70)+theme_bw()+ 
  theme(
    axis.title = element_text(size = 16, hjust = 0.5),  
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),  
    plot.title = element_text(hjust = 0.5, size = 16)
  )


# remove obvious measurement errors
vwc_cleaned <- vwc_cleaned %>%
  mutate_at(vars(VWC_mean), ~ case_when(
    Tree.ID == "B27" & Date == as.Date("2023-07-07") |
      Tree.ID == "E45" & Date == as.Date("2023-07-20") ~ NA,
    TRUE ~ .  # Keep all other cases
  ))


# ad rel. extractable VWC (REW), see also 
# Steger, D.N., Peters, R.L., Blume, T., Hurley, A.G., Balanzategui, D., 
# Balting, D.F., Heinrich, I., 2024. Site matters - Canopy conductance 
# regulation in mature temperate trees diverges at two sites with contrasting 
# soil water availability. Agricultural and Forest Meteorology 345, 109850. 
# https://doi.org/10.1016/j.agrformet.2023.109850

vwc_cleaned$REW <- (vwc_cleaned$VWC_mean/57.03873)


# interpolate vol. SWC linearly to a daily basis

# create full timeframe
full_dates <- seq.Date(min(vwc_cleaned$Date), max(vwc_cleaned$Date), by="day")
full_tree_ids <- unique(vwc_cleaned$Tree.ID)
full_combinations <- expand.grid(Date = full_dates, Tree.ID = full_tree_ids)

# combine data
vwc_full <- full_combinations %>%
  left_join(vwc_cleaned, by = c("Date", "Tree.ID"))

# treatment and tree species
vwc_full <- vwc_full %>%
  group_by(Tree.ID) %>%
  fill(Treatment, Tree.Species, .direction = "downup") %>%
  ungroup()


# Interpolate missing values linearly, so we have values for every day
vwc_interpolated2 <- vwc_full %>%
  group_by(Tree.ID) %>%
  arrange(Date) %>%
  mutate(VWC_mean = zoo::na.approx(VWC_mean, na.rm = FALSE),
         REW = zoo::na.approx(REW, na.rm = FALSE)) %>%
  ungroup()

write.csv(vwc_interpolated2, 
          "./data/processed/GH_VWC_TDR_interpolated.csv")

ggplot(vwc_interpolated2)+
  geom_point(aes(x=Date, y= VWC_mean, color=Treatment))+
  facet_wrap(.~Tree.Species)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  facet_wrap(.~Tree.Species) +
  labs(x = "Date", y = "Vol. SWC in %", 
       title = "Daily Volumetric Soil Water Content") +
  ylim(0, 70)+theme_bw()+ 
  theme(
    axis.title = element_text(size = 16, hjust = 0.5),  
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),  
    plot.title = element_text(hjust = 0.5, size = 16))


# Plot only vol. SWC for dendrometer trees
dendro_info <- read.csv("~/Bewaesserung_Forstkulturen/GH/GH_Analysis/Analysis_dendrometer/data/Tree_Metainfos_dendro.csv")
dendro_trees <- dendro_info[,2]

vwc_dendro_trees <- left_join(dendro_info[,c(2, 5, 6)], vwc_interpolated2,
                              by=c("Tree.ID", "Treatment", "Tree.Species"))


ggplot(vwc_dendro_trees)+
  geom_point(aes(x=Date, y= VWC_mean, color=Treatment))+
  facet_wrap(.~Tree.Species)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  facet_wrap(.~Tree.Species) +
  labs(x = "Date", y = "Vol. SWC in %", 
       title = "Daily volumetric soil water content of dendrometer trees") +
  ylim(0, 70)+theme_bw()+ 
  theme(
    axis.title = element_text(size = 16, hjust = 0.5),  
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),  
    plot.title = element_text(hjust = 0.5, size = 16))



# Plot with mean and standard deviation per Treatment and Tree.Species
dendrostat <- vwc_dendro_trees %>%
  group_by(as.Date(Date), Tree.Species, Treatment) %>%
  summarise(mean_VWC = mean(VWC_mean, na.rm=TRUE),
            min_VWC = min(VWC_mean, na.rm = TRUE),
            max_VWC = max(VWC_mean, na.rm = TRUE),
            sd_VWC = sd(VWC_mean, na.rm =TRUE),
            quan05_VWC = quantile(VWC_mean, 0.05, na.rm = TRUE),
            quan95_VWC = quantile(VWC_mean, 0.95, na.rm = TRUE))
colnames(dendrostat)[1] <- "Date"
dendrostat$Date <- as.Date(dendrostat$Date, format="%m/%d/%Y")
dendrostat$min_VWC[dendrostat$min_VWC== "Inf"] <- NA
dendrostat$max_VWC[dendrostat$max_VWC== "-Inf"] <- NA
dendrostat$mean_VWC[dendrostat$mean_VWC== "NaN"] <- NA

plot_dendrovwc <- ggplot() +
  geom_ribbon(data = dendrostat, aes(x = as.Date(Date), ymin = mean_VWC-sd_VWC, 
                               ymax = mean_VWC+sd_VWC, 
                               fill = Treatment), alpha = 0.2)  +
  geom_line(data=dendrostat, aes(x = as.Date(Date), y = mean_VWC, 
                           color = Treatment), linewidth=1.2)+
  scale_color_manual(values=c("red", "blue"))+
  scale_fill_manual(values=c("red", "blue"))+
  facet_wrap(.~Tree.Species) +
  labs(x = "Date", y = "Vol. SWC [%]", 
       title = "Mean and standard deviation of vol. soil water content (SWC)
       of trees with dendrometer") +
  ylim(0, 70)+theme_bw()+ 
  theme(
    axis.title = element_text(size = 16, hjust = 0.5),  
    axis.text = element_text(size = 12),
    strip.text = element_text(size = 12),  
    plot.title = element_text(hjust = 0.5, size = 16)
  )

print(plot_dendrovwc) 
